"""
    generateLibrary(T,γc,getSuscs...)
Generates libraries in the temperature range `T`, e.g. `T=280:20:2000` using the compression kernel `γc` and the supplied methods `getSuscs...`
"""
function generateLibrary(T,γc,getSuscs...)
    P=1 # add from getSuscs. Needs change in DiaCARS and CO2CARS
    # convert some types:
    γc = convert(Float32, γc);

    # for preallocating some of the arrays, do a dummy simulation for all included species here
    # also, use this to check for potential errors thrown by getSusc
    Ω=[];
    MolecularParameters=[];
    Species = Array{String}(undef,0);
    for getSusc in getSuscs
        try
            χR,χI,ω,species,molecularParameters = getSusc(maximum(T));
            # store relevant information for later checks in arrays:
            push!(Species,species)
            push!(Ω,ω)
            push!(MolecularParameters,molecularParameters)
        catch e
            # if an error occurs, abort the function and return nothing
            println("Error returned by $(getSusc):")
            warn(e)
            return nothing
        end
    end


    # check if all omegas have the same size
    sz = size(Ω[1])
    for i in eachindex(Ω)
        if size(Ω[i]) !=sz 
            error("Inconsistent size of wavenumber arrays. Found $(size(Ω[i])) instead of $(sz)")
        end
    end

    # the new wavenumberarray should be the same ROI, but coarser
    # ROI cannot be used here directly, instead it should be ω from the calculation
    # because this includes the ROI expansion required to avoid convolution edge effects
    # from experience, 1/5 of the preconvolution size is sufficient
    # 5 is the default value for CARSopt.γf. to override it, set a new default value in your script.
    ωcres = γc/CARSopt.γf
    ω = convert(Vector{Float32},Ω[1])
    ωc = convert(Vector{Float32}, collect(ω[1]:ωcres:ω[end]))
    ωres = convert(Float32,round(sum(diff(ω))/length(ω);sigdigits=5))    # get the original resolution

    # construct the gaussian compression kernel
    gc = gaussiankernelT(ωres,γc)
    
    # print some information
    println("Generating a library for species $(join(Species, ',')) in range $(ω[1])-$(ω[end]) cm-1")
    println("Array size before compression: $(length(ω)) (Original resolution: $(ωres) cm-1)")
    println("Array size after compression: $(length(ωc)) (New resolution: $(ωcres) cm-1)")

    # number of cross-terms
    N = length(getSuscs)
    Nc = Int((N^2-N)/2)
    println("Number of required complex cross-terms: $(Nc)")

    # preallocate arrays for real and imaginary components
    χR = Array{Float32}(undef,length(ωc),length(T),N)
    χI = Array{Float32}(undef,length(ωc),length(T),N)
    χM = Array{Float32}(undef,length(ωc),length(T),N)

    # preallocate cross terms if necessary
    χRc,χIc = 0,0
    χRci, χIci = 0,0
    χC = Matrix{Array{Float32}}(undef,N,N)   # convenient way of storing cross terms in a square matrix with just the upper triangle is populated

    if Nc > 0 
        # this is where the cross-terms will be stored
        χRc = Array{Float32}(undef,length(ωc),length(T),Nc)
        χIc = Array{Float32}(undef,length(ωc),length(T),Nc)

        # preallocate scratch arrays for intermediate steps
        Rc = Array{Float32}(undef,length(ω))
        Ic = Array{Float32}(undef,length(ω))
        
        # these are temporary arrays for computing the cross-terms at each temperature
        χRci = Array{Float32}(undef,length(ω),N)
        χIci = Array{Float32}(undef,length(ω),N)
    end

    # estimate library dimension:
    libsize = sizeof(χR)+sizeof(χI)+sizeof(χM)+sizeof(χRc)

    println("Estimated library size: $(round(libsize/1024/1024)) MB")
    println("")

    # now  loop over all temperatures
    ti = 1
    for t in ProgressBar(T)
        # and loop over all species
        sj=1
        for getSusc in getSuscs
            # getSusc is the function to simulate theoretical susceptibility spectra
            # simulate the spectra and calculate the squared sum
            χRl,χIl = getSusc(t);
            χMl = @. χRl^2+χIl^2
            
            # if there are cross-terms, store the uncompressed real and imaginary component for this iteration
            if Nc > 0 
                χRci[:,sj] = χRl
                χIci[:,sj] = χIl
            end

            # now compress the components
            # first, convolve them:
            χRl = conv(χRl,gc)
            χIl = conv(χIl,gc)
            χMl = conv(χMl,gc)

            # now reinterpolate and store in library:
            χR[:,ti,sj] = interp1(ω,ωres,χRl,ωc;extrapolation=false);
            χI[:,ti,sj] = interp1(ω,ωres,χIl,ωc;extrapolation=false);
            χM[:,ti,sj] = interp1(ω,ωres,χMl,ωc;extrapolation=false);
            
            # increase inner counter
            sj+=1
        end

        # compute species cross-terms for this temperature level
        if Nc > 0
            # calculate, convolve, reinterpolated
            Nci=1
            for i in 1:N, k in 1:i-1
                # real cross terms
                @. Rc = χRci[:,i]*χRci[:,k]
                # imaginary cross-terms
                @. Ic = χIci[:,i]*χIci[:,k]
                # convolve, reinterpolate and store in library
                Rc = conv(Rc,gc)
                Ic = conv(Ic,gc)
                χRc[:,ti,Nci] = interp1(ω,ωres,Rc,ωc;extrapolation=false);
                χIc[:,ti,Nci] = interp1(ω,ωres,Ic,ωc;extrapolation=false);
                Nci +=1
            end
        end
        
        # increase outer counter
        ti+=1
    end

    # store the cross-terms in a convenient way:
    if Nc > 0
        # These are NxN matrices where only the upper triangle
        # is filled with the cross-terms      
        Nci=1
        for i in 1:N, k in 1:i-1
            χC[k,i] = 2*χRc[:,:,Nci]+2*χIc[:,:,Nci]
            Nci +=1
        end
    end

    # collect all further relevant information and construct a dictionary for the output
    χNR = [species.CHINR for species in MolecularParameters]

    library = Library(
        Species,
        vec([ωc[1] ωc[end]]),
        0.,  # by default, set a frequency offset of zero. at this point, we do not know it.
        collect(T),
        P,
        ωc,
        ωcres,
        γc,
        χC,
        χM,
        χI,
        χR,
        χNR # convert here to float32 for later spectrum generation
        );


        # do some sanity checks
        if sum(isnan.(library.χM)) != 0
            error("Library spectra contains NaNs. Try to increase resolution (lower ωres).")
        end

    return library
end


"""
This function uses the linewidth and additionalfrequencyoffset in fitparams and applies them to the library.
    This allows to treat these as fixed parameters for the fit and make it much(!) faster
"""
function simplifyLibrary!(fitparams::FitParams,libs::Library...)

    # get the parameters from the current library
    γ = libs[1].γc
    ω = libs[1].ω
    ωres = libs[1].ωres   # wavenumber resolution used in library
    num_pr = length(libs)   # number of pumping regions
    T = eltype(libs[1].χM)

    # iterate through the libraries
    for i in 1:num_pr
        lib = libs[i]
        # calculate the kernel for convolution
        i == 1 ? PR = "P1" : PR = "P2"
        # make sure that it is possible to use
        if fitparams.LineWidth[PR]^2-γ^2 <= 0
            error("Linewidth for compression smaller than library linewidth!")
        end

        # make sure the kernel has the same type as the original susceptibility arrays
        γn = convert(typeof(ωres),sqrt(fitparams.LineWidth[PR]^2-γ^2))  # additional gaussian for convolution
        kern = convert(Vector{T},gaussiankernelT(ωres,γn))
        
        # iterate through species
        si = 1
        for species in lib.species
            # update the laser linewidth entry
            lib.γc = fitparams.LineWidth[PR] # attention: this has to be the absolute linewidth, not the additional one!
            # do the convolutions and interpolations
            if i == 1
                println("Convolving $(species) in pump region $(i)")
            else
                lib.AdditionalFrequencyOffset = fitparams.AdditionalFrequencyOffset
                println("Convolving and adding frequency offset to $(species) in pump region $(i)")
            end
            for j in ProgressBar(1:length(lib.T))
                lib.χM[:,j,si]=lib.χM[:,j,si]∗kern
                lib.χI[:,j,si]=lib.χI[:,j,si]∗kern
                lib.χR[:,j,si]=lib.χR[:,j,si]∗kern
                
                # additionally, interpolate to a new grid!
                if i == 2
                    lib.χM[:,j,si] = interp1(ω .+ fitparams.AdditionalFrequencyOffset,ωres,lib.χM[:,j,si],ω)
                    lib.χI[:,j,si] = interp1(ω .+ fitparams.AdditionalFrequencyOffset,ωres,lib.χI[:,j,si],ω)
                    lib.χR[:,j,si] = interp1(ω .+ fitparams.AdditionalFrequencyOffset,ωres,lib.χR[:,j,si],ω)
                end
            end
        
            # species iterator
            si += 1
        end
        # do not forget cross terms, if there are any
        if length(lib.species)>1
            println("Convolving and adding frequency offset to cross-terms in pump region $(i)")
            # convolve the cross-terms
            for r in 1:length(lib.species), k in 1:r-1 # r instead of i, because i is used in the outer loop
                for j in ProgressBar(1:length(lib.T)), si in 1:size(lib.χC[k,r],3)
                    lib.χC[k,r][:,j,si]=lib.χC[k,r][:,j,si] ∗ kern
                    # if this is in pumping region 2, also do the frequency offset here
                    lib.χC[k,r][:,j,si] = interp1(ω .+ fitparams.AdditionalFrequencyOffset,ωres,lib.χC[k,r][:,j,si],ω)
                end
            end
        end
    end

    # ok, now we have the convolutions and frequency offsets handled.
    # reinterpolate everything on a coarser grid for higher efficiency
    # for this, figure out what the smallest linewidth is:
    γs = minimum([fitparams.LineWidth["P1"] fitparams.LineWidth["P2"]])
    ωresn = convert(Float32,round(γs/CARSopt.γf;sigdigits=1)) # limit precision to avoid OBOE in generation of wavenumber range
    
    # get the destination array
    ωn = collect(ω[1]:ωresn:ω[end])
    # preallocate the new arrays
    # use a for loop to reinterpolate everything
    # replace in library
    # the following arrays need reinterpolation:
    # χM -> ω,T,species
    # χR -> ω,T,species
    # χI -> ω,T,species
    # χC -> species,species -> ω,T
    println("Reinterpolating all spectral components on a coarser grid...")
    for lib in ProgressBar(libs)
        # set the new library resolution and grid:
        lib.ωres = ωresn
        lib.ω = ωn

        lib.χM = interpolateLibrarySpectrum(ω,ωres,ωn,lib.χM)
        lib.χR = interpolateLibrarySpectrum(ω,ωres,ωn,lib.χR)
        lib.χI = interpolateLibrarySpectrum(ω,ωres,ωn,lib.χI) 
        
        # don't forget the crossterms
        for r in 1:length(lib.species), k in 1:r-1 # r instead of i, because i is used in the outer loop
            lib.χC[k,r] = interpolateLibrarySpectrum(ω,ωres,ωn,lib.χC[k,r])
        end
    end

    # Crop the library closer to the experiment
    ωexp = fitparams.Spectra.omega
    # extract the minima and maxima
    ωmin = minimum(ωexp)
    ωmax = maximum(ωexp)
    # these should always be within the bounds of the library.
    # but better to test anyway:
    if ωmin < minimum(ωn) || ωmax > maximum(ωn)
        error("Cannot simplify Library. Wavenumber array from experiment is not within the library bounds!")
    end
    # We want some margin to avoid edge effects when convolving with the apparatus function
    margin = CARSopt.conv_margin*fitparams.Instrumental
    # find the index of the lowest wavenumber array that we need
    idxmin = findfirst(ωn.>(ωmin-margin))
    idxmax = findlast(ωn.<(ωmax+margin))
    # now crop the libraries
    println("Cropping library closer to experimental bounds $(idxmin):$(idxmax) instead of 1:$(length(ωn))...")

    for lib in libs
        lib.χM = lib.χM[idxmin:idxmax,:,:]
        lib.χR = lib.χR[idxmin:idxmax,:,:]
        lib.χI = lib.χI[idxmin:idxmax,:,:]
        lib.ω = lib.ω[idxmin:idxmax]

        # don't forget the crossterms
        for r in 1:length(lib.species), k in 1:r-1 # r instead of i, because i is used in the outer loop
            lib.χC[k,r] = lib.χC[k,r][idxmin:idxmax,:,:]
        end

    end

    
    println("Done.")
end

function interpolateLibrarySpectrum(ω,ωres,ωn,χ)
    # preallocate return array
    χn = zeros(eltype(χ),(length(ωn),size(χ)[2:end]...))
    # iterate over all but the first dimension
    # it is 3 dims tops, so hardcoding is fine
    for i = 1:size(χn,2), j = 1:size(χn,3)
        χn[:,i,j] = interp1(ω,ωres,χ[:,i,j],ωn;extrapolation=false)
    end
    return χn
end