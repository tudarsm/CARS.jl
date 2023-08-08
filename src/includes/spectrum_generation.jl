# per default, assume Float64 as datatype
"""
    generateSpectrum!(params::FitParams,libs...;takesqrt=CARSopt.fittype)

Convenience function to simulate spectra outside fitting mode.
"""
function generateSpectrum!(params::FitParams,libs...;takesqrt=CARSopt.fittype)
    loc = [Dict("fn"=>:nothing)]
    return generateSpectrum!(params,1.,loc,libs...;takesqrt=takesqrt)
end

"""
    generateSpectrum!(params::FitParams,a::T,loc,libs...;takesqrt=CARSopt.fittype) where {T<:Real}

Generates a spectrum based on the parameters provided by ```params``` using the libraries ```libs...```. The loc array is required to identify if some of the degrees of freedom are required, e.g. frequency offset, where the interpolation is not performed if not required.

!!! note
    All methods wihtin this function have to be compatible with ```T<:Real``` to be compatible with FowardDiff.Duals!
"""
function generateSpectrum!(params::FitParams,a::T,loc,libs...;takesqrt=CARSopt.fittype) where {T<:Real}
    ############################################
    # Library parameters
    ############################################
    γ = convert(T,libs[1].γc)
    ω = libs[1].ω
    ωres = convert(T,libs[1].ωres)   # wavenumber resolution used in library
    num_pr = length(libs)   # number of pumping regions
    
    ############################################
    # Preallocations optimized with memoizing
    ############################################
    # Attention: this makes this function not thread-safe!
    χM,χR,χI,I = memoized_preallocate(convert(T,1),num_pr,length(ω))
    
    ############################################
    # Calculate the non-resonant signal
    ############################################
    P = libs[1].P     # pressure from the library
    N = NA*273.15/params.T*(P/molvol);          # number density based on temperature and pressure
    addχNRtoParams!(params,libs...) # does not have to be done every iteration.
    χB = params.chiB/NA*molvol;       # buffer gas susc. per molecule
    XB = 1-sum(values(params.X))    # buffer gas mole fraction
    χNR = convert(T,N*(XB*χB + sum([params.X[key]*params.chiNR[key] for (key,value) in params.X])))  # entire non resonant signal

    ############################################
    # library access and temperature interpolation
    ############################################
    # get the temperature index and values for interpolation
    idx0 = findlower(libs[1].T,params.T)
    idx1 = idx0+1
    idxs = [idx0 idx1]
    T0 = libs[1].T[idxs[1]]
    T1 = libs[1].T[idxs[2]]
    Tx = params.T
    
    # weights for interpolation
    w0 = (T1-Tx)/(T1-T0)
    w1 = (Tx-T0)/(T1-T0)
    #print("Library Access:")
    for pr in 1:size(χM,2)
        si = 1
        for species in libs[pr].species
            
            X = convert(T,params.X[species]);
            # short handle for better readability: these are slices of the library between which to interpolate
            χMv = @view libs[pr].χM[:,idx0:idx1,si]
            χRv = @view libs[pr].χR[:,idx0:idx1,si]
            χIv = @view libs[pr].χI[:,idx0:idx1,si]
            
            # interpolate library between temperatures
            interpSpec!(χM,w0,w1,χMv,X^2,pr)
            interpSpec!(χR,w0,w1,χRv,X,pr)
            interpSpec!(χI,w0,w1,χIv,X,pr)
            
            # species index
            si +=1
        
        end
        # do not forget the cross-terms, if there are any.
        # later: make this optional. it really is only necessary for significant overlap
        # of species inside a given pumping region.
        if si>2     # si will be incremented at least once at the end of the upper loop. si>2 in fact means more than one species
            for i in 1:length(libs[pr].species), k in 1:i-1
                Xki = convert(T,params.X[libs[pr].species[i]]* params.X[libs[pr].species[k]])
                if Xki != 0 # only do this interpolation if Xki is != 0
                    χCv = @view libs[pr].χC[k,i][:,idx0:idx1]
                    interpSpec!(χM,w0,w1,χCv,Xki,pr)
                end
            end
        end

        ############################################
        # convolution to final laser linewidths
        ############################################
        # convolve each pumping region with the "leftover" gaussian kernel
        pr == 1 ? PR = "P1" : PR = "P2"
        # only convolve if necessary
        # simplified libraries have the linewidth baked in, so skip this step
        # check against linewidth from library - if it is the same, do not convolve
        # need to use isapprox() because of potential conversion errors
        # between Float32 and Float64
        γi = convert(T,libs[pr].γc)
        if params.LineWidth[PR] != 0 && !isapprox(params.LineWidth[PR]^2,γi^2;atol=1e-5)
            γk = params.LineWidth[PR]^2-γ^2
            if γk < 0
                error("Linewidth too small for this library.")
            end
            kern = gaussiankernelT(ωres,sqrt(γk));
            χM[:,pr]=dualConv(χM[:,pr],kern)
            χI[:,pr]=dualConv(χI[:,pr],kern)
            χR[:,pr]=dualConv(χR[:,pr],kern)
        end
    end

    ############################################
    # AdditionalFrequencyOffset
    ############################################
    # if AdditionalFrequencyOffset != AdditionalFrequencyOffset(Library),
    # interpolate the second pumping region on a shifted grid
    # check if frequency offset is one of the fitting parameters
    # or if it is non-zero, but fixed.
    if size(χM,2) == 2  # only do this for Dual Pump spectra
        LibraryFrequencyOffset = convert(T,libs[2].AdditionalFrequencyOffset)
        if :AdditionalFrequencyOffset in [loc[i]["fn"] for i in 1:length(loc)] || params.AdditionalFrequencyOffset != LibraryFrequencyOffset
            # χM[:,2] = interp1(ω .+ params.AdditionalFrequencyOffset,ωres,χM[:,2],ω)
            # χI[:,2] = interp1(ω .+ params.AdditionalFrequencyOffset,ωres,χI[:,2],ω)
            # χR[:,2] = interp1(ω .+ params.AdditionalFrequencyOffset,ωres,χR[:,2],ω)
            # Faster implementation, but requires Float64 for the destination array
            ω64=convert(Vector{Float64},ω)
            χM[:,2] = fastLinearInt(ω .+ params.AdditionalFrequencyOffset,χM[:,2],ω64;extrap=true)
            χI[:,2] = fastLinearInt(ω .+ params.AdditionalFrequencyOffset,χI[:,2],ω64;extrap=true)
            χR[:,2] = fastLinearInt(ω .+ params.AdditionalFrequencyOffset,χR[:,2],ω64;extrap=true)
        end
    end

    ############################################
    # Sum the signal
    ############################################
    # output is stored in I
    signalsummation!(I,χM,χI,χR,χNR)

    ############################################
    # Convolution with Gaussian instrumental
    ############################################
    if params.Instrumental > 0
        I = dualConv(I,gaussiankernelT(ωres,convert(T,params.Instrumental)))
    end


    ############################################
    # take square root (if desired)
    ############################################
    #@. I = sqrt(I)
    if takesqrt==:sqrt && CARSopt.fittype == :sqrt
        if minimum(I)<0
            # sometimes Ie can be slightly negative in parts, where the NR signal and the real part interfere
            # this is caused by splitting the convolution in two steps, introducing inaccuracies
            # return a warning if this get's really high
            # the threshold is: the absolute value of the minimum has 0.1% of the peak intensity.
            # everything above that cannot be explained with this convolution artefact
            if abs(minimum(I)/maximum(I)) > 1e-3
                warn("Spectrum generation returned a negative number. If this warning occurs often, something went wrong.")
                show(params)
            end
            @. I[I<0]=0   # as a workaround, set these negative values to zero.
            #Ie = Ie .- minimum(Ie) # this would be an alternative: shift the entire spectrum up. but this may cause issues with the mole fraction result
        end    
        @. I = sqrt(I)
    end


    return ω,I,ωres

end
