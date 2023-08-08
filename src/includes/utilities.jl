function warn(msg)
printstyled("$(msg)\n";color=:yellow)
end

"""
    isFitParam(loc,name)
Returns `true` if `name` is a fit parameter by detecting if it is in the `loc` array created by [`vectorizeParams!`](@ref)
"""
function isFitParam(loc,name)
   return name in [loc[i]["fn"] for i in eachindex(loc)]
end


function getWavenumberPolynomial(IE,x...)
    # because data is collected in wavelength domain and we want wavenumber domain,
    # we need to invert the pixel coordinates
    x1 = length(IE)-x[1][1]
    x2 = length(IE)-x[2][1]

    m = (x[2][2]-x[1][2])/(x2-x1)
    b = x[1][2] - m*x1
    return [b,m,0]
end


"""
vectorizeParams!
"""
function vectorizeParams!(startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Tuple)

    
    # vectorize the parameters and startingsolutions used in fitting
    x0=[]
    xlower=[]
    xupper=[]
    loc=[]
    # variables to be used in fit are those that are different between lower and upper
    # get those:
    fns = fieldnames(FitParams)
    for fn in fns
        # check if it is different. if so, add starting solution
        if getproperty(lower,fn) != getproperty(upper,fn)
            # Some fields should be ignored.
            if fn == :Spectra
                continue
            end
            # the property might be a scalar, an array or a Dict
            # for the latter two, iterate through the fields/indices
            if isa(getproperty(lower,fn),Dict)
                for key in keys(getproperty(lower,fn))
                    if getproperty(lower,fn)[key] != getproperty(upper,fn)[key]
                        push!(x0,getproperty(startingsolution,fn)[key])
                        push!(xlower,getproperty(lower,fn)[key])
                        push!(xupper,getproperty(upper,fn)[key])
                        push!(loc,Dict("fn"=>fn,"idx"=>key))
                        #println("Adding $(fn)[$(key)] as fitting variable.")
                    end
                end
            elseif isa(getproperty(lower,fn), Array)
                for idx in eachindex(getproperty(lower,fn))
                    if getproperty(lower,fn)[idx] != getproperty(upper,fn)[idx]
                        push!(x0,getproperty(startingsolution,fn)[idx])
                        push!(xlower,getproperty(lower,fn)[idx])
                        push!(xupper,getproperty(upper,fn)[idx])
                        push!(loc,Dict("fn"=>fn,"idx"=>idx))
                        #println("Adding $(fn)[$(idx)] as fitting variable.")
                    end
                end
            else
                push!(x0,getproperty(startingsolution,fn))
                push!(xlower,getproperty(lower,fn))
                push!(xupper,getproperty(upper,fn))
                push!(loc,Dict("fn"=>fn))
                #println("Adding $(fn) as fitting variable.")

            end
            
        end
    end


    return x0,xlower,xupper,loc

end

function devectorizeParams!(params::FitParams,x,loc)

    # DeVectorization of parameters
    # takes an array of values x and puts it into the given params struct based on the locations loc provided by vectorizeParams

    for idx in eachindex(loc)
        
        # check if the current location is an array or dict, in this case, add it do the correct location
        if isa(getproperty(params,loc[idx]["fn"]),Dict) || isa(getproperty(params,loc[idx]["fn"]),Array)
            a = getproperty(params,loc[idx]["fn"])
            a[loc[idx]["idx"]]=x[idx]
            setproperty!(params,loc[idx]["fn"],a)
        else
            setproperty!(params,loc[idx]["fn"],x[idx])
        end
    end
end


function findnearest(Ts,t)
    return Int(findmin(abs.(Ts.-t))[2])
end

function findlower(Ts,t)
    idxs = Ts.<=t
    x = 0
    for idx in idxs
        x+=1
        if idx == 0
            break
        end
    end
    return Int(x-1)
end

function csvread(filename::String)
    s = open(filename) do file
        readlines(file)
    end

    ω = zeros(length(s))
    I = similar(ω)
    for i = 1:length(s)
        l = split(s[i],",")
        ω[i]=parse(Float64,l[1])
        I[i]=parse(Float64,l[2])
    end
    return ω,I
end

function permolecule(χNR)
    return χNR/NA*molvol
end

function N(T,P)
    return NA*273.15/T*(P/molvol);
end

function addχNRtoParams!(params,libs...)
    # get list of available species
    species = Vector{String}(undef,0)
    χNRs = Vector{Float32}(undef,0)
    for lib in libs, spec in lib.species
        push!(species,spec)
    end
    for lib in libs, χNR in lib.χNR
        push!(χNRs,χNR)
    end
    # make a dictionary from these vectors
    params.chiNR = Dict{String,Float32}( (species[i] => χNRs[i] for i=1:length(species)) )

end

"""
This function can be used to conveniently store data in either a .MAT or .JLD2 file, depending on what you would like to use for further processing
    If the file exists, the data is added to the file.
"""
function save(filename,varname,data)
if contains(lowercase(filename),".mat")
    return saveMat(filename,varname,data)
else
    error("Currently, only .mat files are supported.")
end
end

function saveMat(filename,varname,data)
    # default: append to file, see https://github.com/JuliaIO/MAT.jl/blob/master/src/MAT.jl
    file = matopen(filename, true , true , true , false, true, true)
    try
        write(file,varname,data)
    catch e
        warn(e)
    end

    close(file)
end