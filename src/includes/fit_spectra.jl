# use blackbox to get a good starting solution and then fit with gradientfitter
"""
combinedFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;blackboxopt::BlackBoxOptions=BlackBoxOptions(MaxFuncEvals=200,TraceMode=:silent),ipoptopt::IpoptOptions=IpoptOptions(print_level=0),showplots=false)
Use [`blackBoxFit`](@ref) first which provides a starting solution for [`gradientFit`](@ref) for final fitting.
"""
function combinedFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;blackboxopt::BlackBoxOptions=BlackBoxOptions(MaxFuncEvals=200,TraceMode=:silent),ipoptopt::IpoptOptions=IpoptOptions(print_level=0),showplots=false)
    params = blackBoxFit(IE,startingsolution,lower,upper,libs...;options=blackboxopt,showplots=false) # show plots only for final fit
    params = gradientFit(IE,params,lower,upper,libs...;options=ipoptopt,showplots=showplots)
    return params
end

function combinedFit(IE::AbstractArray{PreProcSpectrum},moreargs...;blackboxopt::BlackBoxOptions=BlackBoxOptions(MaxFuncEvals=200,TraceMode=:silent),ipoptopt::IpoptOptions=IpoptOptions(print_level=0),showplots=false)
    paramsArray = []
    iter = ProgressBar(1:length(IE))
    for j in iter
        set_description(iter, "Fitting Spectrum $(j) of $(size(IE,2))")
        params = combinedFit(IE[j],moreargs...;blackboxopt=blackboxopt,ipoptopt=ipoptopt,showplots=showplots)
        push!(paramsArray,params)
    end
    
    return paramsArray
end

# fit multiple spectra in one function if IE is supplied as a matrix
function blackBoxFit(IE::AbstractArray{PreProcSpectrum},moreargs...;options::BlackBoxOptions=BlackBoxOptions(),showplots=false)
    paramsArray = []
    iter = ProgressBar(1:length(IE))
    for j in iter
        set_description(iter, "Fitting Spectrum $(j) of $(length(IE))")
        params = blackBoxFit(IE[j],moreargs...;options=options,showplots=showplots)
        push!(paramsArray,params)
    end
    
    return paramsArray
end

function gradientFit(IE::AbstractArray{PreProcSpectrum},moreargs...;options::IpoptOptions=IpoptOptions(),showplots=false)
    paramsArray = []
    iter = ProgressBar(1:length(IE))
    for j in iter
        set_description(iter, "Fitting Spectrum $(j) of $(length(IE))")
        params = gradientFit(IE[j],moreargs...;options=options,showplots=showplots)
        push!(paramsArray,params)
    end
    return paramsArray
end

"""
    blackBoxFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;options::BlackBoxOptions=BlackBoxOptions(),showplots=false)
Use BlackBoxOptim.jl for fitting.
"""
function blackBoxFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;options::BlackBoxOptions=BlackBoxOptions(),showplots=false)
        
    ############################################
    # Check whether the sample should be evaluated or not based on errorcode
    ############################################
    if IE.ErrorCode != 0
        # for these cases, create a FitParams struct containing NaNs for the relevant parts
        # include the experimental spectra, but set fit to zero
        invalidParams = deepcopy(startingsolution)
        invalidParams.T = NaN
        for (key,value) in invalidParams.X
            invalidParams.X[key] = NaN
        end
        invalidParams.Spectra = SpectraStruct(IE.I./sum(IE.I),zeros(size(IE.I)),zeros(size(IE.I)),getWavenumberArray(startingsolution.WavenumberPolynomial,length(IE.I)),IE.ErrorCode,"Ipopt")
        if showplots
            plotSpec(invalidParams)
        end
        return invalidParams
    end

    global Ie
    ############################################
    # make a deep copy of startingsolution to not overwrite it
    ############################################
    params=deepcopy(startingsolution)
    # pass targetfitness as global variable
    TargetFitness = copy(options.TargetFitness)
    
    # extract fitting variables and vectorize them
    x0,xlower,xupper,loc=vectorizeParams!(startingsolution,lower,upper,libs)
    # x0 has to be converted to Float64 for BBOptim to work...
    x0=convert(Vector{Float64},x0)
    xlower=convert(Vector{Float64},xlower)
    xupper=convert(Vector{Float64},xupper)

    ############################################
    # pass library and experiment to residual function using closure
    ############################################
    evoluationaryResidual(x) = calcResidual(IE,libs,params,loc,x...)

    ################## BBOPTIM ##########################
    # get a progressmeter to update in callback
    Prog = Progress(options.MaxFuncEvals; showspeed=true) 
    global Prog, MaxFuncEvals, TargetFitness
    if options.TraceMode != :silent
        show(options)
    end
    result = bboptimize(evoluationaryResidual,x0;
                    SearchRange = [(xlower[i],xupper[i]) for i in 1:length(xlower)],
                    NumDimensions = length(x0),
                    TraceMode=options.TraceMode,
                    MaxFuncEvals=options.MaxFuncEvals,
                    Method = options.Method,
                    PopulationSize=options.PopulationSize,
                    CallbackFunction = bbcallback,  # this callback also implements TargetFitness, because this does not seem wo work directly, see https://github.com/robertfeldt/BlackBoxOptim.jl/issues/209
                    CallbackInterval = 0.5)

    ################## BBOPTIM ##########################
    # store the spectra (experiment, fit, residual and wavenumbervector) in the params struct
    quality = 0
    params.Spectra = SpectraStruct(IE.I./sum(IE.I),Ie,IE.I./sum(IE.I)-Ie,getWavenumberArray(params.WavenumberPolynomial,length(IE.I)),quality,"BBOptim")

    if showplots
        plotSpec(params)
    end

    return params
end

function bbcallback(oc)
    global Prog, MaxFuncEvals, TargetFitness
    ProgressMeter.update!(Prog, BlackBoxOptim.num_func_evals(oc); showvalues = [(:Iteration,BlackBoxOptim.num_func_evals(oc)), (:Fitness,best_fitness(oc))])
    if best_fitness(oc)<TargetFitness
        ProgressMeter.finish!(Prog)
        BlackBoxOptim.shutdown!(oc)
    end
end

"""
    gradientFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;options::IpoptOptions=IpoptOptions(),showplots=false)
Use JuMP+Ipopt for fitting.
"""
function gradientFit(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...;options::IpoptOptions=IpoptOptions(),showplots=false)
    
    ############################################
    # Check whether the sample should be evaluated or not based on errorcode
    ############################################
    if IE.ErrorCode != 0
        # for these cases, create a FitParams struct containing NaNs for the relevant parts
        # include the experimental spectra, but set fit to zero
        invalidParams = deepcopy(startingsolution)
        invalidParams.T = NaN
        for (key,value) in invalidParams.X
            invalidParams.X[key] = NaN
        end
        invalidParams.Spectra = SpectraStruct(IE.I./sum(IE.I),zeros(size(IE.I)),zeros(size(IE.I)),getWavenumberArray(startingsolution.WavenumberPolynomial,length(IE.I)),IE.ErrorCode,"Ipopt")
        if showplots
            plotSpec(invalidParams)
        end
        return invalidParams
    end

    global Ie
    ############################################
    # make a deep copy of startingsolution to not overwrite it
    ############################################
    params=deepcopy(startingsolution)

   
    ############################################
    # extract fitting variables and vectorize them
    ############################################
    x0,xlower,xupper,loc=vectorizeParams!(startingsolution,lower,upper,libs)

    ############################################
    # pass library and experiment to residual function using closure
    ############################################
    jumpResidual(x::T...) where {T<:Real} = calcResidual(IE,libs,params,loc,x...)


    ################## JUMP ##########################
    # IPOPT
    optimizer = Ipopt.Optimizer
    m = JuMP.Model(optimizer)
    for fname in fieldnames(IpoptOptions)
        set_optimizer_attribute(m,string(fname),getfield(options,fname))
    end
       
    ############################################
    # Construct the Model
    ############################################
    # splatting doesn't work if only one var is fitted
    # construct a special case for this.

    if length(x0)>1
        @variable(m, xlower[i]<=x[i=1:length(x0)]<=xupper[i],start=x0[i])
        register(m, :jumpResidual, length(x0), jumpResidual, autodiff=true)
        @NLobjective(m, Min, jumpResidual(x...))
    else
        @variable(m, xlower[1]<=x<=xupper[1],start=x0[1])
        register(m, :jumpResidual, 1, jumpResidual, autodiff=true)
        @NLobjective(m, Min, jumpResidual(x))
    end
    
    ############################################
    # print the Model
    ############################################
    if options.print_level > 1
        show(options)
        print(m)
    end
    
    
    ############################################
    # optimize the Model
    ############################################
    optimize!(m)
    
    ############################################
    # detect exit status of optimizer
    ############################################
    # define unknown quality
    quality = 0
    exit_reason = MOI.get(m, MOI.TerminationStatus())

    if exit_reason == MOI.OPTIMAL
        quality = 1
    end
    if exit_reason == MOI.TIME_LIMIT
        warn("Reached time limit during fit.")
        quality = -100
    end
    if exit_reason == MOI.INVALID_MODEL
        warn("Invalid number in NLP function or derivative detected!")
        quality = -101
    end

    ############################################
    # retrieve the objective value, corresponding x values and the status
    ############################################
    devectorizeParams!(params,JuMP.value.(x),loc)
    if options.print_level > 1
        println("####################")
        println(JuMP.value.(x))
        println("####################")
    end
    ############################################
    # final run to store everything in the params struct
    # and get rid of dual numbers
    ############################################
    jumpResidual(JuMP.value.(x)...)
    # also, put the spectra into the struct
    params.Spectra = SpectraStruct(IE.I./sum(IE.I),Ie,IE.I./sum(IE.I)-Ie,getWavenumberArray(params.WavenumberPolynomial,length(IE.I)),quality,"Ipopt")
    if showplots
       plotSpec(params)
    end

    return params
    
end

"""
    benchmarkGradient(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...)
Function to compute the gradient of the residual with the same argument structure as the fitting operations.
Use e.g with `@benchmark` or simply `@time`.
Helpful when comparing an optimized handling of an operation on duals with a working reference.
"""
function benchmarkGradient(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...)
    params=deepcopy(startingsolution)
    ############################################
    # extract fitting variables and vectorize them
    ############################################
    x0,xlower,xupper,loc=vectorizeParams!(startingsolution,lower,upper,libs)
    ############################################
    # pass library and experiment to residual function using closure
    ############################################
    testGradient(x::T...) where {T<:Real} = calcResidual(IE,libs,params,loc,x...)
    x0 = convert(Vector{Float64},x0)
    return ForwardDiff.gradient(x0 -> testGradient(x0...), x0)
end

"""
    benchmarkValue(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...)
Function to compute the value of the residual with the same argument structure as the fitting operations.
Use e.g with `@benchmark` or simply `@time`.
Helpful when comparing an optimized handling of an operation on duals with a working reference.
"""
function benchmarkValue(IE::PreProcSpectrum,startingsolution::FitParams,lower::FitParams,upper::FitParams,libs::Library...)
    params=deepcopy(startingsolution)
    ############################################
    # extract fitting variables and vectorize them
    ############################################
    x0,xlower,xupper,loc=vectorizeParams!(startingsolution,lower,upper,libs)
    ############################################
    # pass library and experiment to residual function using closure
    ############################################
    testGradient(x::T...) where {T<:Real} = calcResidual(IE,libs,params,loc,x...)
    x0 = convert(Vector{Float64},x0)
    return testGradient(x0...)
end

"""
    calcResidual(spec::PreProcSpectrum,libs,params::FitParams,loc,x::T...) where {T<:Real}

Calculates the resiudal between preprocessed experimental spectrum `spec` and a simulated spectrum defined by `params`.
Requires `loc` to interpret the vector `x` as fitting parameters to store in params struct.
!!! note
    All methods wihtin this function have to be compatible with ```T<:Real``` to be compatible with FowardDiff.Duals!
"""
function calcResidual(spec::PreProcSpectrum,libs,params::FitParams,loc,x::T...) where {T<:Real}

    global Ie
    
    ############################################
    # DeVectorization of fitting parameters
    ############################################
    # print("DeVectorization: ")
    devectorizeParams!(params,x,loc)
    
    ############################################
    # Generate the spectrum without taking square root
    ############################################
    # print("Spectrum generation: ")
    ω,I,ωres = generateSpectrum!(params,convert(T,1.),loc,libs...;takesqrt=false)

    ############################################
    # interpolate on final grid
    ############################################
    # print("Get Wavenumber Array: ")
    ωexp = getWavenumberArray(params.WavenumberPolynomial,length(spec.I))
    # if the wavenumber polynomial is fitted, it may create solutions outside the bounds of the library
    # so if this is the case (which should happen only in the beginning of an experiment anyway),
    # allow extrapolation. otherwise, the data is assumed to be in the bounds of the library and extrapolation can be omitted
    # this makes the code much faster.

    # print("Checking do_extrapolation: ")
    do_extrapolation = isFitParam(loc,:WavenumberPolynomial)

    # use a specialised interpolation function for dual numbers here
    # this is a huge improvement when interpolating on a coarser array (as opposed to just shifting as is done for the second pumping region above)
    # print("Interpolation: ")
    Ie = fastLinearInt(ω .+ params.WavenumberShift,I,ωexp;extrap=do_extrapolation)
    ############################################
    # take square root (if desired)
    ############################################
    #@. I = sqrt(I)
    # print("Square rooting: ")
    if CARSopt.fittype == :sqrt
        if minimum(Ie)<0
            # sometimes Ie can be slightly negative in parts, where the NR signal and the real part interfere
            # this is caused by splitting the convolution in two steps, introducing inaccuracies
            # return a warning if this get's really high
            # the threshold is: the absolute value of the minimum has 0.1% of the peak intensity.
            # everything above that cannot be explained with this convolution artefact
            if abs(minimum(Ie)/maximum(Ie)) > 1e-3
                warn("Spectrum generation returned a negative number. If this warning occurs often, something went wrong.")
                show(params)
            end
            @. Ie[Ie<0]=0   # as a workaround, set these negative values to zero.
            #Ie = Ie .- minimum(Ie) # this would be an alternative: shift the entire spectrum up. but this may cause issues with the mole fraction result
        end    
        @. Ie = sqrt(Ie)
    end

    ############################################
    # normalize to sum
    ############################################
    # print("Normalization: ")
    begin
    Ie = Ie ./ sum(Ie)
    Iexp = spec.I ./ sum(spec.I)    # can be omitted? spectra are already normalized!
    end
    ############################################
    # calculate the residual    
    ############################################
    #res = w .* (Iexp.- Ie)
    # print("Residual: ")
    res = (Iexp.- Ie)./length(Ie) # normalize the residual based on length of wavenumber array
    ## Weighting for fit
    # covariance matrix from Argon measurement // noise model
    # | σ_1^2 ...         |
    # | ...  \ \ \        |
    # |             σ_N^2 |
    # Dim: NxN, N = length(ω)

    # Cholesky decomposition -> w = chol(Γb^-1)
    # b_w => wb    

    ## Uncertainty estimation NxM, N=length(ω), M = length(x)
    # Sensitivity matrix at MLE,
    #       | db_1/dx1, db_1/dx2 ... |
    # J_w = |                        | ↓ ω
    #       |                        |

    # Γ(x) = (J_w^T J_w)^-1 -> convariance matix with all uncertainties
    # check if noise weights are available
    if spec.w[1] != 1.
        return sum(spec.w.*res.^2)  # return the weighted sum of squared residuals
    else
        return sum(res.^2)*1e6  # return the raw sum of squared residuals scaled to higher values
    end  
end

function getWavenumberArray(poly,len)
    idx = len-1:-1:0
    ω = @. poly[1] + poly[2]*idx + poly[3]*idx^2
    # ω = @. poly[1] * ones(len)
    # for j = 1:length(poly)-1
    #     @. ω += poly[j+1]*idx^j
    # end
    return ω
end

function interpSpec!(χ::AbstractArray{<:ForwardDiff.Dual{T,V,K}},w0::ForwardDiff.Dual{T,V,K},w1::ForwardDiff.Dual{T,V,K},χv::AbstractArray,X::ForwardDiff.Dual{T,V,K},pr::Int64) where {T,V,K}   
    # precalculate product of temperature weight and mole fraction
    w0x = w0*X
    w1x = w1*X
    # reinterpret the arrays
    χe = reinterpret(reshape, V, χ)
    w0xe = reinterpret(reshape,V,[w0x])
    w1xe = reinterpret(reshape,V,[w1x])
        
    # χv is an ordinary float32 - so this is simple: just multiple the value and partials with the actual number
    # this is about ~5 times faster than the version below.
    @avxt for j in 1:length(w0xe), i in 1:size(χe,2)
        χe[j,i,pr] += w0xe[j]*χv[i,1]+w1xe[j]*χv[i,2]
    end

    # no need to reconstruct the dual - these operations are inplace on χe → χ!

    # FALLBACK DIRECT IMPLEMENTATION USED FOR VERIFICATION
    #@. χ[:,pr] += w0x*χv[:,1]+w1x*χv[:,2]
end

function interpSpec!(χ::AbstractArray{T},w0::T,w1::T,χv::AbstractArray,X::T,pr::Int64) where {T<:AbstractFloat}   
    @avxt for j in 1:size(χ,1); χ[j,pr] += (w0*χv[j,1]+w1*χv[j,2]) * X;end
end


function signalsummation!(I::AbstractArray{<:ForwardDiff.Dual{T,V,K}},χM::AbstractArray{<:ForwardDiff.Dual{T,V,K}},χI::AbstractArray{<:ForwardDiff.Dual{T,V,K}},χR::AbstractArray{<:ForwardDiff.Dual{T,V,K}},χNR) where {T,V,K}
    
    # same trick as with the convolutions. reinterpret the arrays and handle them as ordinary Floats

    Ire = reinterpret(reshape, V, I)
    χMre = reinterpret(reshape, V, χM)
    χIre = reinterpret(reshape, V, χI)
    χRre = reinterpret(reshape, V, χR)
    
    χNRrev = ForwardDiff.value(χNR)
    χNRrep = reinterpret(reshape, V, [χNR])[2:end]

    # dimensions:
    # Ire[partials,wavenumber]
    # χre[partials,wavenumber,pumpregion]
    # χNRe[partials]
    
    # values
    @tturbo for ω = 1:size(χMre,2)
        #               χM1 +       χM2 +           4 χR1 χNR +         4 χNR χR2 +             2 χR1 χR2 +                 2 χI1 χI2 +             4 χNR^2  
        Ire[1,ω] = χMre[1,ω,1] + χMre[1,ω,2] + 4*χNRrev*χRre[1,ω,1] + 4*χNRrev*χRre[1,ω,2] + 2*χRre[1,ω,1]*χRre[1,ω,2] + 2*χIre[1,ω,1]*χIre[1,ω,2] + 4*χNRrev^2
    end
    
    # partials: (u*v)' = u'v + v'u
    @tturbo for p = 2:size(χMre,1),ω = 1:size(χMre,2)
        Ire[p,ω] = χMre[p,ω,1] + χMre[p,ω,2] + 4*χNRrev*χRre[p,ω,1] + 4*χNRrep[p-1]*χRre[1,ω,1] + 4*χNRrev*χRre[p,ω,2] + 4*χNRrep[p-1]*χRre[1,ω,2] + 2*χRre[1,ω,1]*χRre[p,ω,2] + 2*χRre[p,ω,1]*χRre[1,ω,2] + 2*χIre[1,ω,1]*χIre[p,ω,2] + 2*χIre[p,ω,1]*χIre[1,ω,2] + 8*χNRrep[p-1]*χNRrev
    end

    

    #= FALLBACK IMPLEMENTATION USED FOR VERIFICATION
    #println(I[1])
    #fill!(I,0)
    
    #print("DUAL Summation: ")
        if size(χM,2) == 1  # single pump case: reduced from dual pump case with 1 == 2
            @. I = 2*(χM[:,1] + 2*χNR*χR[:,1] + χNR^2) + 
                2*(χR[:,1]^2 + 2*χNR*χR[:,1] + χNR^2 + χI[:,1]^2)

        else                # dual pump case. divide amplitudes by 2? or is the non resonant background added twice and everything is fine?
            #@. I = χM[1,:] + 2*χNR*χR[1,:] + χNR^2 + χM[2,:] + 2*χNR*χR[2,:] + χNR^2 + 2*(χR[1,:]*χR[2,:] + χNR*χR[1,:] + χNR*χR[2,:] + χNR^2 + χI[1,:]*χI[2,:])
            @. I = χM[:,1] + χM[:,2] + 4*χNR*χR[:,1] + 4*χNR*χR[:,2] + 2*χR[:,1]*χR[:,2] + 2*χI[:,1]*χI[:,2] + 4*χNR^2
        end
    println(I[1])
      =#  
end

function signalsummation!(I::AbstractArray{T},χM::AbstractArray{T},χI::AbstractArray{T},χR::AbstractArray{T},χNR::T) where {T<:AbstractFloat}   
        #print("AVXT Summation: ")
        if size(χM,2) == 1  # single pump case: reduced from dual pump case with 1 == 2
            @tturbo for j in 1:size(I,1)
                I[j] = 2*(χM[j,1] + 2*χNR*χR[j,1] + χNR^2) + 
                    2*(χR[j,1]^2 + 2*χNR*χR[j,1] + χNR^2 + χI[j,1]^2)
            end
        else           
            @tturbo for j in 1:size(I,1)
            #I[j] = χM[j,1] + 2*χNR*χR[j,1] + χNR^2 + 
            #    χM[j,2] + 2*χNR*χR[j,2] + χNR^2 +
            #    2*(χR[j,1]*χR[j,2] + χNR*χR[j,1] + χNR*χR[j,2] + χNR^2 + χI[j,1]*χI[j,2])
                # χM1 + 4 χR1 χNR + 2 χR1 χR2 + 2 χI1 χI2 + 4 χNR^2 + 4 χNR χR2 + χM2
            I[j]=χM[j,1] + 4*χNR*χR[j,1] + 
               χM[j,2] + 4*χNR*χR[j,2] +
               2*χR[j,1]*χR[j,2] +
               2*χI[j,1]*χI[j,2] +
               4*χNR^2
            end
        end

end

#function Base.:*(a::ForwardDiff.Dual{T,V,K},b::ForwardDiff.Dual{T,V,K}) where {T,V,K}
#    va,vb = ForwardDiff.value(a), ForwardDiff.value(b)
#    pa,pb = ForwardDiff.partials(a), ForwardDiff.partials(b)
#    return ForwardDiff.Dual{T,V,K}(va*vb, va*pb+pa*vb)
#end

## functions for preallocations that use caching

function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx, last_dT = nothing, nothing, nothing
    function foo_i(i::T,x...) where {T<:Real}
        if T == Float64
            if x != last_x
                last_x, last_f = x, foo(i,x...)
            else
                fill!(last_f[1],0)
                fill!(last_f[2],0)
                fill!(last_f[3],0)
                fill!(last_f[4],0)
            end
            return last_f
        else
            # if the type changed, also reallocate
            # this is important e.g. when the Hessian is computed
            # or multiple successive runs with different number of 
            # free parameters are used
            if x != last_dx || last_dT != T
                last_dT = T
                last_dx, last_dfdx = x, foo(i,x...)
            else
                fill!(last_dfdx[1],0)
                fill!(last_dfdx[2],0)
                fill!(last_dfdx[3],0)
                fill!(last_dfdx[4],0)
            end
            return last_dfdx
        end
    end
    return (x...) -> foo_i(x...)
end

function preallocate(i::T,num_pr,lω) where {T<:Real} 
    χM,χR,χI,I = _preallocate(i,num_pr,lω)
    return χM,χR,χI,I
end

function _preallocate(i::T,num_pr,lω,) where {T<:AbstractFloat}
    χM = zeros(T,lω,num_pr)
    χR = zeros(T,lω,num_pr)
    χI = zeros(T,lω,num_pr)
    I = zeros(T,lω)
    return χM,χR,χI,I
end

function _preallocate(i::T,num_pr,lω) where {T<:ForwardDiff.Dual}
    χM = zeros(T,lω,num_pr)
    χR = zeros(T,lω,num_pr)
    χI = zeros(T,lω,num_pr)
    I = zeros(T,lω)
    return χM,χR,χI,I
end

memoized_preallocate = memoize(_preallocate,4)

function gaussiankernelT(ωres::T,wid::T) where {T<:Real}
    # get a grid 4 times the width
    # at this width, g[1]/maximum(g) ≈ 1.5e-5 which should be close enough
    x = collect(T,0:ωres:2*wid)
    # force it to be symmetric. not necessarily true if -2:ωres:2
    x = vcat(reverse(x[2:end]),x)
    #x = collect(T,-2:ωres:2)
    # get gaussian
    c = convert(T,0.6006)
    g = @. exp(-((x)./(c.*wid)) .^2);
    # normalize
    g = g./sum(g);
    return g
end
