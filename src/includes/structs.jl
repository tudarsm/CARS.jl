
"""
    Base.@kwdef mutable struct IpoptOptions
        max_iter::Int = 1000                        # Abort after a certain number of iterations
        max_wall_time::Float64 = 120.               # Abort after reaching a certain wall time
        tol::Float64 = 1e-8                         # Abort when reaching this tolerance (see Ipopt Documentation, 1e-7 is the default value)
        print_level::Int = 3                        # Amount of information during optimization.
        mu_strategy::String = "adaptive"            # "monotone" or "adaptive".
    end
Specify options for [`gradientFit`](@ref)
Some options passed to the Ipopt Optimizer.
If you wish to extend this, check out the official [Ipopt documentation](https://coin-or.github.io/Ipopt/OPTIONS.html)
Values in this struct are all passed to the optimizer by name.
"""
Base.@kwdef mutable struct IpoptOptions
    max_iter::Int = 1000                        # Abort after a certain number of iterations
    max_wall_time::Float64 = 120.               # Abort after reaching a certain wall time
    tol::Float64 = 1e-8                         # Abort when reaching this tolerance (see Ipopt Documentation, 1e-7 is the default value)
    print_level::Int = 3                        # Amount of information during optimization.
    mu_strategy::String = "adaptive"            # "monotone" or "adaptive".
end

"""
    Base.@kwdef mutable struct BlackBoxOptions
        MaxFuncEvals::Int = 1000                    # Stop optimizing after reaching this number of iterations
        TargetFitness::Float64 = 1e-7                 # Stop optimizing after reachting this Fitness level
        # FitnessTolerance not yet implemented
        FitnessTolerance::Float64 = 1e-7            # If StoppingCriterion==:tol, stop if change between iterations is lower than this value
        TraceMode::Symbol = :silent                 # Output of blackboxoptim during optimization, see docs of BlackBoxOptim.jl for choices
        Method::Symbol = :adaptive_de_rand_1_bin    # Method used for optimization, see docs of BlackBoxOptim.jl for choices
        StoppingCriterion::Symbol = :best           # If set to :best, iteration will stop if BestFitness is reached. If set to :tol, it will stop if difference between iterations is smaller than FitnessTolerance
        PopulationSize::Int = 5                    # Population Size for methods which actually use a population. No effect for others.
    end
Specify options for [`blackBoxFit`](@ref)
"""
Base.@kwdef mutable struct BlackBoxOptions
    MaxFuncEvals::Int = 1000                    # Stop optimizing after reaching this number of iterations
    TargetFitness::Float64 = 1e-7                 # Stop optimizing after reachting this Fitness level
    # FitnessTolerance not yet implemented
    FitnessTolerance::Float64 = 1e-7            # If StoppingCriterion==:tol, stop if change between iterations is lower than this value
    TraceMode::Symbol = :silent                 # Output of blackboxoptim during optimization, see docs of BlackBoxOptim.jl for choices
    Method::Symbol = :adaptive_de_rand_1_bin    # Method used for optimization, see docs of BlackBoxOptim.jl for choices
    StoppingCriterion::Symbol = :best           # If set to :best, iteration will stop if BestFitness is reached. If set to :tol, it will stop if difference between iterations is smaller than FitnessTolerance
    PopulationSize::Int = 5                    # Population Size for methods which actually use a population. No effect for others.
end

"""
    mutable struct CARSOptions
Specifies global options that control certain aspects of CARS.jl
- `fittype::Symbol`: either `:sqrt` to take the square root of the signal intensity (after preprocessing) or anything else (use the signal directly). Affects [`preprocess_spectra`](@ref), [`calcResidual`](@ref), [`generateSpectrum!`](@ref) and [`noiseModel!`](@ref).
- ```γf::AbstractFloat```: Controls rediscretization after convolution with the compression kernel in [`generateLibrary`](@ref) and [`simplifyLibrary!`](@ref). Higher values means higher resolution of wavenumber axis. Defaults to 5.
- ```interp```: Control interpolation scheme for some of the operations. Stick with Linear().
- ```conv_margin```: Specifies a safety margin for the convolution with the apparatus function to avoid edge effects after simplifying a library. Applied as a factor to the instrumental linewidth.
"""
Base.@kwdef mutable struct CARSOptions
    fittype::Symbol = :sqrt                   # instruct the fit to either fit the intensity directly or the sqrt of the intensity
    γf::AbstractFloat = 5                     # this factor determines how many points are used inside the laser linewidth for rediscretization. a higher value leads to less interpolation artifacts, a lower value leads to faster fitting.
    interp = Linear()                    # specifies the order of interpolation used when interpolating on another wavenumber grid. Possible options: Constant(), Linear(), Quadratic(), Cubic(). Note: This is not used for spectral interpolations between temperatures!
    conv_margin = 5                      # specifies a safety margin for the convolution with the apparatus function to avoid edge effects after simplifying a library. applied as a factor to the instrumental linewidth.
end


mutable struct SpectraStruct
    Exp::Vector
    Fit::Vector
    Res::Vector
    omega::Vector
    quality::Int
    optimizer::String
end

"""
    mutable struct PreProcSpectrum
        I::AbstractVector{Float64}          # area normalized intensity
        w::AbstractArray{Float64}           # weight used for fit
        A::Float64                          # original area (used to scale the variance for noise weighting)
        ErrorCode::Int8                     # specify errorcode after preprocessing (0: no error, 1:saturated, 2: signal too low, 3: breakdown)
    end
"""
mutable struct PreProcSpectrum
    I::AbstractVector{Float64}          # area normalized intensity
    w::AbstractArray{Float64}           # weight used for fit
    A::Float64                          # original area (used to scale the variance for noise weighting)
    ErrorCode::Int8                     # specify errorcode after preprocessing (0: no error, 1:saturated, 2: signal too low, 3: breakdown)
end

# Base.@kwdef allows to set defaults for some of the parameters
    # still, at least T, X and LineWidth have to be set.
    # if adding or removing parameters here, don't forget to add them to the Base.copy(s::FitParams) implementation below!
    
"""
    Base.@kwdef mutable struct FitParams
        T                          # Temperature as value
        X::Dict{String,Any}        # Dictionary of Mole Fractions of resonant species
        chiNR=0                      # Non-resonant signal of resonant species. No need to set this manually, will be filled with correct values automatically
        chiB=8.5                     # Non-resonant susceptibility of the buffer gas    
        AdditionalFrequencyOffset=0  # For DP CARS Spectra. Additional shift (to what's already in the library) of Pump region 2 relative to pump region 1
        LineWidth::Dict{String,Any} # Dict of Laser linewidth(s), "P1" for pump region 1, "P2" for pump region 2
        Instrumental=0             # Instrumental LineWidth of spectrometer
        WavenumberPolynomial::Vector{Any}=[2200.,0.5,0.]     # Vector of polynomial parameter to fit the wavenumber axis, Syntax [a, b, c] will be interpreted as a+x^1*b+x^2*c
        WavenumberShift=0          # Shift of calibrated wavenumber axis (zero order)
        Spectra::SpectraStruct=SpectraStruct([0.],[0.],[0.],[0.],0,"unknown")     # After fitting, this will contain the original experimental spectrum (after area normalization) along with the fit
    end
Struct containing all relevant degrees of freedom for simulatin/fitting spectra.
Every parameter in here (with the exception of `Spectra`) can be used in the fit.
"""
    Base.@kwdef mutable struct FitParams
        T                          # Temperature as value
        X::Dict{String,Any}        # Dictionary of Mole Fractions of resonant species
        chiNR=0                      # Non-resonant signal of resonant species. No need to set this manually, will be filled with correct values automatically
        chiB=8.5                     # Non-resonant susceptibility of the buffer gas    
        AdditionalFrequencyOffset=0  # For DP CARS Spectra. Additional shift (to what's already in the library) of Pump region 2 relative to pump region 1
        LineWidth::Dict{String,Any} # Dict of Laser linewidth(s), "P1" for pump region 1, "P2" for pump region 2
        Instrumental=0             # Instrumental LineWidth of spectrometer
        WavenumberPolynomial::Vector{Any}=[2200.,0.5,0.]     # Vector of polynomial parameter to fit the wavenumber axis, Syntax [a, b, c, ..] will be interpreted as a*x^1+b*x^2+c*x^3...
        WavenumberShift=0          # Shift of calibrated wavenumber axis (zero order)
        Spectra::SpectraStruct=SpectraStruct([0.],[0.],[0.],[0.],0,"unknown")     # After fitting, this will contain the original experimental spectrum (after area normalization) along with the fit
    end

"""
    mutable struct Library
        species::Vector{String}  # name of resonant species in this library
        ROI::Vector                 # Region of Interest (excluding extension)
        AdditionalFrequencyOffset     # FrequencyOffset relative to PumpRegion 1
        T::Vector           # Vector of temperatures
        P::Float32          # Pressure
        ω::Vector{Float32}  # Wavenumberarray after compression
        ωres::Float32       # Wavenumber resolution used in library
        γc::Float32          # width of preconvolution kernel
        χC::Matrix{Array{Float32}} # compressed complex cross-terms
        χM::Array{Float32} # Compressed magnitude
        χI::Array{Float32} # Compressed imaginary components
        χR::Array{Float32} # Compressed real components
        χNR::Vector{Float32}  # (Purely real) non-resonant susceptibility of species in library. Already converted in cm3/molecule
    end

Struct containing all relevant data for a single pumping region.
"""
    mutable struct Library
        species::Vector{String}  # name of resonant species in this library
        ROI::Vector                 # Region of Interest (excluding extension)
        AdditionalFrequencyOffset     # FrequencyOffset relative to PumpRegion 1
        T::Vector           # Vector of temperatures
        P::Float32          # Pressure
        ω::Vector{Float32}  # Wavenumberarray after compression
        ωres::Float32       # Wavenumber resolution used in library
        γc::Float32          # width of preconvolution kernel
        χC::Matrix{Array{Float32}} # compressed complex cross-terms
        χM::Array{Float32} # Compressed magnitude
        χI::Array{Float32} # Compressed imaginary components
        χR::Array{Float32} # Compressed real components
        χNR::Vector{Float32}  # (Purely real) non-resonant susceptibility of species in library. Already converted in cm3/molecule
    end

    # also define what copy means here
    # because it has to be an immutable struct, copying has to be defined
    # to simply create upper and lower bounds from a template
    # use copy everywhere. for scalars, it is not required (but doesn't hurt),
    # for dicts and arrays and other immutable containers, it is.
    Base.copy(s::SpectraStruct) = SpectraStruct(deepcopy(s.Exp),deepcopy(s.Fit),deepcopy(s.Res),deepcopy(s.omega),deepcopy(s.quality),deepcopy(s.optimizer))
    Base.copy(s::FitParams) = FitParams(deepcopy(s.T),deepcopy(s.X),deepcopy(s.chiNR),deepcopy(s.chiB),deepcopy(s.AdditionalFrequencyOffset),deepcopy(s.LineWidth),deepcopy(s.Instrumental),deepcopy(s.WavenumberPolynomial),deepcopy(s.WavenumberShift),deepcopy(s.Spectra))
    Base.copy(s::Library) = Library(deepcopy(s.species),deepcopy(s.ROI),deepcopy(s.AdditionalFrequencyOffset),deepcopy(s.T),deepcopy(s.P),deepcopy(s.ω),deepcopy(s.ωres),deepcopy(s.γc),deepcopy(s.χC),deepcopy(s.χM),deepcopy(s.χI),deepcopy(s.χR),deepcopy(s.χNR))

    # when showing a variable, make the output nice and pretty:
    function Base.show(params::FitParams)
        for fname in fieldnames(FitParams)
            if fname == :Spectra
                println("Not displaying $fname.")
                continue
            end
            if isa(getfield(params,fname),Dict)
                for key in keys(getproperty(params,fname))
                    println("$fname[\"$key\"] = $(getfield(params,fname)[key])")
                end
            else
                println("$fname = $(getfield(params,fname))")
            end
        end
    end

    function Base.show(opt::IpoptOptions)
        for fname in fieldnames(IpoptOptions)
            println("$fname = $(getfield(opt,fname))")
        end
    end

    function Base.show(opt::BlackBoxOptions)
        for fname in fieldnames(BlackBoxOptions)
            println("$fname = $(getfield(opt,fname))")
        end
    end