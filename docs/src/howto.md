```@meta
CurrentModule = CARS
```

# How to ...
This section gives hints for developers on how to extend functionality to CARS.jl and how to perform typical tasks

## ... cite CARS.jl
You can cite [this](https://doi.org/10.1002/jrs.6639) publication which describes the method.

## ... contribute to CARS.jl
- Clone the repository
- Open the root folder (the one containing `src/`) in VS Code with the Julia extension installed
- Open a REPL with ALT+J+O (without releasing ALT)
- Enter `] dev .` which activates the local repository for development
- Start modifying, commit, push/create PR.

## ... build the docs locally before pushing to the repo
Open a command line and navigate to the ```CARS.jl``` folder and start julia:
```console
cd yourLocalCARS.jlFOLDER
julia --project=docs/
```
Then, in Julia, execute the following commands:
```julia
using Revise
using Pkg;Pkg.develop(PackageSpec(path=pwd()));Pkg.instantiate();
include("docs/make.jl")
```
This will build the docs in the build directory. Never push this to the repository. It is in the .gitignore file for a reason. Documentations are published to Gitlab/Github pages after pushing to the repository. Leave the command window open and repeat the ```include("docs/make.jl")``` command to update the pages.
!!! note 
    The first run will probably be very slow. Subsequent runs are blazing fast.

To view the newly build documentation, I recommend to use Live Server addon by Ritwick Dey in VSCode. Right click on the index.html (in VSCode) and click Open with Live Server. Upon recreating the docs using the above method, the page automatically reloads.

## ... add more species
In terms of how to generate theoretical susceptibilities, this is way beyond the scope of this documentation and in fact not the job of CARS.jl.
If it is a simple diatomic molecule, check out DiaCARS.jl on how these are modelled. For more complex molecules, such as CO2, H2O, NH3, check the relevant literature or come up with your own idea. This is not trivial.

If you do have a model for the molecule you are interested in, it can be used in CARS.jl simply using the following:

```julia
using CARS
import YourFancyMolecule

# generate a function that is only depending on temperature
yourfancymolecule(t) = YourFancyMolecule(t,yourmodelparameters,...)

# generate the library (if it is the only species) by passing the function to generateLibrary
T = 280:5:2000              # temperature range of interest
library = generateLibrary(T,γc,yourfancymolecule)
```

For this to work, the function has to return the following parameters in exactly that order:

```julia
struct MolParams
    CHINR::Float64
end

function YourFancyMolecule(t,yourmodelparameters,...)

    # do your simulation here

    # molecularParameters is a struct with at least CHINR per molecule
    # e.g. for nitrogen, the value is 8.5e-18 cm3/(erg*amagat)
    # the value that has to be stored is 8.5 / NA * Molvol = 3.1635e-19 (note that the 1e-18 is omitted)
    # NA = 6.022...e23
    # Molvol = 2.2413e4
    # check out DiaCARS.jl for reference
    molecularParameters = MolParams(1.12345e-19)
    return χR,χI,ω,species,molecularParameters      
end

```

## ... add fitting parameters
To understand how CARS.jl interprets something as a fitting parameter, look at the ```FitParams``` struct:
```julia
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
```
Furthermore, fit parameters are defined at runtime by specifying a startingsolution as well as upper and lower bounds:
```julia
startingsolution = FitParams(T=800,X=Dict("N2"=>0.7,"CO2"=>0.08),LineWidth=Dict("P1"=>0.2,"P2"=>.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=1);
lower=copy(startingsolution);upper=copy(startingsolution);
lower.X["N2"]=0.6
upper.X["N2"]=0.75
```

When running in fitting mode, CARS.jl runs these three structs through [`vectorizeParams!`](@ref) and identifies all fields that have a lower and upper bound **differing** from the starting solution as a fit parameter.
It returns the starting value, lower and upper bound in a vector alongside a ```loc``` array which instructs the spectrum generator how to interpret that value, i.e. if it is temperature, mole fraction or another degree of freedom.

So in short, to add a new fit parameter, simply add it to the [`FitParams`](@ref) struct and use it in  [`calcResidual`](@ref) or  [`generateSpectrum!`](@ref).
Make sure that the added functionality can handle ```FowardDiff.Dual```, i.e. ```T<:Real``` datatypes for the input. Check out the next section for examples how to do that efficiently.

## ... add new preprocessing
Recall how currently the scheme of loading and preprocessing works:
```julia
signal = load_spectra("./src/examples/data/210622A163.SPE")     # load the signal
nr = load_spectra("./src/examples/data/210622A144.SPE")         # load the non-resonant signal (e.g. from Argon)
preprocessed_signals = preprocess_spectra(signal,nr,[600 1000]) # preprocess using defaults
```
[`preprocess_spectra`](@ref) assumes the datatype [`RawSpectra`](@ref) and returns an array of [`PreProcSpectrum`](@ref).
If you want to do fancier preprocessing, such as referencing the signal to single-shot non-resonant signals or different background removals, the only requirement is that the return value of your new preprocessing is also an array of [`PreProcSpectrum`](@ref)

As a default, `PreProcSpectrum.w` can be `[1.]`, which disables noise weighting. `PreProcSpecttrum.ErrorCode` can be any value of type `Int8`. Only those spectra that have a value of `0` here will be evaluated.

## ... implement efficient handling of Duals for the gradient based fitter
!!! warning 
    This requires some experience with Julia and an understanding of how automatic mode differentiation works. Please make sure you understand what you are planning on doing before you actually do so. Benchmark and validate your implementation before including it in the production version!
    Before proceeding, read and understand the [documentation](https://juliadiff.org/ForwardDiff.jl/stable/) of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) and [this section](https://jump.dev/JuMP.jl/stable/manual/nlp/#Common-mistakes-when-writing-a-user-defined-function) of the JuMP documentation.

Ok. You now know the requirements of your method to accept `T<:Real` datatypes and decided to include a new degree of freedom into either one of [`calcResidual`](@ref) or [`generateSpectrum!`](@ref).
The need for some trickery really depends on how expensive the computation is.
For example, if you just include a vertical shift with a constant number to the entire spectrum, the following is likely totaly fine:
```julia
function addVerticalShift!(spectrum:AbstractArray{T},shift::T) where {T<:Real}
    @. spectrum += shift
end
```
Note, that for efficiency, this already uses an inplace operation on the variable `spectrum`, as indicated by the exclamation point in the name of the function.
The reason why this naive implementation is unproblematic, is that the sum of two duals is simply summing the value and partials independently.
However, more complex functions that operate on multiple arrays, such as convolutions, it is not that straightforward.
Using the default `DSP.conv()` method for example, does not work for Dual numbers. Also, using `@tturbo` or `@avxt` from the package [LoopVectorizations.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) does not allow using Duals directly.
For an example on how to work around this issue, check out this implementation of [`fastLinearInt`](@ref):

```julia
function fastLinearInt(wlib::AbstractArray,data::AbstractArray{T},wexp::AbstractArray{T};extrap=false)  where {T<:AbstractFloat} 
    right,left,w = findNeighbours(wexp,wlib,extrap=extrap)
    interpolated = @avxt @. (w)*data[left]+(1-w)*data[right]
    return interpolated
end

function fastLinearInt(wlib::AbstractArray{R},data::AbstractArray{<:ForwardDiff.Dual{T,V,K}},wexp::AbstractArray{<:AbstractFloat};extrap=false)  where {R,T,V,K}  
   # in this function data, is always of dual type
   # wlib is dual type if wavenumber shift is fitted
   # wexp is always float here, because dual case is handled in another implementation below
   
   # reinterpret values and partials
   dataE = reinterpret(reshape, V, data)

   # preallocate output array
   out = Array{eltype(data)}(undef,length(wexp))
   # and reinterpret it here to fill with data.
   outE = reinterpret(reshape, V, out)

   if R<:ForwardDiff.Dual
        wlibE = reinterpret(reshape,V,wlib)
        # get index where wlibE has a partial derivative of 1.
        idx = findlast(wlibE[:,1].==1.)
        # remove the partials, we just need the wavenumber array for interpolation
        wlibE = wlibE[1,:]
   else
        wlibE = wlib
   end
   
   # find neighbors and weights
   right,left,w,dw = findNeighbours(wexp,wlibE,extrap=extrap)

   # interpolate values and partials
   @tturbo for i in axes(dataE,1), j in axes(outE,2)
        outE[i,j] = (w[j])*dataE[i,left[j]]+(1-w[j])*dataE[i,right[j]]
   end

   # special case: w depends on partial of wlib
   # so if wlib is a Dual (ie. WavenumberShift is a FitParameter)
   # account for the partial derivative
    if R<:ForwardDiff.Dual
        for j in axes(outE,2)
                outE[idx,j]=dw*dataE[1,left[j]]-dw*dataE[1,right[j]]
        end
    end

   return out   # out has contents of outE!
end
```

Note a couple of things:
- The implementation makes use of [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/) to decide when to use the method for interpolating an array of `Dual` or an array of `AbstractFloat`.
- We reinterpret the array `data::AbstractArray{<:ForwardDiff.Dual}` which essentially means, we can manipulate the values and partials of this array explicitly.
- We prellocate the output `out` to be the same datatype as `data`, but with the new length determined by `wexp`
- The interpolation is done for both values and partials using this for-loop where the index `i` iterates over the values and partials and `j` is the wavenumber direction:

```julia
# interpolate values and partials
@tturbo for i in axes(dataE,1), j in axes(outE,2)
    outE[i,j] = (w[j])*dataE[i,left[j]]+(1-w[j])*dataE[i,right[j]]
end
```
But: one of the partials is the partial derivative with respect to a number that `w` depends on because we use this function sometimes to fit a wavenumber shift which has an effect on the weights for the interpolation.
If we look at the definition of [`findNeighbours`](@ref):
```julia
function findNeighbours(wexp,wlib;extrap=false)
    r=1
    qp = Array{Int32}(undef,length(wexp))
    w = Array{Float64}(undef,length(wexp))
    dw = wlib[2]-wlib[1]
    # auto detect reversal
    if wexp[2]>wexp[1]
        iterator = 1:length(wexp)
    else
        iterator = length(wexp):-1:1
    end
    for i in iterator
        for j in r:length(wlib)
            if wexp[i]<wlib[j]
                r=j
                qp[i]=j
                w[i]=(wlib[j]-wexp[i])/dw
                break
            else
                qp[i]=length(wlib)
                w[i]=0
            end
        end
    end

    if extrap
        qp[qp.==1] .= 2
    end
    return qp,qp.-1,w,1/dw
end
```
and have a closer look at this line in the body of the for loop:
```julia
                w[i]=(wlib[j]-wexp[i])/dw
```
we realize, that the derivative of `w` with respect to `wlib` (and thus to the wavenumber shift) is simply `1/dw`.
We can use this to compute the partial derivative as a special case when interpolating our array `data` for cases where the wavenumber shift is fitted:
```julia
if R<:ForwardDiff.Dual
    for j in axes(outE,2)
            outE[idx,j]=dw*dataE[1,left[j]]-dw*dataE[1,right[j]]
    end
end
```

!!! tip "More References" 
    This is, in the beginning, a little confusing. But it makes a night and day difference for computationally expensive operations. In this example, the speed up between using a naive linear interpolation based on [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) and this approach is about a factor of 30 when handling duals.
    If you would like to see more examples, check out `includes/convolutions.jl` or `includes/fit_spectra.jl`. Specifically `signalsummation!` provides a nice entry point for this way of handling partial derivatives semi-analytically.

!!! tip
    Always have a fallback solution implemented, at least for verification of your implementation!

## ... benchmark and validate your implementation
Since CARS.jl is all about efficiency, you might want to benchmark the adaptions or optimizations you were implementing.
You can use the methods [`benchmarkGradient`](@ref) and [`benchmarkValue`](@ref) for this, which take the same arguments as [`gradientFit`](@ref) and [`blackBoxFit`](@ref) without the options. For example, you can do the following:

```julia
startingsolution = FitParams(T=800.,X=Dict("N2"=>0.7,"CO2"=>0.08,"O2"=>0.1),LineWidth=Dict("P1"=>0.7,"P2"=>.7),AdditionalFrequencyOffset=-1.,Instrumental=1.,WavenumberShift=-5);
startingsolution.WavenumberPolynomial = [2260, 0.285, 0]
lower=copy(startingsolution);upper=copy(startingsolution);
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.0
upper.X["CO2"]=0.5
lower.X["O2"]=0.
upper.X["O2"]=0.5
lower.T=300
upper.T=2200
lower.Instrumental = 0.5
upper.Instrumental = 1.5
lower.LineWidth["P1"]=0.25
upper.LineWidth["P1"]=0.7
lower.LineWidth["P2"]=0.25
upper.LineWidth["P2"]=1.6
lower.AdditionalFrequencyOffset = -5
upper.AdditionalFrequencyOffset = 5
lower.WavenumberPolynomial = [2258, 0.45, 0]
upper.WavenumberPolynomial = [2263, 0.5, 0]
lower.WavenumberShift =-6
upper.WavenumberShift=1
# run either this:
@benchmark CARS.benchmarkGradient(spectra[1],startingsolution,lower,upper,lib1,lib2)
# or this:
@benchmark CARS.benchmarkValue(spectra[1],startingsolution,lower,upper,lib1,lib2)
```

!!! tip "Validate results"
    The functions return the gradient or value respectively. Use this to ensure that your implementation is done correctly by comparing the results of your optimized implementation and a known reference.
