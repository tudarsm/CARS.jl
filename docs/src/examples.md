# Examples
This section covers typical examples. More can be found in the `src/examples` directory of this repository. Just copy them somewhere, [install the required packages](installation.md) and run the code.


## Generate Libraries
### [Dual Pump, Single Species per Pump Region](@id dp_ss_example)
In this example, we use DiaCARS to simulate Dual-Pump CARS spectra containing only N2 and O2:

```julia
using CARS
import DiaCARS
## Dual Pump N2/O2-CARS
ROI1 = [2200 2700]              # Pump region 1
FrequencyOffset = 985;          # Offset for pump region two. Pump1-FrequencyOffset=Region2
ROI2 = ROI1 .- FrequencyOffset  # Pump region 2
T=280:20:2400                   # Temperature range for library
γc = 0.1                        # preconvolution kernel
ωres = 5e-4                     # resolution of theoretical susceptibility
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres) # N2, Pump 1
o2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres) # O2, Pump 2

pump1_library = generateLibrary(T,γc,n2_pump1); 
pump2_library = generateLibrary(T,γc,o2_pump2); 
```

You should see an output like this:
```julia
Generating a library for species N2 in range 2200.0-2700.0 cm-1
Array size before compression: 1000001 (Original resolution: 0.0005 cm-1)
Array size after compression: 25000 (New resolution: 0.020000000298023225 cm-1)
Number of required complex cross-terms: 0
Estimated library size: 31.0 MB

100.0%┣████████████████████████████████┫ 107/107 [00:07<00:00, 15it/s]

Generating a library for species O2 in range 1215.0-1715.0 cm-1
Array size before compression: 1000001 (Original resolution: 0.0005 cm-1)
Array size after compression: 25000 (New resolution: 0.020000000298023225 cm-1)
Number of required complex cross-terms: 0
Estimated library size: 31.0 MB

100.0%┣████████████████████████████████┫107/107 [00:07<00:00, 14it/s]
```

### Generate Multi-Species Libraries

```julia
using CARS
import DiaCARS
import CO2CARS
## Dual Pump N2/O2-CARS
ROI1 = [2200 2700]              # Pump region 1
FrequencyOffset = 985;          # Offset for pump region two. Pump1-FrequencyOffset=Region2
ROI2 = ROI1 .- FrequencyOffset  # Pump region 2
T=280:20:2400                   # Temperature range for library
γc = 0.1                        # preconvolution kernel
ωres = 5e-4                     # resolution of theoretical susceptibility
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres) # N2, Pump 1
o2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres) # O2, Pump 2
co2_pump2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres)         # CO2, Pump 2

pump1_library = generateLibrary(T,γc,n2_pump1); 
pump2_library = generateLibrary(T,γc,o2_pump2,co2_pump2); # Pass as many species here as you like
```

## Simulate Single Spectra
Generating spectra is done with ```generateSpectrum!()``` together with the ```FitParams``` struct and the libraries.
This struct contains all physical parameters required for simulating a spectrum and will also be used to provide starting solutions as well as lower and upper bounds for the fitting process later.
For this example, use the libraries generated in [this example](@ref dp_ss_example).

### Single Pump
As stated in the [limitations](@ref limitations), single pump spectra are not directly supported yet. But we can workaround this by acknowleding that single-pump CARS is the same as dual-pump CARS with the same region of interest. Simply pass the `pump1_library` twice to the functions for generating or fitting a spectrum. Make sure to use the same linewidths for both pumping regions.

```julia
exp_params = FitParams(T=2000,X=Dict("N2"=>0.7,),LineWidth=Dict("P1"=>0.7,"P2"=>0.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=1);
ωE,IE = generateSpectrum!(exp_params,pump1_library,pump1_library)
plotSpec(ωE,IE)
```
### Dual Pump
```julia
exp_params = FitParams(T=2000,X=Dict("N2"=>0.7,"O2"=>0.2),LineWidth=Dict("P1"=>0.2,"P2"=>1),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=1);
ωE,IE = generateSpectrum!(exp_params,pump1_library,pump2_library)
plotSpec(ωE,IE)
```

## Fit Experimental Spectra
!!! danger
    Create simulated spectra for this example to avoid dependence on non public CO2CARS.jl

Supplied in the `src/examples/data` directory are two files:
- `210622A163.SPE` contains N2/CO2 CARS spectra
- `210622A144.SPE` is the non-resonant signal measured before the experiment

These samples were recorded at 20 Hz with intermittent capturing of the background for subtraction.

At first, we load the relevant modules:
```julia
using CARS
using Plots
import DiaCARS
import CO2CARS
```

Then, we load the data and plot it for visualization:
```julia
signal = load_spectra("./src/examples/data/210622A163.SPE")
nr = load_spectra("./src/examples/data/210622A144.SPE")     
l = @layout [a b;c d]
p1=plot(signal.I)
p2=plot(nr.I)
p3=plot(signal.BG)
p4=plot(nr.BG)
plot(p1,p2,p3,p4,layout=l,legend=:none)
```

Afterwards, we can preprocess the signal which involves background correction with the intermittently recorded background, normalization to the non-resonant background and error checking.
For the latter, CARS.jl has reasonable defaults but allows for fine grained control to exclude spectra with too high or too low signal counts and detection of breakdowns.
The region of interest `[600 1000]` in pixel units is used to crop out the desired spectral region.
!!! note
    All spectra are normalized to unity area in the region of interest.

```julia
preprocessed_signals = preprocess_spectra(signal,nr,[600 1000];errorCriteria=ErrorCriteria(max_signal=64000,min_signal=1000,breakdown_threshold=1000,breakdown_location=200))
showErrorCodes(preprocessed_signals)
# show the preprocessed spectra
plot(preprocessed_signals)
```

[`preprocess_spectra`](@ref) returns an array of [`PreProcSpectrum`](@ref) containing the preprocessed signals, a vector that can be used for residual weighting during the fit, the ErrorCode and the original area under the spectrum.

If desired, a noise model can be used to extract the weights for the residual used during the fit to increase precision and accuracy. Check out [`noiseModel!`](@ref) for possible options:
```julia
noiseModel!(:argon,preprocessed_signals,nr,[600 1000])
```

Now we need to generate the libraries:
```julia
ROI1 = [2240 2460]                              # Define Region of Interest for Pump 1
T = 280:20:2400                                 # Temperature range for which the library is created
FrequencyOffset = 985                           # Frequency Offset between Pump 1 and Pump 2
ROI2 = ROI1 .- FrequencyOffset                  # Define Region of Interest for Pump 2
γc = 0.1                                        # preconvolution kernel. has to be smaller than expected linewidth
ωres = 7e-4                                     # resolution for computation of theoretical susceptibility. has to be equal for all species in a given pumping region.
# Define short handles to functions that generate the desired species' theoretical susceptibilities
# These basically set everything fixed with the exception of temperature:
n2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres)     # N2, Pump 1
co2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres)             # CO2, Pump 2
# generate the libraries
lib1 = generateLibrary(T,γc,n2)    # Create library for pumping region 1 containing N2
lib2 = generateLibrary(T,γc,co2)   # ... and pumping region2 containing CO2
lib2org = deepcopy(lib2)
## hack CO2 model: CARSFT model has severely wrong amplitudes
factor = 1.5
lib2.χM = lib2org.χM * factor^2
lib2.χR = lib2org.χR * factor
lib2.χI = lib2org.χI * factor
```

Now we can start fitting. As proposed in the publication, it makes sense to determine the instrumental parameters on a spectrum recorded at known conditions or at least conditions, that are steady and allow for averaging out the noise. 
```julia
Is = mean(preprocessed_signals)   # mean spectrum, disables noise weighting
```
CARS.jl is able to fit a second order polynomial to the data, given a decent starting solution which can be extracted by plotting the data and a simulation under approximated/guessed conditions:
```julia
display(plot(Is))     # plot the mean experimental spectrum
sim = FitParams(T=1000,X=Dict("N2"=>0.7,"CO2"=>0.05),LineWidth=Dict("P1"=>0.2,"P2"=>.7),Instrumental=1)  
ωE,IE = generateSpectrum!(sim,lib1,lib2);
display(plot(ωE,IE))  # plot the simulated spectrum
```
Now pick two points you would like to closely match pixels to wavenumbers.
It does not have to be perfect, we just need a starting solution for the fit.
Make sure to just use points from a single pumping region.
If your frequency offset is not correctly determined, you will get a bad starting solution.
In addition to the points, you need to pass an example spectrum to the function in order to be able to invert the polynomial when converting wavelength to wavenumber.

```julia
wavenumberpolynomial = getWavenumberPolynomial(Is.I,(303,2302),(211,2328),303)
```

Now we need to provide a startingsolution, including the generated starting solution for the wavenumberpolynomial:

```julia
startingsolution = FitParams(T=800,X=Dict("N2"=>0.7,"CO2"=>0.08),LineWidth=Dict("P1"=>0.2,"P2"=>.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=1);
startingsolution.WavenumberPolynomial = wavenumberpolynomial
```

To specify, which values are to be fitted, we need to provide lower and upper bounds for the desired parameter. CARS.jl interprets the parameter it is supposed to use for fitting by detecting if the lower and upper bounds differ from the starting solution. If the value is the same, it will not be fitted. If it is not, it will be used as a fit parameter.

```julia
startingsolution = FitParams(T=800,X=Dict("N2"=>0.7,"CO2"=>0.08),LineWidth=Dict("P1"=>0.2,"P2"=>.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=1);
startingsolution.WavenumberPolynomial = wavenumberpolynomial
lower=copy(startingsolution);upper=copy(startingsolution);
# define lower and upper bounds
# set bounds fairly widely for first fit
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.02
upper.X["CO2"]=0.5
lower.T = 300
upper.T = 2400
lower.Instrumental = 0.5
upper.Instrumental = 1.5
lower.LineWidth["P1"]=0.2
upper.LineWidth["P1"]=0.5
lower.LineWidth["P2"]=0.2
upper.LineWidth["P2"]=1.6
lower.AdditionalFrequencyOffset = -2
upper.AdditionalFrequencyOffset = 2
startingsolution.WavenumberPolynomial = wavenumberpolynomial
# choosing the bounds of the wavenumber polynomial requires some fiddling in the beginning
# a good starting point is +- 4 cm-1 for the first element and +- 10% for the second value.
# the second order term (last element of the array) should be small in the beginning and only
# increased if the fit hits one of the bounds. start with +- 1e-5
lower.WavenumberPolynomial = [2270, 0.26, -1e-5]
upper.WavenumberPolynomial = [2278, 0.3, 1e-5]
```
Great. Now we have the libraries, a reference spectrum that we would like to fit and specified the starting solution and the fit parameter bounds.
Let's get to it!
CARS.jl provides three mechanisms for fitting:
- [`blackBoxFit`](@ref) uses [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) to do a gradient-free optimization.
- [`gradientFit`](@ref) uses [JuMP](https://jump.dev/JuMP.jl/stable/)+[Ipopt](https://github.com/jump-dev/Ipopt.jl) to do a gradient-based optimization.
- [`combinedFit`](@ref) uses several executions of [`blackBoxFit`](@ref) to provide a starting solution for [`gradientFit`](@ref). This helps avoiding local minima.

The syntax for all of these is quite similar, and differ only in the optional arguments. Check their documentation for details.

```julia
params=CARS.blackBoxFit(Is,startingsolution,lower,upper,lib1,lib2,options=BlackBoxOptions(MaxFuncEvals=5000),showplots=true)
params=CARS.gradientFit(Is,startingsolution,lower,upper,lib1,lib2,options=IpoptOptions(print_level=5,max_wall_time=15),showplots=true);
```

!!! tip
    The key to accurate results when fitting single shots is to get the laser linewidths and frequency offset as accurately as possible from these fits.
    It is possible - however only acceptable when the conditions are well controlled and known - to fix the mole fractions and temperature for this fit to the known values to retrieve the instrumental parameters as close as possible.

As we know gained the relevant knowledge about our CARS instrument, we can specialize the libraries we generated to include the final laser linewidths and frequency offsets. We do so by creating a copy first and then fixing these values in the library:

```julia
lib1n = deepcopy(lib1)
lib2n = deepcopy(lib2)
simplifyLibrary!(params,lib1n,lib2n)
```

Now we can start fitting the single shots and all other experimental data, were these instrumental parameters are applicable to.
We again need to define a starting solution. This time, we have to use the results from the first fit, because the correct linewidths and frequency offset are set in this struct. If using another one, the fit might try to refit these parameters. So just use temperature, mole fractions and wavenumber shift:

```julia
startingsolution = deepcopy(params)
# define upper and lower bounds
lower=deepcopy(startingsolution);upper=deepcopy(startingsolution);
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.0
upper.X["CO2"]=0.5
lower.T = 300
upper.T = 2400
lower.WavenumberShift = -1
upper.WavenumberShift = 1
# throw off the starting solution for good measure
startingsolution.T = 2000
startingsolution.X["N2"] = 0.7
startingsolution.X["CO2"] = 0.
```

Now we can start fitting all spectra. The syntax is exactly the same as before. CARS.jl automatically detects the size of the input and uses multiple dispatch to fit the single shots when an array of `PreProcSpectrum` is passed as the first argument:
```julia
## BLACKBOX
results=blackBoxFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,options=BlackBoxOptions(MaxFuncEvals=10000));
## GRADIENT
results=gradientFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,options=IpoptOptions(print_level=1,tol=1e-6,max_wall_time=2));
## COMBINATION
results=combinedFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,ipoptopt=IpoptOptions(print_level=0,tol=1e-6,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=100));
```

Typically, `combinedFit` will give the best results both in terms of runtime and accurarcy. Note, that ```showplots=true``` adds significant overhead to the runtime per spectrum. In many cases, plotting actually is slower than fitting.