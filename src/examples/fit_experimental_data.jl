## this is an example workflow to fit experimental spectra

# at first, import the required packages
using CARS
using Plots
import DiaCARS
import CO2CARS

##########################################
## LOAD EXAMPLE DATA
##########################################
# check out load_spectra for hints what to do if you do not use WinSpec's SPE format.
signal = load_spectra("./src/examples/data/210622A163.SPE")     # example file for measurement
nr = load_spectra("./src/examples/data/210622A144.SPE")         # example non-resonant signal acquired in Argon

# the data returned is a struct with fields I, Backround and SingleShotNR. The latter is not yet implemented but can be used if the non-resonant signal is mesured simultaneously with the CARS signal
@show fieldnames(CARS.RawSpectra)
# If you plan to do so, check comments in load_spectra.jl on some assistance for what to do.
# Very likely, another fit parameter for a horizontal shift of the measured NR signal is required though, because beam-steering might be different. As there is yet no data available where this could be tested, this is not yet implemented. Adding this should be straight forward by omitting the normalization to the non resonant background (and area) in the preprocessing step and instead doing it during the fit.
# A new fit parameter is simply added by adding it to the FitParams struct in structs.jl. This should make it available for the fit automagically.
# A starting solution for the horizontal fit could be obtained by measuring a non-resonant signal with the main CARS setup and with the reference beams simultaneously and horizontally shifting the signals until they fit. then the fit just has to wiggle it around by a couple of pixels to adjust for single-shot variations.

# show the spectra
l = @layout [a b;c d]
p1=plot(signal.I)
p2=plot(nr.I)
p3=plot(signal.BG)
p4=plot(nr.BG)
plot(p1,p2,p3,p4,layout=l,legend=:none)

# NOTE: Possibly you received a warning about saturated samples in the background signal.
# this is because the CCD camera might have accumulated charge on it for the first acquisition
# the default background extractor detects this and removes it from the data before computing the 
# average background. look at the plots to see if the backgrounds are correct.
# the NR signal in the upper right corner may still have it, but this will be checked during preprocessing

##########################################
## PREPROCESSING
##########################################
# preprocessing involves background correction, normalization to the non-resonant Argon signal
# in addition, errorcodes are returned for invalid spectra (saturation, signal too low, breakdowns,...)
# check out preprocess_spectra.jl for definition of thresholds for error checking
# samples are not removed in order to not mess up the order of samples, but non-zero errorcodes will be skipped while fitting
# to use the defaults, just run preprocess_spectra like this:
preprocessed_signals = preprocess_spectra(signal,nr,[600 1000])
showErrorCodes(preprocessed_signals)
# if you want to adjust the thresholds, supply a ErrorCriteria struct to preprocess_spectra() like this:
preprocessed_signals = preprocess_spectra(signal,nr,[600 1000];errorCriteria=ErrorCriteria(max_signal=64000,min_signal=1000,breakdown_threshold=1000,breakdown_location=200))
showErrorCodes(preprocessed_signals)

# show the preprocessed spectra
plot(preprocessed_signals)

# NOTE:
# The preprocessing is extremely case dependent. If you change anything about the current workflow,
# it is likely required to implement a new preprocessing scheme.
# The only requirement for the fitting procedure is that the output of your preprocessing function
# returns a PreProcSpectrum struct compatible with the one above.
# N: Number of spectra, M: spectral direction
# A: 1xN
# I: MxN
# ErrorCode: Vector with N elements
# w: not added at this point. Optionally added by the noise model code below.

##########################################
## NOISE MODEL EXTRACTION (OPTIONAL)
##########################################
@time noiseModel!(:argon,preprocessed_signals,nr,[600 1000])

##########################################
## LIBRARY GENERATION
##########################################
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


##########################################
## FIRST FIT
##########################################
# at first, we need a good starting solution for the wavenumber array. for this, it may be useful to look at the mean spectrum of your data
plotlyjs()
Is = mean(preprocessed_signals)     # this is a dedicated function to extract the mean spectrum of all preprocessed_signals which already ignores bad samples for you (e.g. saturated etc, check preprocess_spectra.jl for details)
display(plot(Is))
# .. and compare it to a spectrum at similar conditions
sim = FitParams(T=1000,X=Dict("N2"=>0.7,"CO2"=>0.05),LineWidth=Dict("P1"=>0.2,"P2"=>.7),Instrumental=1)
ωE,IE = generateSpectrum!(sim,lib1,lib2);
display(plot(ωE,IE))
# now pick two points you would like to closely match pixels to wavenumbers
# it does not have to be perfect, we just need a starting solution for the fit.
# make sure to just use points from a single pumping region
# if your frequency offset is not correctly determined, you will get a bad starting solution
# In addition to the points, you need to pass an example spectrum to the function in order
# to be able to invert the polynomial when converting wavelength to wavenumber
wavenumberpolynomial = getWavenumberPolynomial(Is.I,(303,2302),(211,2328),303)
##
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
# if your starting solution is bad, generate some runs of blackboxoptim below to get a better starting value
# and repeat the process until it converges. the gradient based fit has trouble fitting the wavenumber polynomial
# when the bounds are too wide.
# also note: do not fit the wavenumber shift at this point! this is essentially the same degree of freedom as the first element of the wavenumber polynomial.
# the intended use for the wavenumber shift is for correcting shot-to-shot fluctuations when doing single shot fitting

# either use blackboxfitting or gradient based fitting to fit the initial spectra
## BlackBoxOptim
params=CARS.blackBoxFit(Is,startingsolution,lower,upper,lib1,lib2,options=BlackBoxOptions(MaxFuncEvals=5000),showplots=true)
## Gradient
params=CARS.gradientFit(Is,startingsolution,lower,upper,lib1,lib2,options=IpoptOptions(print_level=5,max_wall_time=15),showplots=true);
## you can also run a blackboxfit to generate a good starting solution for the gradient-based fit by passing the result
# of blackboxoptim to gradientFit as a startingsolution:
params=CARS.blackBoxFit(Is,startingsolution,lower,upper,lib1,lib2,options=BlackBoxOptions(MaxFuncEvals=1000),showplots=true);
params=CARS.gradientFit(Is,params,lower,upper,lib1,lib2,options=IpoptOptions(print_level=5,max_wall_time=60),showplots=true);


##########################################
## SIMPLIFY LIBRARY
##########################################
# From the first fit, we now "know" the laser linewidths, wavenumber axis, frequency offset
# Assuming these are constant, we can incorporate them in a simplified version of the library to make the fit much faster
# We are then left with mole fractions, temperature and wavenumber shift as fit parameters

# at first, make a copy of our library
lib1n = deepcopy(lib1)
lib2n = deepcopy(lib2)
# and now simplify the libraries. this sets a new preconvolution (and instructs the fit to avoid the second convolution)
# also, pump region2 is reinterpolated to account for the new frequencyoffset
simplifyLibrary!(params,lib1n,lib2n)

# if you want to store it for later use, execute the following line.
# do not forget to store the params struct alongside the libraries as it contains useful information!
# note that the non-simplified libraries are not that critical, because they can be recreated any time
# but the simplified library depends on the quality of your initial fit which may depend on your preprocessing,
# bounds for the fit and many more. it is probably a good idea to store it.
# @save "library.jld2" lib1n lib2n params

##########################################
## FIT THE SAME SPECTRUM WITH THE SIMPLIFIED LIBRARY
##########################################
# use the previous result as starting solution
startingsolution = deepcopy(params)
# define upper and lower bounds
lower=deepcopy(startingsolution);upper=deepcopy(startingsolution);
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.0
upper.X["CO2"]=0.5
lower.T = 300
upper.T = 2400
# throw of the starting solution for good measure
startingsolution.T = 2000
startingsolution.X["N2"] = 0.7
startingsolution.X["CO2"] = 0.

# fit all single shot spectra

## BLACKBOX
results=blackBoxFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,options=BlackBoxOptions(MaxFuncEvals=10000));
## GRADIENT
results=gradientFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,options=IpoptOptions(print_level=1,tol=1e-6,max_wall_time=2));
results=gradientFit(preprocessed_signals[6],startingsolution,lower,upper,lib1n,lib2n;showplots=true,options=IpoptOptions(print_level=5,tol=1e-6,max_wall_time=2));

## COMBINATION
results=combinedFit(preprocessed_signals,startingsolution,lower,upper,lib1n,lib2n;showplots=true,ipoptopt=IpoptOptions(print_level=0,tol=1e-6,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=100));

##########################################
## FIT ALL FILES IN A GIVEN FOLDER AND STORE RESULTS
##########################################

##########################################
## ONLINE SUPERVISION OF FOLDER (OPTIONAL)
##########################################
# currently, there is only the option to do a quasi live fit
# for this case, the code waits for a new file in the requested folder
# and starts processing immediately after it is available
# you can choose whether to just use the newest file
# or all files in the folder. beware, that the latter than might lag behind by more than one measurement
# if choosing to skip data in between, data processing is only done for the newest file (that was not processed yet)

# to do the live fitting, you need to construct some helper functions for loading, preprocessing and fitting
# these functions are passed as arguments to the live fit function to instruct it on how to do these operations
# this approach makes it very easily extendable to other file formats and other types of preprocessing

# choose function to load data
load(file) = load_spectra(file)         
# function to preprocess data. nr is the one from above and ROI is defined to crop the data
preprocess(signal) = preprocess_spectra(signal,nr,[600 1000]) 
# function for the fit. only used compressed libraries.
# don't use showplots because livefitting has it's own fancy plot!
# also, plotting reduces efficiency of the fitting procedure.
fit(preprocessed) = combinedFit(preprocessed,startingsolution,lower,upper,lib1n,lib2n;showplots=false,ipoptopt=IpoptOptions(print_level=1,tol=1e-6,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=100));

# start watching the current folder
# to test it, copy files into the folder and see the fit results immediately afterwards
# you can specify a skip factor to be faster by skipping spectra
# e.g. a skipping factor of 2 only evaluates half the spectra in the file
# if you choose the option newest_only = true, only the most recent file (after finishing processing of another file) will be used and files stored during processing will be ignored
quasiLiveFit(pwd(),load,preprocess,fit;skip=10,newest_only=false)

# warning: if this is running for a long time, you might want to delete the plots in vscode
# otherwise, you might run out of memory ;)
# to do so, press shift+ctrl+p and type delete plots