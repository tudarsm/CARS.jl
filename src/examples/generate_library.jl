# ok, you planned your experiment and you know your region of interest(s)
# and would like to generate libraries for fitting
# great, let's go over some often used cases:
# 1. Single Pump N2 CARS
# 2. Dual Pump N2+O2 CARS

# For all cases, it is important to understand, which degrees of freedom of the physical model
# are fixed during library generation and which are still available during fitting.
# 1. Temperature is only available at the values you specify at library generation.
#    1 K intervals are more than fine enough, you're experimental uncertainty will be much higher.
# 2. Pressure. Pressure broadening affects the Lorentzian width of each individual transition.
#    As such, this is fixed in the library because only summed up spectra are tabulated, not individual transitions.
# 3. The region of interest has to be broader than your experimental spectrum to allow shifting it a little.
# 4. The lineshape has to be preselected.
#    In principle, one library could include several lineshapes for different spectra,
#    but there is no real reason to do so. Implement it if you need it :)
# 5. The preconvolution kernel for compression should be:
# UPDATE THIS SECTION WHEN FINISHED! Use only one linewidth for compression!
#    a) smaller than the lowest laser linewidth
#    b) large enough to actually compress the library sufficiently.
#    For single pump CARS with Nd:YAG lasers without injection seeding, 0.7 cm-1 is a good value
#    For dual pump CARS, the two pumping regions can have different preconvolution kernels.
#    But: be mindful, which pump laser convolves which spectral region!
#    E.g. when using a 532+565+607 nm combination for N2/O2/CO2,
#    the green beam convolves the O2/CO2 spectrum, the yellow beam convolves the N2 spectrum.
#    This may be confusing at first, but makes perfect sense when you look at the CARS energy scheme (do it!).
#    If you have very good knowledge of the laser linewidths (e.g. from previous fits with the same setup),
#    you can compress the library with these fixed values. Later in the fit, set the starting solution and boundaries to the fixed value.
#    This essentially deactivates two degrees of freedom and omits the costly convolutions when generating a spectrum, making the fit much faster.

#    In principle, even within the same pumping regions, different species could have different
#    Preconvolution Kernels. However, this severely increases complexity of the code and reduces
#    performance. So this is not implemented.

# but now, let's get started.
# include required packages
using CARS
import DiaCARS
using Plots
using JLD2
import CO2CARS


#-----------------------------------
## Single Pump N2-CARS
T=300:100:2400      # temperature range for library generation.
ROI1 = [2200 2700]  # Region of interest in wavenumbers
γc = 0.1            # preconvolution kernel
ωres = 1e-3         # resolution of theoretical susceptibility

# generate the library for N2. This will take a couple of minutes.
n2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres)   # create a function that returns χ depending only(!) temperature. all other parameters are fixed
library = CARS.generateLibrary(T,γc,n2)

# store the library on disk (if you want to - for some libraries it is actually faster to generate them)
@save "single_pump_n2.jld2" library

# load the library
@load "single_pump_n2.jld2"

#-----------------------------------
## Dual Pump N2/O2-CARS
ROI1 = [2200 2700]      # Pump region 1
FrequencyOffset = 985;  # Offset for pump region two. Pump1-FrequencyOffset=Region2
ROI2 = ROI1 .- FrequencyOffset  # Pump region 2
T=280:20:2400 
γc = 0.1            # preconvolution kernel
ωres = 1e-3         # resolution of theoretical susceptibility

# Just as above, construct functions that create χ only as a function of temperature!
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres) # N2, Pump 1
#o2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI1,ωres) # O2, Pump 1
#n2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI2,ωres) # N2, Pump 2
o2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres) # O2, Pump 2
# Generate the libraries for both pumping regions
#pump1_library = CARS.generateLibrary(T,γc,n2_pump1,o2_pump1)    # This will fail, because O2 has no transitions in pump region 1
pump1_library = CARS.generateLibrary(T,γc,n2_pump1)             # This will work
#pump2_library = CARS.generateLibrary(T,γc,n2_pump2,o2_pump2)    # This will work, because N2 actually has transitions in that region of interest
pump2_library = CARS.generateLibrary(T,γc,o2_pump2)             # ... but these are weak, so ignore them for now.

# That's it. Store for later use:
@time @save "dp_n2_o2_library_large.jld2" pump1_library pump2_library

#-----------------------------------------------
## Now generate a spectrum from the newly generated libraries!
params = FitParams(T=2100,X=Dict("N2"=>0.7,"O2"=>0.2),LineWidth=Dict("P1"=>1.,"P2"=>2.),AdditionalFrequencyOffset=0,chiB=20,Instrumental=1);
@time ω,I,χM,χI,χR = CARS.generateSpectrum!(params,lib1,lib2)
#params = FitParams(T=2000,X=Dict("N2"=>0.7),LineWidth=[1],Instrumental=1);
#@time ω,I,χM,χI,χR = CARS.generateSpectrum!(params,lib1)
# once preallocated, you can also use this function to reuse the memory:
# the ! is a convention to notify the user that a function overwrites some of its arguments
# in this case, the result is stored in I.
@time CARS.generateSpectrum!(ω,I,χM,χI,χR,params,lib1,lib2)
plot(ω,I)

#-----------------------------------------------
## DP CARS CO2 N2 O2
ROI1 = [2200 2700]      # Pump region 1
FrequencyOffset = 985;  # Offset for pump region two. Pump1-FrequencyOffset=Region2
ROI2 = ROI1 .- FrequencyOffset  # Pump region 2
T=280:20:2400 
γc = 0.1            # preconvolution kernel
ωres = 1e-3         # resolution of theoretical susceptibility
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres) # N2, Pump 1
o2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres) # O2, Pump 2
co2_pump2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres) # O2, Pump 2

pump1_library = CARS.generateLibrary(T,γc,n2_pump1)             # Create library for pumping region 1 containing N2
pump2_library = CARS.generateLibrary(T,γc,o2_pump2,co2_pump2)   # ... and pumping region2 containing O2 and CO2

@time @save "dp_n2_o2_co2_library_large.jld2" pump1_library pump2_library

##
ROI1 = [2250 2435]
T = 280:20:2400
FrequencyOffset = 960
ROI2 = ROI1 .- FrequencyOffset  # Pump region 2
γc = 0.1            # preconvolution kernel
ωres = 6e-4         # resolution of theoretical susceptibility
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres) # N2, Pump 1
co2_pump2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres) # O2, Pump 2

lib1 = CARS.generateLibrary(T,γc,n2_pump1)             # Create library for pumping region 1 containing N2
#lib2 = CARS.generateLibrary(T,γc,co2_pump2)   # ... and pumping region2 containing O2 and CO2

@show sum(isnan.(lib1.χM))