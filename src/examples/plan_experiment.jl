# This may be helpful if you would like to plan your experiment
# You can enter the laser wavelengths, bandwidth and such and see, what types of spectra to expect and which ranges to include in the library generation

# include required packages
import DiaCARS
using Plots

#
# Enter your wavelengths here
λ1 = 532e-9;    # pump 1
λ2 = 565e-9;    # pump 2
λ3 = 611e-9;    # stokes
Δλ3 = 6e-9;    # stokes bandwidth FWHM

# Enter your desired species here:
species = ["N2" "O2"]
# Enter the highest temperature you expect here
T = 2400

# pumping regions
ROI1 = @. round([1/λ1-1/(λ3-Δλ3) 1/λ1-1/(λ3+Δλ3)]/100) # in wavenumbers
ROI2 = @. round([1/λ2-1/(λ3-Δλ3) 1/λ2-1/(λ3+Δλ3)]/100) # in wavenumbers

println("------------------------")
println("Pumping region1: $(ROI1)")
println("Pumping region2: $(ROI2)")
println("------------------------")
println("Trying to simulate spectra for these ranges and species...")

# create a plot
plot();

for mol in species
    # ROI 1
    print("$(mol), ROI1: ")
    # Get the theoretical Susceptbility
    try
        χR,χI,ω = DiaCARS.simulateTheoreticalSusceptibility(T,mol,1,"I",ROI1,1e-3);
        # Get the resonant part and plot it
        χRes = @. sqrt(χR^2+χI^2)
        display(plot!(ω,χRes,label="$(mol), ROI1")) # plot! adds to the existing plot
        println("Found transitions")
    catch e
        println(e)
    end
    try
        print("$(mol), ROI2: ")
        χR,χI,ω = DiaCARS.simulateTheoreticalSusceptibility(T,mol,1,"I",ROI2,1e-3);
        # Get the resonant part and plot it
        χRes = @. sqrt(χR^2+χI^2)
        display(plot!(ω,χRes,label="$(mol), ROI2")) # plot! adds to the existing plot
        println("Found transitions")
    catch e
        println(e)
    end
    
end

println("Look at the plots and decide which transitions are relevant for you.")
println("E.g. O2 has no transitions between 2150 and 2700 cm-1, where the N2 transitions are. So this range does not have to be accounted for in the library.")
println("Similarly, N2 has only very weak transitions between 1050 and 1800 cm-1, so these may probably be ignored.")