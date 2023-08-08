# a simple example on how to simulate various different spectra
# similar to matlab, ## denotes code cells. enter a cell and hit alt+enter to run it
# to run a single line, press ctrl+enter
# first run may take a while.

## include required packages
using CARS
import DiaCARS
import CO2CARS
using Plots

# some constants
ωres = 1e-3     # wavenumber resolution for calculating susceptibilities

## N2 Theoeretical Susceptbility
# Spectrum in the region 2200-2390 cm-1 using the isolated lines model at 1 bar
χR,χI,ω = DiaCARS.simulateTheoreticalSusceptibility(2000,"N2",1,"I",[2200 2390],ωres);
# Plot it
plot(ω,[χR,χI],title="Complex Susc. of N2 Spectrum at 2000 K",label=["Real" "Imaginary"])

## O2 Theoeretical Susceptbility
# Spectrum in the region 1300-1500 and using the isolated lines model at 1 bar
χR,χI,ω = DiaCARS.simulateTheoreticalSusceptibility(2000,"O2",1,"I",[1300 1800],ωres);
# Plot it
plot(ω,[χR,χI],title="Complex Susc. of O2 Spectrum at 2000 K",label=["Real" "Imaginary"])

## Create a single pump N2 Spectrum from scratch
# get the susceptibility
T = 2000    # Temperature
X = 0.6     # Mole Fraction
P = 1       # Pressure
ROI = [2200 2700]
N = CARS.N(T,P) # Number density
γ = 0.7     # Laser linewidth in cm-1
γi = 1.      # instrumental
# simulate χ
χR,χI,ω,species,MolPar = DiaCARS.simulateTheoreticalSusceptibility(T,"N2",P,"I",ROI,ωres);
# buffer gas Susceptbility
χNRBuffer = CARS.permolecule(8.5)
χNR = N*((1-X)*χNRBuffer + X*MolPar.CHINR)
# scale χ with mole fraction and number density. add purely real non-resonant signal
χR = X*χR .+ χNR
χI = X*χI
# laser linewidth kernel
g = CARS.gaussiankernelT(ωres,γ)
gi = CARS.gaussiankernelT(ωres,γi)
# CARS signal: magnitude squared, then convolved + magnitude convolved, then squared
I = 2*CARS.conv((χR.^2 + χI.^2),g) + 2*(CARS.conv(χR,g)).^2 + 2*(CARS.conv(χI,g)).^2
I = CARS.conv(I,gi)
# interpolate on a coarser grid
ωc = collect(ω[1]:γ/5:ω[end])
I = CARS.interp1(ω,ωres,I,ωc;extrapolation=false)
I = sqrt.(I./sum(I))
plot(ωc,I)


## Create a dual-pump O2-N2 Spectrum from scratch
# get the susceptibility
T = 2000    # Temperature
XN2 = 0.7     # Mole Fraction
XO2 = 0.3     # Mole Fraction
P = 1       # Pressure
ROI = [2200 2700]
ROI2 = ROI1 .- 985
N = CARS.N(T,P) # Number density
γ1 = 0.2     # Laser linewidth in cm-1
γ2 = 0.7     # Laser linewidth in cm-1
γi = 1.      # instrumental
# simulate χ
χRN2,χIN2,ω,species,MolParN2 = DiaCARS.simulateTheoreticalSusceptibility(T,"N2",P,"I",ROI,ωres);
χRO2,χIO2,ωO2,species,MolParO2 = DiaCARS.simulateTheoreticalSusceptibility(T,"O2",P,"I",ROI2,ωres);
# buffer gas Susceptbility
χNRBuffer = CARS.permolecule(8.5)
χNR = N*((1-XN2-XO2)*χNRBuffer + XN2*MolPar.CHINR + XO2*MolParO2.CHINR)
# scale χ with mole fraction and number density. add purely real non-resonant signal
χR1 = XN2*χRN2 .+ χNR
χI1 = XN2*χIN2
χR2 = XO2*χRO2 .+ χNR
χI2 = XO2*χIO2

# laser linewidth kernel
g1 = CARS.gaussiankernelT(ωres,γ1)
g2 = CARS.gaussiankernelT(ωres,γ2)
gi = CARS.gaussiankernelT(ωres,γi)
# CARS signal: magnitude squared, then convolved + magnitude convolved, then squared
I = CARS.conv((χR1.^2 + χI1.^2),g1)+ CARS.conv((χR2.^2 + χI2.^2),g2)+ 2*(CARS.conv(χR1,g1).*CARS.conv(χR2,g2))+ 2*(CARS.conv(χI2,g2).*CARS.conv(χI2,g2))
I = CARS.conv(I,gi)
# interpolate on a coarser grid
ωc = collect(ω[1]:γ/5:ω[end])
I = CARS.interp1(ω,ωres,I,ωc;extrapolation=false)
I = sqrt.(I./sum(I))
plot(ωc,I)

## Create a dual-pump CO2-O2-N2 Spectrum from scratch
# get the susceptibility
ωres = 1e-3
T = 1000    # Temperature
XN2 = 0.75     # Mole Fraction
XO2 = 0.     # Mole Fraction
XCO2 = 0.1     # Mole Fraction
P = 1       # Pressure
ROI = [2200 2700]
ROI2 = ROI .- 985
N = CARS.N(T,P) # Number density
γ1 = 0.2     # Laser linewidth in cm-1
γ2 = 0.7     # Laser linewidth in cm-1
γi = 0.1      # instrumental
# simulate χ
χRN2,χIN2,ω,species,MolParN2 = DiaCARS.simulateTheoreticalSusceptibility(T,"N2",P,"I",ROI,ωres);
χRO2,χIO2,ωO2,species,MolParO2 = DiaCARS.simulateTheoreticalSusceptibility(T,"O2",P,"I",ROI2,ωres);
χRCO2,χICO2,ωCO2,species,MolParCO2 = CO2CARS.simulateTheoreticalSusceptibility(T,P,ROI2,ωres);
# buffer gas Susceptbility
χNRBuffer = CARS.permolecule(18.5)
χNR = N*((1-XN2-XO2-XCO2)*χNRBuffer + XN2*MolParN2.CHINR + XO2*MolParO2.CHINR+XCO2*MolParCO2.CHINR)
# scale χ with mole fraction and number density. add purely real non-resonant signal
χR1 = XN2*χRN2 .+ χNR
χI1 = XN2*χIN2
χR2 = XCO2*χRCO2 + XO2*χRO2 .+ χNR
χI2 = XCO2*χICO2 + XO2*χIO2
# laser linewidth kernel
g1 = CARS.gaussiankernelT(ωres,γ1)
g2 = CARS.gaussiankernelT(ωres,γ2)
gi = CARS.gaussiankernelT(ωres,γi)
# CARS signal: magnitude squared, then convolved + magnitude convolved, then squared
I = CARS.conv((χR1.^2 + χI1.^2),g1)+ CARS.conv((χR2.^2 + χI2.^2),g2) + 2*(CARS.conv(χR1,g1).*CARS.conv(χR2,g2))+ 2*(CARS.conv(χI1,g1).*CARS.conv(χI2,g2))
I = CARS.conv(I,gi)
# interpolate on a coarser grid
ωc = collect(2200:γi/5:2600)
I = CARS.interp1(ω,ωres,I,ωc;extrapolation=false)
I = sqrt.(I)
plot(ωc,I)

# Now do the same with the library.
T = [1000]
n2_pump1(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI,ωres) # N2, Pump 1
o2_pump2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres) # O2, Pump 2
co2_pump2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres) # CO2, Pump 2
γc = 0.1
pump1_library = CARS.generateLibrary(T,γc,n2_pump1)             # Create library for pumping region 1 containing N2
pump2_library = CARS.generateLibrary(T,γc,o2_pump2,co2_pump2)   # ... and pumping region2 containing O2 and CO2

exp_params = FitParams(T=1000,X=Dict("N2"=>XN2,"O2"=>XO2,"CO2"=>XCO2),LineWidth=Dict("P1"=>γ1,"P2"=>γ2),AdditionalFrequencyOffset=0,chiB=18.5,Instrumental=γi);
ωE,IE = CARS.generateSpectrum!(exp_params,pump1_library,pump2_library)

IE = CARS.interp1(ωE,pump1_library.ωres,IE,ωc;extrapolation=false)

plot(ωc,I)
plot!(ωc,IE)
plot!(ωc,I-IE)
xlims!(2370,2400)

# plot(ω,χICO2)
# plot!(pump1_library.ω,pump2_library.χI[:,1,2])

##
gc = CARS.gaussiankernelT(0.001f0,0.1f0)
plot(ω,CARS.conv(χICO2.^2 + χRCO2.^2,gc))
plot!(pump1_library.ω,pump2_library.χM[:,1,2])


##
sig = rand(Float32,30000)
sig = χRCO2
##
γ1 = .1f0
γ = .3f0
γ2 = sqrt(γ^2-γ1^2)
ωres = 0.001f0
g1 = CARS.gaussiankernelT(ωres,γ1)
g = CARS.gaussiankernelT(ωres,γ)
g2 = CARS.gaussiankernelT(ωres,γ2)

using DSP
a = CARS.conv(sig,g)
b = CARS.conv(CARS.conv(sig,g1),g2)

plot(a[10000:end-10000])
plot(a[10000:end-10000]-b[10000:end-10000])
