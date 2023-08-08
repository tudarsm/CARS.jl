##########################################
## GENERATE STANDARD LIBRARY FOR THE PAPER
##########################################
ROI1 = [2240 2560]                              # Define Region of Interest for Pump 1
T = 280:20:2400                                 # Temperature range for which the library is created
FrequencyOffset = 985                           # Frequency Offset between Pump 1 and Pump 2
ROI2 = ROI1 .- FrequencyOffset                  # Define Region of Interest for Pump 2
γc = 0.1                                        # preconvolution kernel. has to be smaller than expected linewidth
ωres = 5e-4                                     # resolution for computation of theoretical susceptibility. has to be equal for all species in a given pumping region.
# Define short handles to functions that generate the desired species' theoretical susceptibilities
# These basically set everything fixed with the exception of temperature:
n2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres)     # N2, Pump 1
o2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres)     # N2, Pump 1
co2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres)             # CO2, Pump 2
# generate the libraries
lib1 = generateLibrary(T,γc,n2)    # Create library for pumping region 1 containing N2
lib2 = generateLibrary(T,γc,co2,o2)   # ... and pumping region2 containing CO2

## FUNCTIONS FOR CONVENIENCE
∗(a,b) = CARS.dualConv(a,b)
function fromScratch(sim,n2,co2,o2;ω=nothing)
    # get n2, co2 and o2 susceptibilities at temperature T
    n2_χR,n2_χI,n2_ω,n2_species,n2_molecularParameters = n2(sim.T)
    co2_χR,co2_χI,co2_ω,co2_species,co2_molecularParameters = co2(sim.T)
    o2_χR,o2_χI,o2_ω,o2_species,o2_molecularParameters = o2(sim.T)

    # get the non-resonant signal
    N = CARS.NA*273.15/sim.T*(1/CARS.molvol);          # number density based on temperature and pressure
    @show XB = 1-sim.X["N2"]-sim.X["O2"]-sim.X["CO2"]        # buffer gas mole fraction
    chiB = sim.chiB/CARS.NA*CARS.molvol                     # buffer gas nr susc per molecule
    @show χNR = N*(n2_molecularParameters.CHINR * sim.X["N2"] +
                co2_molecularParameters.CHINR * sim.X["CO2"] + 
                o2_molecularParameters.CHINR * sim.X["O2"] + 
                chiB * XB)                                  # ... χNR already multiplied with N (just like spectra)

    # convolution kernels
    g2 = CARS.gaussiankernelT(ωres,sim.LineWidth["P1"])                
    g1 = CARS.gaussiankernelT(ωres,sim.LineWidth["P2"])                
    gI = CARS.gaussiankernelT(ωres,sim.Instrumental)        
    
    # χ for pumping regions
    χ1R = sim.X["N2"]*n2_χR.+χNR
    χ1I = sim.X["N2"]*n2_χI
    χ2R = sim.X["O2"]*o2_χR+sim.X["CO2"]*co2_χR.+χNR
    χ2I = sim.X["O2"]*o2_χI+sim.X["CO2"]*co2_χI
    
    χ1M = @. χ1R^2 + χ1I^2
    χ2M = @. χ2R^2 + χ2I^2
    # CARS signal according to Eq. 6 but without refactoring χNR
    I₄ = χ1M∗g2 + χ2M∗g1 +
         2*(χ1R∗g2).*(χ2R∗g1)+
         2*(χ1I∗g2).*(χ2I∗g1)
    


    # convolution with apparatus
    I₄ = I₄∗gI

    # interpolate on the library grid for comparison
    if ω == nothing
        ωexp = collect(n2_ω[1]:lib1.ωres:n2_ω[end])
    else
        ωexp = ω
    end
    I₄ = CARS.interp1(n2_ω,ωres,I₄,ωexp;extrapolation=false)

    @. I₄ = sign(I₄)*sqrt(abs(I₄))
    return ωexp,I₄

end