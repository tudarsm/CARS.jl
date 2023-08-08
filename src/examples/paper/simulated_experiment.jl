using Plots
using CARS
import CO2CARS
import DiaCARS
using BenchmarkTools
using Statistics
using Distributions

##
ωlim = (2260,2545)
ramshift="Raman shift in cm⁻¹"
theme(:ggplot2, linewidth=1,xlabel=ramshift,ylabel="√I",legend=false,dpi=600,xlims=ωlim,fontfamily="sans-serif",size=(600,300))

##########################################
## GENERATE STANDARD LIBRARY FOR THE PAPER
##########################################
CARSopt.γf = 5
ROI1 = [2240 2560]                              # Define Region of Interest for Pump 1
T = 280:5:2400                                 # Temperature range for which the library is created
FrequencyOffset = 985                           # Frequency Offset between Pump 1 and Pump 2
ROI2 = ROI1 .- FrequencyOffset                  # Define Region of Interest for Pump 2
γc = 0.2                                        # preconvolution kernel. has to be smaller than expected linewidth
ωres = 5e-4                                     # resolution for computation of theoretical susceptibility. has to be equal for all species in a given pumping region.
# Define short handles to functions that generate the desired species' theoretical susceptibilities
# These basically set everything fixed with the exception of temperature:
n2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"N2",1,"I",ROI1,ωres)     # N2, Pump 1
o2(t) = DiaCARS.simulateTheoreticalSusceptibility(t,"O2",1,"I",ROI2,ωres)     # N2, Pump 1
co2(t) = CO2CARS.simulateTheoreticalSusceptibility(t,1,ROI2,ωres)             # CO2, Pump 2
# generate the libraries
lib1 = generateLibrary(T,γc,n2)    # Create library for pumping region 1 containing N2
lib2 = generateLibrary(T,γc,co2,o2)   # ... and pumping region2 containing CO2

##

# base conditions
sim = FitParams(T=1000,X=Dict("N2"=>0.6,"CO2"=>0.1,"O2"=>0.2),LineWidth=Dict("P1"=>0.3,"P2"=>.7),Instrumental=1.,AdditionalFrequencyOffset=-3)
ωEx,IEx = generateSpectrum!(sim,lib1,lib2);
plot(ωEx,IEx,ylims=(0,maximum(IEx)))

# experimental conditions
ω = collect(LinRange(2260,2545,1000))
function randomSpectrum(sim,ω,libs...)
    par = deepcopy(sim)
    par.T = rand(Uniform(300.,2200.),1)[1]
    par.X["N2"] = rand(Uniform(0.6,0.7),1)[1]
    par.X["O2"] = rand(Uniform(0.,0.2),1)[1]
    par.X["CO2"] = rand(Uniform(0.,0.1),1)[1]
    par.WavenumberShift = rand(Uniform(-0.8,0.8),1)[1]
    ωa,Ia = generateSpectrum!(par,lib1,lib2);
    return CARS.interp1(ωa,lib1.ωres,Ia,ω;extrapolation=false), par
end

function randomSpectra(count,sim,ω,libs...)
    spectra = Array{PreProcSpectrum}(undef,count)
    truth = Array{FitParams}(undef,count)
    for i = 1:count
        I,truth[i] = randomSpectrum(sim,ω,libs...)
        spectra[i] = PreProcSpectrum(reverse(vec(I))./sum(I),[1.],sum(I),0)
    end
    return spectra,truth
end

# generate random spectra
spectra,truth = randomSpectra(10,sim,ω,lib1,lib2);
plot(ω,spectra[1].I)

## try to fit the simulated experiment
startingsolution = FitParams(T=800.,X=Dict("N2"=>0.7,"CO2"=>0.1,"O2"=>0.1),LineWidth=Dict("P1"=>0.5,"P2"=>.5),AdditionalFrequencyOffset=0.,Instrumental=.8);
startingsolution.WavenumberPolynomial = [2258, 0.275, 0]
lower=copy(startingsolution);upper=copy(startingsolution);
# define lower and upper bounds
# set bounds fairly widely for first fit
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.02
upper.X["CO2"]=0.5
lower.X["O2"]=0.02
upper.X["O2"]=0.5
lower.T = 300
upper.T = 2300
lower.Instrumental = 0.8
upper.Instrumental = 1.2
lower.LineWidth["P1"]=0.25
upper.LineWidth["P1"]=0.5
lower.LineWidth["P2"]=0.5
upper.LineWidth["P2"]=1.2
lower.AdditionalFrequencyOffset = -5
upper.AdditionalFrequencyOffset = 5
lower.WavenumberPolynomial = [2258, 0.275, 0]
upper.WavenumberPolynomial = [2261, 0.295, 0]

t = CARS.getWavenumberArray(startingsolution.WavenumberPolynomial,1000)
plot(t,spectra[1].I)
##
# params=CARS.blackBoxFit(spectra[1],startingsolution,lower,upper,lib1,lib2,options=BlackBoxOptions(MaxFuncEvals=10000),showplots=true)
##
# params=CARS.gradientFit(spectra[1],startingsolution,lower,upper,lib1,lib2,options=IpoptOptions(print_level=5,max_wall_time=60,tol=1e-6),showplots=true);
params=CARS.combinedFit(spectra[1],startingsolution,lower,upper,lib1,lib2,blackboxopt=BlackBoxOptions(MaxFuncEvals=1000),ipoptopt=IpoptOptions(print_level=5,max_wall_time=60,tol=1e-7),showplots=true);
truth[1]
params

## PLOT THE FIT OF THE REFERENCE SPECTRUM
p1=plot(params.Spectra.omega,params.Spectra.Exp,label="Library",xlims=ωlim,xlabel="",legend=true,xformatter=Returns(""))
p1=plot!(params.Spectra.omega,params.Spectra.Fit,label="Fit",xlims=ωlim,xlabel="",legend=true,xformatter=Returns(""),linestyle=:dash)
p2=plot(params.Spectra.omega,params.Spectra.Res,label="Difference",xlims=ωlim,ylims=(-0.0001,0.0001),yticks=-0.0001:0.0001:0.0001,ylabel="Δ√I",color=:green)
l = @layout [a;b{0.2h}]
plot(p1,p2,layout=l)
savefig("fit_reference.png") 
## simplify the library
# at first, make a copy of our library
CARSopt.γf = 2
lib1n = deepcopy(lib1)
lib2n = deepcopy(lib2)
# and now simplify the libraries. this sets a new preconvolution (and instructs the fit to avoid the second convolution)
# also, pump region2 is reinterpolated to account for the new frequencyoffset
@time simplifyLibrary!(params,lib1n,lib2n)


##
startingsolution = deepcopy(params)
# define upper and lower bounds
lower=deepcopy(startingsolution);upper=deepcopy(startingsolution);
lower.X["N2"]=0.6
upper.X["N2"]=0.75
lower.X["CO2"]=0.0
upper.X["CO2"]=0.5
lower.X["O2"]=0.0
upper.X["O2"]=0.5
lower.T = 300
upper.T = 2400
lower.WavenumberShift=-1
upper.WavenumberShift=1
# throw of the starting solution for good measure
startingsolution.T = 300
startingsolution.X["N2"] = 0.7
startingsolution.X["CO2"] = 0.
startingsolution.X["O2"] = 0.1

# fit all single shot spectra
##

spectra,truth = randomSpectra(1000,sim,ω,lib1,lib2);


## BLACKBOX
# results=blackBoxFit(spectra,startingsolution,lower,upper,lib1n,lib2n;showplots=false,options=BlackBoxOptions(MaxFuncEvals=1000));
# ## GRADIENT
# results=gradientFit(spectra,startingsolution,lower,upper,lib1n,lib2n;showplots=false,options=IpoptOptions(print_level=1,tol=1e-5,max_wall_time=2));

## COMBINATION
results_clean=combinedFit(spectra,startingsolution,lower,upper,lib1n,lib2n;showplots=false,ipoptopt=IpoptOptions(print_level=0,tol=1e-6,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=300));

## with noise
function addNoise(spectra::AbstractArray{PreProcSpectrum})
    noisyspectra=deepcopy(spectra)
    for i in eachindex(spectra)
        noisyspectra[i].I .*= (1 .+rand(Normal(0,0.075),length(noisyspectra[i].I)))
    end
    return noisyspectra
end
spectra_noise = addNoise(spectra)
plot(ω,spectra_noise[1].I)
##
results=blackBoxFit(spectra_noise,startingsolution,lower,upper,lib1n,lib2n;showplots=false,options=BlackBoxOptions(MaxFuncEvals=1000));
##
results=gradientFit(spectra_noise,startingsolution,lower,upper,lib1n,lib2n;showplots=false,options=IpoptOptions(print_level=1,tol=1e-5,max_wall_time=2));
##
results=combinedFit(spectra_noise,startingsolution,lower,upper,lib1n,lib2n;showplots=false,ipoptopt=IpoptOptions(print_level=0,tol=5e-6,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=100));
##
theme(:ggplot2, linewidth=1,xlabel=ramshift,legend=false,dpi=600,fontfamily="sans-serif",size=(600,300))
function scatterResults(results,truth,field)
    trueT = [truth[i].T for i in eachindex(truth)]
    trueN2 = [truth[i].X["N2"] for i in eachindex(truth)]
    trueCO2 = [truth[i].X["CO2"] for i in eachindex(truth)]
    trueO2 = [truth[i].X["O2"] for i in eachindex(truth)]

    fitT = [results[i].T for i in eachindex(results)]
    fitN2 = [results[i].X["N2"] for i in eachindex(results)]
    fitCO2 = [results[i].X["CO2"] for i in eachindex(results)]
    fitO2 = [results[i].X["O2"] for i in eachindex(results)]

    if field == :T
        p=scatter(trueT,fitT,alpha=0.1,color=:black,xlabel="True T in K",ylabel="Fit T in K")
        @show mean(fitT.-trueT),std(fitT.-trueT)
    elseif field == :N2
        p=scatter(trueN2,fitN2,alpha=0.1,color=:black,xlabel="True N₂", ylabel="Fit N₂")
        @show mean(fitN2.-trueN2),std(fitN2.-trueN2)
    elseif field == :CO2
        p=scatter(trueCO2,fitCO2,alpha=0.1,color=:black,xlabel="True CO₂", ylabel="Fit CO₂",ylims=(0,0.2))
        @show mean(fitCO2.-trueCO2),std(fitCO2.-trueCO2)
    elseif field == :O2
        p=scatter(trueO2,fitO2,alpha=0.1,color=:black,xlabel="True O₂", ylabel="Fit O₂")
        @show mean(fitO2.-trueO2),std(fitO2.-trueO2)
    end
    # display(scatter!())
    return p
end

p9 = scatterResults(results,truth,:T)
p10 = scatterResults(results,truth,:N2)
p11 = scatterResults(results,truth,:CO2)
p12 = scatterResults(results,truth,:O2)
plot(p9,p10,p11,p12)
##
savefig("comparison_noise.png")

##
p9 = scatterResults(results_clean,truth,:T)
p10 = scatterResults(results_clean,truth,:N2)
p11 = scatterResults(results_clean,truth,:CO2)
p12 = scatterResults(results_clean,truth,:O2)
plot(p9,p10,p11,p12)

##
savefig("comparison_clean.png")

##
##### GRADIENT TESTING
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
# lower.Instrumental = 0.5
# upper.Instrumental = 1.5
# lower.LineWidth["P1"]=0.25
# upper.LineWidth["P1"]=0.7
# lower.LineWidth["P2"]=0.25
# upper.LineWidth["P2"]=1.6
# lower.AdditionalFrequencyOffset = -5
# upper.AdditionalFrequencyOffset = 5
lower.WavenumberPolynomial = [2258, 0.45, 0]
upper.WavenumberPolynomial = [2263, 0.5, 0]
lower.WavenumberShift =-6
upper.WavenumberShift=1
@benchmark CARS.benchmarkGradient(spectra[1],startingsolution,lower,upper,lib1,lib2)
# @benchmark CARS.benchmarkValue(spectra[1],startingsolution,lower,upper,lib1,lib2)
# @time CARS.testGradient(spectra[1],startingsolution,lower,upper,lib1,lib2);
# params=CARS.gradientFit(spectra[1],startingsolution,lower,upper,lib1,lib2,options=IpoptOptions(print_level=5,max_wall_time=60),showplots=true);
# CARS.benchmarkGradient(spectra[1],startingsolution,lower,upper,lib1,lib2)
## compare accuracy/speed noise spectra
# results=combinedFit(spectra_noise,startingsolution,lower,upper,lib1n,lib2n;showplots=false,ipoptopt=IpoptOptions(print_level=0,tol=1e-5,max_wall_time=2),blackboxopt=BlackBoxOptions(MaxFuncEvals=500));