using CARS
using JLD2
using Plots
using BenchmarkTools

@time @load "dp_n2_o2_co2_library_large.jld2"

##
# ok, now let's go
# First, we create a single spectrum:

exp_params = FitParams(T=1200,X=Dict("N2"=>.71,"O2"=>.1,"CO2"=>.1),LineWidth=Dict("P1"=>.2,"P2"=>.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=0.5);
ωE,IE = CARS.generateSpectrum!(exp_params,pump1_library,pump2_library);
# # add some noise
# IE = IE#.+(rand(Float32,length(IE)).-0.5)
# # show the spectrum that we would like to fit
# # limit the experimental spectrum to 2200 to 2600 cm-1.
# # obviously the experimental spectrum has to be narrower than what the fit provides
IE = IE[(ωE.>2250) .* (ωE .< 2600)]
ωE = ωE[(ωE.>2250) .* (ωE .< 2600)]
IE = IE./sum(IE)


plot(ωE,IE)
#ylims!(-10,10)
display(plot!())
##
# final step: make it a Nx2 Matrix:
#
# now we start the fitting.
# we need to provide a starting solution
startingsolution = FitParams(T=2400,X=Dict("N2"=>0.7,"O2"=>0.2,"CO2"=>0.05),LineWidth=Dict("P1"=>0.2,"P2"=>1),AdditionalFrequencyOffset=-2.5,chiB=19.74,Instrumental=0.92);
# we also need upper and lower bounds
# providing upper and lower bounds automagically makes them available for the fit
# in other words, if e.g. lower.T = 300 and upper.T = 2000, temperature will be fitted
# if lower.T == upper.T, it will be fixed. this can be done with every variable.
# to conveniently select the fitting vars, copy startingsolution and modify the fields you want to fit
# note that the following won't work:
# lower = startingsolution
# lower.T = 300
# This will modify the temperature field in lower AND startingsolution because these are mutable structs.
# If you are new to Julia, this might confuse you. Read here: https://docs.julialang.org/en/v1/base/base/#=
# instead, we need to copy them and just modify the values we want to:
lower=copy(startingsolution)
upper=copy(startingsolution)
lower.X["N2"]=0.
upper.X["N2"]=1.0
lower.X["O2"]=0.
upper.X["O2"]=1.0
lower.X["CO2"]=0.
upper.X["CO2"]=1.0
lower.T = 300 
upper.T = 2400
lower.Instrumental = 0.1
upper.Instrumental = 1.5
lower.WavenumberShift = -5
upper.WavenumberShift = 5
lower.AdditionalFrequencyOffset = -10
upper.AdditionalFrequencyOffset = 10
lower.LineWidth["P1"]=0.105
upper.LineWidth["P1"]=1.5
lower.LineWidth["P2"]=0.3
upper.LineWidth["P2"]=1.

# ok, now do the actual fitting!
@time params=CARS.blackBoxFit(ωE,IE,startingsolution,lower,upper,pump1_library,pump2_library;options=BlackBoxOptions(MaxFuncEvals=500,TargetFitness=1e-9))
# @time params=CARS.gradientFit(ωE,IE,startingsolution,lower,upper,pump1_library,pump2_library,options=IpoptOptions(print_level=5,tol=1e-5))
plotSpec(ωE,IE,params)
plot!(ωE,params.Spectra.Fit)
display(plot!(ωE,params.Spectra.Res))
show(params)

#end

##
# ok, we now the linewidth and frequency offset now, simplify the library with this information!
# attention, this overwrites the old libraries. run only once!
pump1_libraryn = deepcopy(pump1_library)
pump2_libraryn = deepcopy(pump2_library)
CARS.simplifyLibrary!(params,pump1_libraryn,pump2_libraryn)

## use the result from blackboxoptim as a starting solution for JuMP
# NOTE!: the library now has the frequencyoffset and linewidth "baked in"
# keeping the values constant instructs the fit not to perform the convolutions
# and interpolations, which makes it much(!) faster.
startingsolution = copy(params)
# get lower and upper bounds. but do not change frequency offset or linewidth. these are fixed now!
lower=copy(startingsolution)
upper=copy(startingsolution)
lower.X["N2"]=0.5
upper.X["N2"]=0.8
lower.X["O2"]=0.
upper.X["O2"]=0.3
lower.X["CO2"]=0.
upper.X["CO2"]=0.3
lower.T = 300 
upper.T = 2400

exp_params = FitParams(T=1800,X=Dict("N2"=>.51,"O2"=>.2,"CO2"=>.1),LineWidth=Dict("P1"=>.2,"P2"=>.7),AdditionalFrequencyOffset=0,chiB=19.74,Instrumental=0.5);
ωE,IE = CARS.generateSpectrum!(exp_params,pump1_library,pump2_library);
# # add some noise
# IE = IE#.+(rand(Float32,length(IE)).-0.5)
# # show the spectrum that we would like to fit
# # limit the experimental spectrum to 2200 to 2600 cm-1.
# # obviously the experimental spectrum has to be narrower than what the fit provides
IE = IE[(ωE.>2250) .* (ωE .< 2600)]
ωE = ωE[(ωE.>2250) .* (ωE .< 2600)]
IE = IE./sum(IE)

@time paramsi=CARS.combinedFit(ωE,IE,startingsolution,lower,upper,pump1_libraryn,pump2_libraryn;ipoptopt=CARS.IpoptOptions(print_level=5,tol=1e-7));

plot(ωE,IE)
plot!(ωE,paramsi.Spectra.Fit)
display(plot!(ωE,paramsi.Spectra.Res))

#println(params)