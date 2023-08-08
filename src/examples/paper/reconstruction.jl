# this script generates plots for section 2.2: reconstruction of spectra
using Plots
using CARS
import CO2CARS
import DiaCARS

## fonts: https://github.com/JuliaPlots/Plots.jl/blob/ec4b7a42e1eff7c4484ade7dfc1e53be175f5356/src/backends/gr.jl#L64
figpath = "C:/Users/greifenstein/Documents/armstrong/Documents/paper/2304_JRS_CARS.jl/figures/"
ωlim = (2252,2548)
ramshift="Raman shift in cm⁻¹"

theme(:ggplot2, linewidth=1,xlabel=ramshift,ylabel="√I",legend=false,dpi=600,xlims=ωlim,fontfamily="sans-serif",size=(600,300))

##
include("gen_lib.jl")

## spectrum generation with library
sim = FitParams(T=1600,X=Dict("N2"=>0.6,"CO2"=>0.1,"O2"=>0.2),LineWidth=Dict("P1"=>0.2,"P2"=>.7),Instrumental=1.)
@time ωE,IE = generateSpectrum!(sim,lib1,lib2);
p1=plot(ωE,IE,label="Library",xlims=ωlim,xlabel="",legend=true,xformatter=Returns(""))

# spectrum generation from scratch
@time ω4,I4=fromScratch(sim,n2,co2,o2;ω=ωE)
p1=plot!(ω4,I4,linestyle=:dash,label="Direct")
# residual
p2=plot(ω4,IE.-I4,label="Difference",xlims=ωlim,ylims=(-0.05,0.05),yticks=-0.05:0.05:0.05,ylabel="Δ√I",color=:green)

l = @layout [a;b{0.2h}]
plot(p1,p2,layout=l)

# get the maximum deviation but without edge effects
@show maximum(IE[500:end-500].-I4[500:end-500]), maximum(I4), maximum(IE[500:end-500].-I4[500:end-500])/maximum(I4)

savefig(figpath.*"reconstruction.png") 