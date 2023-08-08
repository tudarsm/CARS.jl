# this script generates plots for section 2.2: reconstruction of spectra
using Plots
using CARS
import CO2CARS
import DiaCARS
using BenchmarkTools
using Statistics
## fonts: https://github.com/JuliaPlots/Plots.jl/blob/ec4b7a42e1eff7c4484ade7dfc1e53be175f5356/src/backends/gr.jl#L64
figpath = "C:/Users/greifenstein/Documents/armstrong/Documents/paper/2304_JRS_CARS.jl/figures/"
ωlim = (2252,2548)
ramshift="Raman shift in cm⁻¹"

theme(:ggplot2, linewidth=1,xlabel=ramshift,ylabel="√I",legend=false,dpi=600,xlims=ωlim,fontfamily="sans-serif",size=(600,300))

##
include("gen_lib.jl")

#####
## NOTE ON BENCHMARKS:
## To avoid measuring the compilation time, it is run multiple times in a row and the last iteration is measured

## compare times
lib1 = generateLibrary(T,γc,n2)
lib1 = generateLibrary(T,γc,n2)
lib1 = generateLibrary(T,γc,n2)
lib1time=@timed lib1 = generateLibrary(T,γc,n2)    # Create library for pumping region 1 containing N2
generateLibrary(T,γc,co2,o2)
generateLibrary(T,γc,co2,o2)
generateLibrary(T,γc,co2,o2)
lib2time=@timed lib2 = generateLibrary(T,γc,co2,o2)   # ... and pumping region2 containing CO2

## underlying models
for t in T;n2(t);end
for t in T;n2(t);end
for t in T;n2(t);end
n2time = @timed for t in T    
    n2(t)
end

for t in T;o2(t);end
for t in T;o2(t);end
for t in T;o2(t);end
o2time = @timed for t in T    
    o2(t)
end

for t in T;co2(t);end
for t in T;co2(t);end
for t in T;co2(t);end
co2time = @timed for t in T    
    co2(t)
end

## compare times
println("Time required for generating pump library1: $(lib1time.time) s")
println("Time required for generating pump library2: $(lib2time.time) s")
println("Time required for N2(t): $(n2time.time) s")
println("Time required for O2(t): $(o2time.time) s")
println("Time required for CO2(t): $(co2time.time) s")
println("Time spent in CARS.jl for lib1: $(lib1time.time-n2time.time) s")
println("Time spent in CARS.jl for lib2: $(lib2time.time-o2time.time-co2time.time) s")


## benchmarks for γc variation
γcs = [0.05 0.1 0.2 0.3 0.4 0.5]
lib1s = Array{CARS.Library}(undef,length(γcs))
lib2s = Array{CARS.Library}(undef,length(γcs))
for i in eachindex(γcs)
    lib1s[i] = generateLibrary(T,γcs[i],n2)    # Create library for pumping region 1 containing N2
    lib2s[i] = generateLibrary(T,γcs[i],co2,o2)   # ... and pumping region2 containing CO2
end
##
sim = FitParams(T=1000,X=Dict("N2"=>0.6,"CO2"=>0.1,"O2"=>0.2),LineWidth=Dict("P1"=>0.6,"P2"=>.7),Instrumental=1.)
timings = zeros(length(γcs),10)
p1=plot()
for i in eachindex(γcs)
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    for j in axes(timings,2)
        t = @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
        timings[i,j] = t.time
    end
    p1=plot!(ωE,IE)
end

p2=plot(vec(γcs),mean(timings.*1000,dims=2),ylims=(0,1.75),xlims=(γcs[1],γcs[end]),ylabel="Runtime in ms",xlabel="Width of compression kernel in cm⁻¹",markershape=:circle)

l = @layout [a;b]
@show timings
plot(p1,p2,layout=l)

##
png("benchmark_spectra.png")
meantime=mean(timings.*1000,dims=2)
@show minimum(meantime)
# julia> timings
# 6×10 Matrix{Float64}:
#  0.0014794  0.0014654  0.0014669  0.0014683  0.0014685  0.0015062  0.0016719  0.0016035  0.0014412  0.0014193
#  0.000539   0.0005406  0.0005418  0.0005473  0.0005491  0.0005548  0.0005502  0.0005849  0.0005475  0.0005469
#  0.0002552  0.0002537  0.0002754  0.0002399  0.0002445  0.0002394  0.0002359  0.0002395  0.0002314  0.0002311
#  0.0001632  0.0001595  0.0001657  0.0001629  0.0001602  0.0001673  0.0001656  0.0001625  0.0001655  0.0001687
#  0.000127   0.0001272  0.0001262  0.0001253  0.0001261  0.0001265  0.0001262  0.0001261  0.0001267  0.0001251
#  0.0001121  0.0001106  0.000109   0.0001094  0.0001068  0.0001096  0.000108   0.0001105  0.0001133  0.0001083


## Benchmark of CARSopt variation
γfs = [20 10 5 2 1 0.5]
lib1s = Array{CARS.Library}(undef,length(γfs))
lib2s = Array{CARS.Library}(undef,length(γfs))

for i in eachindex(γfs)
    CARSopt.γf = γfs[i]
    lib1s[i] = generateLibrary(T,0.2,n2)    # Create library for pumping region 1 containing N2
    lib2s[i] = generateLibrary(T,0.2,co2,o2)   # ... and pumping region2 containing CO2
end
##
sim = FitParams(T=1000,X=Dict("N2"=>0.6,"CO2"=>0.1,"O2"=>0.2),LineWidth=Dict("P1"=>0.6,"P2"=>.7),Instrumental=1.)
timings = zeros(length(γfs),10)
ωEx,IEx = generateSpectrum!(sim,lib1s[1],lib2s[1]);
p1=plot(ylabel="Δ√I",ylims=(-10,10))
for i in eachindex(γfs)
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
    for j in axes(timings,2)
        t = @timed ωE,IE = generateSpectrum!(sim,lib1s[i],lib2s[i]);
        timings[i,j] = t.time
    end
    # interpolate on best resolved Array
    IEi = CARS.interp1(ωE,lib1s[i].ωres,IE,ωEx;extrapolation=true)
    # p1=plot!(ωEx,IEi)
    # compute residual and plot it
    res = IEi.-IEx
    @show γfs[i],maximum(res[500:end-500]),maximum(res[500:end-500]./IEx[500:end-500])
    p1=plot!(ωEx,IEi.-IEx)
end

p2=plot(vec(γfs),mean(timings.*1000,dims=2),ylims=(0,1.5),xlims=(γfs[end],γfs[1]),ylabel="Runtime in ms",xlabel="d",markershape=:circle)

l = @layout [a;b]
@show timings
meantime=mean(timings.*1000,dims=2)
@show minimum(meantime),maximum(meantime)
plot(p1,p2,layout=l)

##
png("benchmark_spectra_d.png")

## simulate from scratch
@time ω4,I4=fromScratch(sim,n2,co2,o2;ω=ωE)


