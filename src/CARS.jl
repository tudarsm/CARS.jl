module CARS

using DSP
using Interpolations
using ProgressBars
using ProgressMeter
using JLD2
using MAT
using Plots
using JuMP
using Ipopt
using Base.Cartesian
using BlackBoxOptim
using ForwardDiff
using LoopVectorization
using BenchmarkTools
using Statistics
using FileWatching

include("includes/interpolations.jl")
include("includes/structs.jl")
include("includes/library_generation.jl")
include("includes/spectrum_generation.jl")
include("includes/fit_spectra.jl")
include("includes/utilities.jl")
include("includes/convolutions.jl")
include("includes/load_spectra.jl")
include("includes/noise_model.jl")
include("includes/preprocess_spectra.jl")
include("includes/plots.jl")
include("includes/live_fit.jl")

# exported functions
# library stuff
export generateLibrary
export assembleLibraries
export simplifyLibrary!

# structs
export FitParams
export BlackBoxOptions
export IpoptOptions
export ErrorCriteria
export RawSpectra
export PreProcSpectrum

# fitting
export gradientFit
export blackBoxFit
export combinedFit
export quasiLiveFit

# utilities
export dualConv
export plotSpec
export plot         # for convenience, export the plot function
export generateSpectrum!
export @save
export @load
export getWavenumberPolynomial
export mean
export save
export load

# handling of spectra
export load_spectra
export preprocess_spectra
export noiseModel!
export showErrorCodes

# some constants
const NA = 6.02214086e23  # Avogadro
const molvol = 2.2413e4 # ideal gas molar volume at 101325 Pa and 273.15 K in L/kmol

# instantiate CARSOptions as global variable
CARSopt = CARSOptions();
# export CARSopt defaults as a global variable
export CARSopt
# also export some interpolators:
export Constant
export Linear
export Quadratic
export Cubic

end
