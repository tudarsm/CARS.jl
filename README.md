# CARS.jl
Welcome. This is [CARS.jl](https://github.com/tudarsm/cars.jl).
If you haven't read the publication yet, please do so. It will answer many of the questions that you may have about what this code does and why.

## What is CARS.jl?
CARS.jl is a library-based fitting algorithm for coherent anti-Stokes Raman spectra that features:
- Dual-Pump/Single-Pump evaluation
- High processing speed
- High accuracy
- Automatic determination of wavenumber axis using a 2nd order polynomial fitted to the actual data
- In the libraries, the only tabulated parameter is temperature, thus the size grows (nearly) linearly - not with the power of - the number of species


## What CARS.jl does
- generate final spectra from theoretical susceptibilities supplied by other packages, such as DIACars.jl
- generate libraries of complex susceptibilities using a loss-less compression scheme for fitting
- fit spectra to experimental data
- supply several noise models based on shot-noise, Stokes noise and possibly custom noise-models for weighting the residual in the non-linear least squares fit
- store evaluated data in either JLD2 or MAT format for further data evaluation

## What CARS.jl does not do
- Generate theoretical susceptibilities. This has to be done using other packages such as DIACARS.jl. This modular structure allows to easily include additional species.

## Limitations
- CARS.jl assumes independency of theoretical susceptibilities of the mole fractions of other species in the probe volume. This implies that foreign gas broadening is neglected. This is a common approximation and is e.g. also used by CARSFT
- At present, only fitting at constant pressure are possible. Including pressure in the fit requires introducing pressure variations in the library, because it affects the underlying line shapes.
- Single Pump CARS is not directly supported at the moment. Simply pass the same library twice (as Single-Pump is actually a special case of Dual-Pump) to the fit/simulation to workaround this limitation. This is likely to be resolved in a future version but does not have the highest priority.

## Documentation
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tudarsm.github.io/CARS.jl/stable/)