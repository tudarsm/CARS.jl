# CARS.jl
This is CARS.jl. Welcome. This is heavy work in progress. Stay tuned.

## What is CARS.jl?
CARS.jl is a library-based fitting algorithm for coherent anti-Stokes Raman spectra that features:
- Dual-Pump/Single-Pump evaluation
- High processing speed
- High accuracy
- Automatic determination of wavenumber axis (planned)
- In the libraries, the only tabulated parameter is temperature, thus the size grows only linearly, not with the power of the number of species


## What CARS.jl does
CARS.jl does the following things:
- generate final spectra from theoretical spectra supplied by other packages, such as DIACars.jl
- generate libraries for fitting
- fit spectra to experimental data
- store evaluated data in HDF5 format for compatibility with MATLAB and data management

## What CARS.jl does not do
- Generate theoretical susceptibilities. This has to be done using other packages such as DIACARS.jl.

## How to use it
- Install julia https://julialang.org/downloads/
- Install VSCode and the julia extension
- Clone this repository
- Clone at least one package such as DIACARS.jl in the folder next to it. You should end up with the following folder structure:
```
yourfolder/
├─ cars.jl/
│  ├─ docs/
│  ├─ src/
│  ├─ .../
├─ diacars.jl/
│  ├─ docs/
│  ├─ src/
│  ├─ .../
``` 
- in VSCode, open the folder cars.jl. If there is a popup, trust the authors for all files in this parent folder.
- open a julia REPL (ALT+J,O)
- open setup.jl in the src folder and run it (ALT+ENTER)
- look in the REPL until the precompilation is finished
- look at the examples to learn how to use the code


[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tuda_rsm.pages.rwth-aachen.de/cross-sections/cars/jcars/cars.jl/dev/)
[![Build Status](https://git.rwth-aachen.de/tuda_rsm/cross-sections/cars/CARS.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/tuda_rsm/cross-sections/cars/CARS.jl/pipelines)
[![Coverage](https://git.rwth-aachen.de/tuda_rsm/cross-sections/cars/CARS.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/tuda_rsm/cross-sections/cars/CARS.jl/commits/master)
