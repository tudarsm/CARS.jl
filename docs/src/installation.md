```@meta
CurrentModule = CARS
```
# Installation instructions
- Install Julia https://julialang.org/downloads/
- Install VSCode and the Julia extension
- Open VSCode and create a new file with ```.jl``` extension and saven in a folder that makes sense to you
- Open a REPL with ALT+J+O (without releasing ALT)
- add at least CARS.jl and one module to simulate theoretical susceptibilities (such as DiaCARS.jl) with the Package manager. This will install and precompile all necessary packages for you:

```julia
] add https://github.com/tudarsm/cars.jl.git#master
] add https://github.com/tudarsm/diacars.jl.git#master
```

- Have a look at the examples provided in the code or this documentation to get started
- It is planned to add these packages to the Julia registry to simplify the installation

!!! note "Note"
    There is no need to clone the packages, unless you actually want to modify CARS.jl itself.



For more information on using VS Code with Julia, check out https://code.visualstudio.com/docs/languages/julia


# How to check/update the currently installed version?
## Check Version
```julia
] st CARS
```
This will return e.g. the following:
```julia
(CARS) pkg> st CARS
[a83ec7ef] CARS v0.2.4 `https://git.rwth-aachen.de/tuda_rsm/cross-sections/cars/jcars/cars.jl.git#master`
```

## Update Version
```julia
] up CARS
```