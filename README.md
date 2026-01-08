# PDEBase.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/PDEBase/stable/)

[![codecov](https://codecov.io/gh/SciML/PDEBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/PDEBase.jl)
[![Build Status](https://github.com/SciML/PDEBase.jl/workflows/Tests/badge.svg)](https://github.com/SciML/PDEBase.jl/actions?query=workflow%3ATests)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

PDEBase.jl provides common types and an interface for building discretizers of ModelingToolkit PDESystems. It serves as a foundation for PDE discretization packages in the SciML ecosystem, such as [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl).

## Overview

This package defines:

- Abstract types for discretization (`AbstractEquationSystemDiscretization`, `AbstractOptimizationSystemDiscretization`)
- Abstract types for discrete spaces (`AbstractDiscreteSpace`, `AbstractCartesianDiscreteSpace`)
- Boundary condition types (`LowerBoundary`, `UpperBoundary`, `InterfaceBoundary`)
- Variable mapping utilities (`VariableMap`, `PeriodicMap`)
- Symbolic utilities for parsing and transforming PDE systems
- Default interface functions for implementing custom discretizers

## Installation

```julia
using Pkg
Pkg.add("PDEBase")
```

## Usage

PDEBase.jl is primarily intended as a dependency for PDE discretization packages. If you are building a discretizer for ModelingToolkit PDESystems, you can use PDEBase.jl as a foundation by implementing its interface functions.

For solving PDEs directly, consider using [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl), which builds on PDEBase.jl.

## Contributing

If you are using this package or want to use it, please open an issue so we are aware. Contributions and feedback are welcome!


