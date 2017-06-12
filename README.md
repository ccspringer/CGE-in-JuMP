# CGE-with-JuMP

This is a simple CGE model implemented in Julia for Mathematical Optimization (JuMP).

## Requirements

The minimum requirements to run the model are [Julia](http://julialang.org/) and the [JuMP](https://github.com/JuliaOpt/JuMP.jl) package. To read in data, you should also install the [DataFrames](https://github.com/JuliaStats/DataFrames.jl) package.

## Summary ##

This model uses the JuMP package with the Ipopt solver to demonstrate an open-source alternative to GAMS.

The model has firms and a representative household. The firms each produce a single commodity for their sector. The household is the sole consumer of these three commodities. Firms maximize profit according to a constant elasticity of substitution (CES) production function and the household maximizes utility according to a Linear Expenditure System (LES) utility function. The LES utility function is subject to a household budget constraint.  

For ease of use and expansion, the model economy is roughly divided into two main activities – production and consumption - which are annotated as comments within the parameter-variable-equation model structure commonly used in GAMS.

Parameters were constructed to match the sectoral index used to organize input data. The model specifies a sector index for variables that matches the index in the parameter data. The model is calibrated to example data with three sectors, which users can change or expand.

As specified, the model replicates baseline data. Users can specify shocks to the economy and check results using the ``getvalue(variable_name)`` function.

Feedback on this project is welcome.
