using ModelingToolkit, PDEBase

## Approximation of function 1D
println("Approximation of function 1D")

@parameters x
@variables u(..)

func(x) = @. 2 + abs(x - 0.5)

eq = [u(x) ~ func(x)]
bc = [u(0) ~ u(0)]

x0 = 0
x_end = 2
domain = [x ∈ (x0, x_end)]

@named pde_system = PDESystem(eq, bc, domain, [x], [u(x)])

ū = [u(x)]
x̄ = [x]
ps = []
time = nothing
intervals = Dict(x => (x0, x_end))
depvar_ops = [operation(safe_unwrap(u(x)))]
args = Dict(depvar_ops[1] => [x])
x2i = Dict(x => 1)
i2x = Dict(1 => x)
replaced_vars = Dict()

v_generated = VariableMap(pde_system) 
v_true = VariableMap(ū, x̄, ps, time, intervals, args, depvar_ops, x2i, i2x, replaced_vars)
@test v_generated == v_true

## Approximation of function 2D
println("Approximation of function 2D")

@parameters x, y
@variables u(..)
func(x, y) = -cos(x) * cos(y) * exp(-((x - pi)^2 + (y - pi)^2))
eq = [u(x, y) ~ func(x, y)]
bc = [u(0, 0) ~ u(0, 0)]

x0 = -10
x_end = 10
y0 = -10
y_end = 10
d = 0.4

domain = [x ∈ (x0, x_end), y ∈ (y0, y_end)]
@named pde_system = PDESystem(eq, bc, domain, [x, y], [u(x, y)])

ū = [u(x, y)]
x̄ = [x, y]
ps = []
time = nothing
intervals = Dict([x => (x0, x_end), y => (y0, y_end)])
depvar_ops = [operation(safe_unwrap(u(x,y)))]
args = Dict(depvar_ops[1] => [x, y])
x2i = Dict([x => 2, y => 1])
i2x = Dict([2 => x, 1 => y])
replaced_vars = Dict()

v_generated = VariableMap(pde_system) 
v_true = VariableMap(ū, x̄, ps, time, intervals, args, depvar_ops, x2i, i2x, replaced_vars)
@test v_generated == v_true

## Approximation of function 1D 2
println("Approximation of function 1D 2")

@parameters x
@variables u(..)

func(x) = @. abs(x - 0.5) - 0.5

eq = [u(2x - u(x-3)) - u(0) ~ func(2x)]
bc = [u(0) ~ 2.5]

x0 = 0
x_end = 2
domain = [x ∈ (x0, x_end)]

@named pde_system = PDESystem(eq, bc, domain, [x], [u(x)])

ū = [u(2x - u(x-3)), u(x-3)]
x̄ = [x]
ps = []
time = nothing
intervals = Dict(x => (x0, x_end))
depvar_ops = [operation(safe_unwrap(u(x)))]
args = Dict(depvar_ops[1] => [x])
x2i = Dict(x => 1)
i2x = Dict(1 => x)
replaced_vars = Dict()

v_generated = VariableMap(pde_system) 
v_true = VariableMap(ū, x̄, ps, time, intervals, args, depvar_ops, x2i, i2x, replaced_vars)
@test v_generated == v_true
