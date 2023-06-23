struct VariableMap
    """
    Calls of the dependent variable, not including boundary terms
    """
    ū
    """
    Independent variables, excluding time
    """
    x̄
    """
    Parameters
    """
    ps
    """
    Time independent variable (nothing, if not applicable)
    """
    time
    """
    Dictionary mapping independent variables to their domain
    """
    intervals
    """
    Dictionary mapping dependent variables to the independent variables that appear in their arguments
    """
    args
    """
    ??? Dependent variables ???
    """
    depvar_ops
    """
    Dictionary mapping independent variables to indices
    """
    x2i
    """
    Dictionary mapping indices to independent variables
    """
    i2x
    """
    ???
    """
    replaced_vars
end

function VariableMap(pdesys, disc; replaced_vars = Dict())
    time = safe_unwrap(get_time(disc))
    eqs = pdesys.eqs
    bcs = pdesys.bcs
    depvars = pdesys.dvs
    domain = pdesys.domain
    ps = pdesys.ps

    if ps isa SciMLBase.NullParameters
        ps = []
    end
    ps = map(ps) do p
        safe_unwrap(p.first)
    end
    depvar_ops = get_ops(depvars)
    # Get all dependent variables in the correct type
    alldepvars = get_all_depvars([eqs; bcs], depvar_ops)
    # Filter out boundaries
    ū = filter(u -> !any(x -> x isa Number, arguments(u)), alldepvars)
    # Get all independent variables in the correct type
    allivs = get_all_indvars(alldepvars)
    x̄ = remove(allivs, time)
    intervals = Dict(map(allivs) do x
        xdomain = domain[findfirst(d -> isequal(x, d.variables), domain)]
        x => (DomainSets.infimum(xdomain.domain), DomainSets.supremum(xdomain.domain))
    end)
    nspace = length(x̄)
    args = [
            u => get_all_indvars(
                    filter(_u -> isequal(operation(_u), u), alldepvars)
                ) 
            for u in depvar_ops
        ] # [operation(u) => arguments(u) for u in ū]
    x̄2dim = [x̄[i] => i for i in 1:nspace]
    dim2x̄ = [i => x̄[i] for i in 1:nspace]
    return VariableMap(ū, x̄, ps, time, Dict(intervals), Dict(args), depvar_ops, Dict(x̄2dim), Dict(dim2x̄), replaced_vars)
end

VariableMap(pdesys) = VariableMap(pdesys, nothing)

function update_varmap!(v::VariableMap, newdv)
    push!(v.ū, newdv)
    merge!(v.args, Dict(operation(newdv) => arguments(newdv)))
    push!(v.depvar_ops, operation(safe_unwrap(newdv)))
end


ivs(u, v::VariableMap) = remove(v.args[operation(u)], v.time)

Base.ndims(u, v::VariableMap) = length(ivs(u, v))

all_ivs(v::VariableMap) = v.time === nothing ? v.x̄ : v.x̄ ∪ [v.time]

all_ivs(u, v::VariableMap) = v.args[operation(u)]

depvar(u, v::VariableMap) = operation(u)(v.args[operation(u)]...)

depvars(v::VariableMap) = v.ū

indvars(v::VariableMap) = v.x̄

x2i(v::VariableMap, u, x) = findfirst(isequal(x), remove(v.args[operation(u)], v.time))

@inline function axiesvals(v::VariableMap, u_, x_, I)
    u = depvar(u_, v)
    map(ivs(u, v)) do x
        x => (I[x2i(v, u, x)] == 1 ? v.intervals[x][1] : v.intervals[x][2])
    end
end

import Base.==

function ==(a::VariableMap, b::VariableMap)
    isequal(Set(a.ū), Set(b.ū)) &&
    isequal(Set(a.x̄), Set(b.x̄)) &&
    isequal(Set(a.ps), Set(b.ps)) &&
    isequal(a.time, b.time) &&
    isequal(a.intervals, b.intervals) &&
    isequal(keys(a.args), keys(b.args)) &&
    all([isequal(Set(a.args[k]), Set(b.args[k])) for k in keys(a.args)]) &&
    isequal(a.depvar_ops, b.depvar_ops) &&
    isequal(a.x2i, b.x2i) &&
    isequal(a.i2x, b.i2x) &&
    isequal(a.replaced_vars, b.replaced_vars)
end