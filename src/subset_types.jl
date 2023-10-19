### subset variable update data structures.
# see updatewandderivatives!()


"""
```
struct SubsetVarsIndices
    shift::Vector{Int}
    phase::Vector{Int}
    T2::Vector{Int}
```
"""
struct SubsetVarsIndices
    shift::Vector{Int}
    phase::Vector{Int}
    T2::Vector{Int}
end

"""
```
function exportindices(A::SubsetVarsIndices)::Vector{Int}
```
    Returns [A.shift; A.phase; A.T2]
"""
function exportindices(A::SubsetVarsIndices)::Vector{Int}
    return [A.shift; A.phase; A.T2]
end

####

"""
```
abstract type VariableSetTrait end
```
"""
abstract type VariableSetTrait end

"""
```
struct AllVars <: VariableSetTrait end
```
"""
struct AllVars <: VariableSetTrait end

function getNvars(::AllVars, MSS::MixtureSpinSys)::Int
    return getNvars(MSS)
end

"""
```
struct SubsetVars <: VariableSetTrait
    indices::SubsetVarsIndices
    active_systems::Vector{Tuple{Int,Int}}
```
"""
struct SubsetVars <: VariableSetTrait
    indices::SubsetVarsIndices
    active_systems::Vector{Tuple{Int,Int}}
end

"""
``` function exportindices(A::SubsetVars)::Vector{Int} ```
    Returns index in A.
"""
function exportindices(A::SubsetVars)::Vector{Int}
    return exportindices(A.indices)
end

"""
``` function getindices(A::AllVars)::AllVars ```
    returns A.
"""
function getindices(A::AllVars)::AllVars
    return A
end

"""
``` function getindices(A::SubsetVars)::SubsetVarsIndices ```
    returns A.indices.
"""
function getindices(A::SubsetVars)::SubsetVarsIndices
    return A.indices
end

"""
``` function getactivesystems(A::AllVars)::AllVars ```
    returns A.
"""
function getactivesystems(A::AllVars)::AllVars
    return A
end

"""
``` function getactivesystems(A::SubsetVars)::Vector{Tuple{Int,Int}} ```
    returns A.active_systems.
"""
function getactivesystems(A::SubsetVars)::Vector{Tuple{Int,Int}}
    return A.active_systems
end

function getNvars(A::SubsetVars, args...)::Int
    return length(A.active_systems)
end

"""
```
function setupsubsetvars(
    indicies_input::Vector{Int},
    mapping::ParamsMapping;
)::SubsetVars
```

Convinence constructor for `SubsetVars`.
"""
function setupsubsetvars(
    indicies_input::Vector{Int},
    mapping::ParamsMapping;
    )::SubsetVars
    
    shift_range = getshiftrange(mapping)
    phase_range = getphaserange(mapping)
    T2_range = getT2range(mapping)

    var_inds = sort(unique(indicies_input))

    shift_inds = Vector{Int}(undef, 0)
    phase_inds = Vector{Int}(undef, 0)
    T2_inds = Vector{Int}(undef, 0)
    active_systems = Vector{Tuple{Int,Int}}(undef, 0)

    n = 0
    i = 0
    for k in eachindex(var_inds)
        ind = var_inds[k]
        if ind in shift_range
            n, i = findspinsystem(ind, mapping.shift)
            push!(shift_inds, ind)
            
        elseif ind in phase_range
            n, i = findspinsystem(ind, mapping.phase)
            push!(phase_inds, ind)

        elseif ind in T2_range
            n, i = findspinsystem(ind, mapping.T2)
            push!(T2_inds, ind)

        else
            "return error"
            n, i = 0, 0
        end
        push!(active_systems, (n,i))

    end

    indices = SubsetVarsIndices(shift_inds, phase_inds, T2_inds)

    return SubsetVars(indices, active_systems)
end

function findspinsystem(ind::Integer, mapping::MoleculeParamsMapping)::Tuple{Int,Int}
    for n in eachindex(mapping.st)
        for i in eachindex(mapping.st[n])

            if ind in (mapping.st[n][i]):(mapping.fin[n][i])
                return (n,i)
            end
        end
    end

    "return error"
    return n, i
end