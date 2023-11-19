

### import model from p.var_flat. variables: shift (ζ), phase (β), T2 (λ).
"""
```
function importmodel!(p::MixtureModelParameters)::Nothing
```

Updates `p.MSS` with the contents of `p.var_flat`.
"""
function importmodel!(p::MixtureModelParameters)

    updatespinsystems!(p.MSS, p.var_flat, p.systems_mapping)
    resolvespinsystems!(p.MSS)

    #updatesinglets!(p.MS, p.var_flat, p.singlets_mapping)

    return nothing
end

"""
```
function importmodel!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing
```

Only update the entries as specified by inds.
"""
function importmodel!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)

    k = 0
    for i in inds.shift
        k += 1
        p.var_flat[i] = x[k]
    end

    for i in inds.phase
        k += 1
        p.var_flat[i] = x[k]
    end

    for i in inds.T2
        k += 1
        p.var_flat[i] = x[k]
    end

    return importmodel!(p)
end

# set the other parameters to the default value. inefficient, but avoids uninitialized errors.
"""
```
function importmodelreset!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing
```

This is a version of `importmodel!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing` that resets the other variables not mentioned in `inds`.

First, resets the model to default parameter values; 0 for shift and phase parameters, 1 for T2 multiplier parameters.
Then update the model parameters as specified by inds with the entries in x.
"""
function importmodelreset!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)

    resetvarflat!(p)
    return importmodel!(p, x, inds)
end

function importmodelreset!(p::MixtureModelParameters, x, ::AllVars)

    return importmodel!(p, x)
end

"""
```
function importmodel!(p::MixtureModelParameters, x, ::AllVars)::Nothing
```

Calls importmodel!(p, x).
"""
function importmodel!(p::MixtureModelParameters, x, ::AllVars)
    return importmodel!(p, x)
end

"""
```
function importmodel!(p::MixtureModelParameters, x)::Nothing
```

Update the model `p` using the contents of `x`.
"""
function importmodel!(p::MixtureModelParameters, x)

    p.var_flat[:] = x

    return importmodel!(p)
end

function resetvarflat!(p::MixtureModelParameters{T, ST}) where {T, ST}

    fill!(p.var_flat, zero(T))
    
    mapping = getvarrange(p.systems_mapping.T2)
    for i in mapping
        p.var_flat[i] = one(T)
    end

    return nothing
end

"""
```
function exportmodel!(p::MixtureModelParameters)::Nothing
```

Flattens the model parameters into `p.var_flat`.
"""
function exportmodel!(p::MixtureModelParameters)

    exportspinsystems!(p.var_flat, p.MSS, p.systems_mapping)
    #exportsinglets!(p.var_flat, p.MS, p.singlets_mapping)

    return nothing
end

#### single evaluation.

function evalmodel!(
    p::MixtureModelParameters,
    u_rad::T,
    x,
    )::Complex{T} where T <: AbstractFloat
    
    # # check.
    #@assert length(x) == getNvars(p)
    
    # parse.
    p.var_flat[:] = x
    importmodel!(p)
    
    # eval.
    out = evalmodel(u_rad, p)

    return out
end


#### version for multiple evaluations
function evalmodel!(
    q_U::Vector{Complex{T}}, # output.
    p::MixtureModelParameters,
    U_rad,
    x,
    ) where T <: AbstractFloat
    
    # parse.
    p.var_flat[:] = x
    importmodel!(p)
    
    # eval.
    resize!(q_U, length(U_rad))
    for m in eachindex(q_U)
        q_U[m] = evalmodel(U_rad[m], p)
    end

    return nothing
end

# no error checking on the lengths of fieldnames in MSS.
function evalmodel(
    u_rad::T,
    p::MixtureModelParameters,
    )::Complex{T} where T
    
    MSS, w = p.MSS, p.w

    # intermetdiate buffer.
    out_sys = zero(Complex{T})
    #out_singlets = zero(Complex{T})

    # eval.
    out = zero(Complex{T})

    for n in eachindex(w)
        out_sys = evalsystems(
            u_rad,
            MSS.srs[n],
            MSS.sis[n],
            MSS.shifts[n],
            MSS.phases[n],
            MSS.T2s[n],
        )

        out += w[n]*out_sys
    end

    return out
end

function evalsystems(
    u_rad::T,
    sr,
    si,
    shift::SharedShift{T}, # contains d.
    phase::CoherencePhase{T}, # contains β.
    T2::SharedT2{T}, # contains ξ.
    ) where T <: AbstractFloat

    out = zero(Complex{T})

    # # tmp storage.
    r = zero(T)
    λ = one(T)

    for i in eachindex(phase.cos_β)
        
        #ξ = T2.var[i]
        λ = T2.λ[i]
        r = u_rad - shift.var[i]

        for k in eachindex(phase.cos_β[i])

            out += evalq(
                sr[i][k],
                si[i][k],
                r,
                λ,
                phase.cos_β[i][k],
                phase.sin_β[i][k],
            )
        end
    end

    return out
end

function evalsystems(
    u_rad::T,
    sr,
    si,
    shift::CoherenceShift{T}, # contains d.
    phase::CoherencePhase{T}, # contains β.
    T2::SharedT2{T}, # contains ξ.
    ) where T <: AbstractFloat

    out = zero(Complex{T})

    r = zero(T)

    for i in eachindex(phase.cos_β)

        # the number of resonance groups should be the same.
        #@assert length(phase.cos_β[i]) == length(phase.sin_β[i]) == length(shift.ζ[i])

        #ξ = T2.var[i]
        λ = T2.λ[i]

        for k in eachindex(phase.cos_β[i])
            r = u_rad - shift.ζ[i][k]

            # 300 ms.
            out += evalq(
                sr[i][k],
                si[i][k],
                r,
                λ,
                phase.cos_β[i][k],
                phase.sin_β[i][k],
            )

        end
    end

    return out
end


