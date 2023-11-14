"""
```
function recoverclproxies(
    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},
    ss_params_set::Vector{SpinSysParams{ST,PT,T2T}},
    op_range_set::Vector{CLOperationRange{T}},
    As::Vector{HAM.SHType{T}},
    λ0::T,
    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},CLOperationRange{T}}},
    CLMixtureSpinSys{T,ST,PT,T2T}} where {T,ST,PT,T2T}
```

Takes in deserialized quantities to output a surrogate model and its corresponding mixture spin system data structure.
"""
function recoverclproxies(
    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},
    ss_params_set::Vector{SpinSysParams{ST,PT,T2T}},
    op_range_set::Vector{CLOperationRange{T}},
    As::Vector{HAM.SHType{T}},
    λ0::T,
    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},CLOperationRange{T}}},
    CLMixtureSpinSys{T,ST,PT,T2T}} where {T,ST,PT,T2T}

    N = length(As)

    Bs = Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},CLOperationRange{T}}}(undef, N)

    srs = Vector{Vector{Vector{Function}}}(undef, N)
    sis = Vector{Vector{Vector{Function}}}(undef, N)

    ∇srs! = Vector{Vector{Vector{Function}}}(undef, N)
    ∇sis! = Vector{Vector{Vector{Function}}}(undef, N)

    for n in eachindex(As)

        Bs[n], srs[n], sis[n], ∇srs![n], ∇sis![n] = recoverclproxy(
            itp_samps_set[n],
            ss_params_set[n],
            op_range_set[n],
            As[n],
            λ0,
        )
    end

    # alternative specification for derivatives.
    shifts, phases, T2s, Δc_bars = setupSSvars(As, Bs)
    mixSS = CLMixtureSpinSys(
        srs,
        sis,
        ∇srs!,
        ∇sis!,
        shifts,
        phases,
        T2s,
        Δc_bars,
        Bs[begin].λ0,
    )

    return Bs, mixSS
end

function recoverclproxy(
    itp_samps::Vector{InterpolationSamples{T}},
    ss_params::SpinSysParams,
    op_range::CLOperationRange{T},
    A::HAM.SHType{T},
    λ0::T,
    ) where T <: AbstractFloat

    qs, sr, si, ∇sr!, ∇si! = setupclmoleculepartitionitp(
        itp_samps,
        ss_params.phase.cos_β,
        ss_params.phase.sin_β,
        CLSurrogateSpinSysInputs(A.parts, A.αs, A.Ωs, op_range.u_min, op_range.u_max),
    )

    core = MoleculeType(
        qs,
        ss_params,
        op_range,
        λ0,
    )

    return core, sr, si, ∇sr!, ∇si!
end