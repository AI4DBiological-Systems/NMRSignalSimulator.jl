
# based on serializclproxies()
function exportmixtureproxy(
    Bs::Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},OT}},
    ) where {T,ST,PT,T2T,OT}

    ss_params_set = collect( Bs[n].ss_params for n in eachindex(Bs) )
    op_range_set = collect( Bs[n].op_range for n in eachindex(Bs) )

    λ0 = first(Bs).λ0

    return ss_params_set, op_range_set, λ0
end

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

###### FID case

function recoverfidproxies(
    itp_samps_set::Vector{Vector{FIDInterpolationSamples{T}}},
    ss_params_set::Vector{SpinSysParams{ST,PT,T2T}},
    op_range_set::Vector{FIDOperationRange{T}},
    As::Vector{HAM.SHType{T}},
    λ0::T,
    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},FIDOperationRange{T}}},
    FIDMixtureSpinSys{T,ST,PT,T2T}} where {T,ST,PT,T2T}

    N = length(As)

    Cs = Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},FIDOperationRange{T}}}(undef, N)

    for n in eachindex(As)

        Cs[n] = recoverfidproxy(
            itp_samps_set[n],
            ss_params_set[n],
            op_range_set[n],
            As[n],
            λ0,
        )
    end

    # MSS_fid.
    shifts, phases, T2s, Δc_bars = setupSSvars(As, Cs)
    mixSS = FIDMixtureSpinSys(
        shifts,
        phases,
        T2s,
        Δc_bars,
        λ0,
    )
    return Cs, mixSS
end

function recoverfidproxy(
    itp_samps::Vector{FIDInterpolationSamples{T}},
    ss_params::SpinSysParams,
    op_range::FIDOperationRange{T},
    A::HAM.SHType{T},
    λ0::T,
    ) where T <: AbstractFloat

    qs = setupfidmoleculepartitionitp(
        itp_samps,
        ss_params.phase.cos_β,
        ss_params.phase.sin_β,
        CLSurrogateSpinSysInputs(A.parts, A.αs, A.Ωs, zero(T), zero(T)),
    )

    core = MoleculeType(
        qs,
        ss_params,
        op_range,
        λ0,
    )

    return core, sr, si, ∇sr!, ∇si!
end