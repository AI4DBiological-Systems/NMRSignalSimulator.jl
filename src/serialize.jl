
saveasJSON = NMRHamiltonian.saveasJSON

abstract type SerializationTrait end
struct BasicSerialization <: SerializationTrait end # non-Vector types, easy to map to JSON format.
struct CustomSerialization <: SerializationTrait end 


# Default behavior is to use custom routine.
SerializationTrait(::Type{T}) where T = CustomSerialization()

# abstract types that can be serialized via the basic template.
SerializationTrait(::Type{T}) where T <: MoleculeParams = BasicSerialization()

# specific types that can be serialized via the basic template.
SerializationTrait(::Type{CLOperationRange{T}}) where T = BasicSerialization()
#SerializationTrait(::Type{InterpolationSamplesSerialSafe{T}}) where T = BasicSerialization()
#SerializationTrait(::Type{InterpolationSamples{T}}) where T = CustomSerialization() # JSON needs explicit separation of real and imaginary numbers.

SerializationTrait(::Type{InterpolationSamples{T}}) where T = BasicSerialization()


### default, catch-all.

function serializetype(x::DT) where DT
    return serializetype(SerializationTrait(DT), x)
end

function deserializetype(X::DT) where DT

    return deserializetype(recovertype(X[:own_type]), X)
end

function deserializetype(type_tag::String, args...)
    deserializetype(recovertype(type_tag), args...)
end

# vector versions.
#function serializetype(Xs::VT) where VT <: Union{Vector,JSON3.Array}
function serializetype(Xs::Vector)

    return collect( serializetype(Xs[n]) for n in eachindex(Xs) )
    #return Dict( vector => serializetype(Xs[n]) for n in eachindex(Xs) )
end

#function deserializetype(Ws::VT) where VT <: Union{Vector,JSON3.Array}
function deserializetype(Ws::Vector)

    return collect( deserializetype(Ws[n]) for n in eachindex(Ws) )
    #return collect( deserializetype(Ws[n]) for n in eachindex(Ws) )
end

# For JSON3 or a serialized dictionary-like data structure that stores `DataType` variables as `String`.

function recovertype(tag::String)::DataType
    return eval(Meta.parse(tag))
end

function recovertype(tag::DataType)::DataType
    return tag
end


### basic serialization: fields are all of same type or (nested)-arrays of the same type.

function serializetype(::BasicSerialization, A)

    vals = Dict(
        x => getfield(A, x)
        for x in propertynames(A)
    )

    types = Dict(
        Symbol("$(String(x))_type") => typeof(getfield(A, x))
        for x in propertynames(A)
    )
    own_type = Dict(:own_type => typeof(A))
    return merge(vals, types, own_type)
end

function deserializetype(::Type{DT}, W) where DT

    return DT(
        collect(
            convert( recovertype(W[Symbol("$(String(x))_type")]), W[x] ) # recovertype() needed in case `DataType` variables are stored as `String` or other types, and need to be converted to `DataType`. This is the case for JSOn3.jl.
            for x in fieldnames(DT)
        )...
    )
end

function deserializetype(::Type{Vector{T}}, W) where T <: Number

    return collect(
        convert( T, W[i] )
        for i in eachindex(W)
    )
end

# function deserializetype(::Type{Vector{T}}, W) where T

#     return collect(
#         deserializetype( T, W[i] )
#         for i in eachindex(W)
#     )
# end

### custom serialization routines.

function serializetype(::CustomSerialization, A::SpinSysParams)

    return Dict(
        :shift => serializetype(A.shift),
        :shift_type => typeof(A.shift),
        
        :phase => serializetype(A.phase),
        :phase_type => typeof(A.phase),
        
        :T2 => serializetype(A.T2),
        :T2_type => typeof(A.T2),

        :own_type => typeof(A),
    )
end

function deserializetype(
    ::Type{SpinSysParams{ST,PT,T2T}},
    W,
    )::SpinSysParams{ST,PT,T2T} where {ST,PT,T2T}

    shift = deserializetype(W[:shift_type], W[:shift])
    phase = deserializetype(W[:phase_type], W[:phase])
    T2 = deserializetype(W[:T2_type], W[:T2])

    return SpinSysParams(shift, phase, T2)
end

###### models.

### mixture: As::SHType{T}.

"""
```
function serializclproxies(
    Bs::Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},OT}},
) where {T,ST,PT,T2T,OT}
```
    
    Returns a dictionary containing the surrogate model.
"""
function serializclproxies(
    Bs::Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T},OT}},
    ) where {T,ST,PT,T2T,OT}

    if isempty(Bs)
        println("input is an empty collection. Cannot serialize.")
    else

        # prefix S_ stands for serialized (in a Dict{:Symbol, .} format) version of a variable.

        S_ss_params_set = collect( serializetype(Bs[n].ss_params) for n in eachindex(Bs) )
        S_op_range_set = collect( serializetype(Bs[n].op_range) for n in eachindex(Bs) )
  
        # S_itp_samps_set = serializetype(itp_samps_set)

        λ0 = first(Bs).λ0

        return Dict(
            # values.
            :S_ss_params_set => S_ss_params_set,
            :S_op_range_set => S_op_range_set,
            # :S_itp_samps_set => S_itp_samps_set,
            :λ0 => λ0,

            # types.
            :ST => ST,
            :PT => PT,
            :T2T => T2T,
            :T => T,

            :S_ss_params_set_type => typeof(S_ss_params_set),
            :S_op_range_set_type => typeof(S_op_range_set),
            # :S_itp_samps_set_type => typeof(S_itp_samps_set),
            :λ0_type => typeof(λ0),
        )
    end

    return nothing
end

"""
```
function deserializclproxies(W)
```

W is of a data type that can be addressed via a key (of type Symbol) and returns a blue.
Examples: type Dict{Symbol, Any} or JSON.Object.
"""
function deserializclproxies(W)

    N = length(W[:S_ss_params_set])
    @assert length(W[:S_op_range_set]) == N

    λ0 = convert(recovertype(W[:λ0_type]), W[:λ0])

    T = recovertype(W[:T])
    ST = recovertype(W[:ST])
    PT = recovertype(W[:PT])
    T2T = recovertype(W[:T2T])

    @assert T <: AbstractFloat
    @assert ST <: MoleculeParams
    @assert PT <: MoleculeParams
    @assert T2T <: MoleculeParams

    ss_params_set = collect( deserializetype(W[:S_ss_params_set][n]) for n in eachindex(W[:S_ss_params_set]))
    op_range_set = collect( deserializetype(W[:S_op_range_set][n]) for n in eachindex(W[:S_op_range_set]))

    return ss_params_set, op_range_set, λ0
end

"""
```
function serializitpsamples(
    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},
) where T
```
    
Returns a dictionary of the interpolation samples used to fit the surrogate.
"""
function serializitpsamples(
    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},
    ) where T

    if isempty(itp_samps_set)
        println("input is an empty collection. Cannot serialize.")
    else

        # prefix S_ stands for serialized (in a Dict{:Symbol, .} format) version of a variable.
        S_itp_samps_set = serializetype(itp_samps_set)

        return Dict(
            # values.
            :S_itp_samps_set => S_itp_samps_set,
            :S_itp_samps_set_type => typeof(S_itp_samps_set),
            :own_type => typeof(itp_samps_set),
        )
    end

    return nothing
end

"""
```
function deserializitpsamples(W)
```

W is of a data type that can be addressed via a key (of type `Symbol`).
Examples: type `Dict{Symbol, Any}` or `JSON.Object`.
"""
function deserializitpsamples(W)

    itp_samps_set = collect( deserializetype(W[:S_itp_samps_set][n]) for n in eachindex(W[:S_itp_samps_set]))

    return itp_samps_set
end


### comparison between serialized objects. For round trip testing.

function compareserialized(x::Vector, y::Vector, args...)

    @assert length(x) == length(y)

    return all(
        compareserialized(x[n], y[n], args...)
        for n in eachindex(x)
    )
end

function compareserialized(x::DataType, y::DataType, args...)::Bool
    return x == y
 end

function compareserialized(x::Array{T,D}, y::Array{T,D}, zero_tol)::Bool where {T <: Number, D}
   return norm(x-y) < zero_tol
end

function compareserialized(x::T, y::T, zero_tol)::Bool where T <: Number
    return abs(x-y) < zero_tol
end

# """
# ```
# compareserialized(
#     X::Dict{Symbol, Any},
#     Y::Dict{Symbol, Any},
#     args...
#     )::Bool
# ```
# For testing `serializetype()` and `deserializetype()`.

# For now, `args...` is just a floating point `zero_tol` input that should be a very small number. More documentation here as we develop the serialization routines.

# Example:
# ```
# using Test

# op1 = NMRSignalSimulator.CLOperationRange(rand(), rand(), randn(4), rand(), rand(), rand(), rand())
# S_op_range = NMRSignalSimulator.serializetype(op1)

# op2 = NMRSignalSimulator.deserializetype(typeof(op1), S_op_range)

# W1 = NMRSignalSimulator.serializetype(op1)
# W2 = NMRSignalSimulator.serializetype(op2)

# zero_tol = 1e-12
# @test NMRSignalSimulator.compareserialized(W1, W2, zero_tol)
# ```
# """
function compareserialized(
    X::Dict{Symbol, Any},
    Y::Dict{Symbol, Any},
    args...
    )::Bool

    @assert keys(X) == keys(Y)

    return all(
        compareserialized(val, Y[key], args...)
        for (key,val) in X
    )
end

# assumes the function output is numerical.
function compareserialized(x::Function, y::Function, zero_tol, U_rad)

    return norm(x.(U_rad) - y.(U_rad)) < zero_tol
end