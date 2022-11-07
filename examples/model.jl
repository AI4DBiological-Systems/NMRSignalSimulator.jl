# run lorentzian.jl first.

include("./helpers/utils.jl")

q2 = uu->forwardmodelmolecule(uu, As, Bs)

q2_U = q2.(U_rad)
@show norm(q2_U-q_U)

#vals = As[1].αs
#vals = Bs[1].ss_params.d

κds, κβs = NMRSignalSimulator.extractshifts(Bs)

vals = κds

println("flatten state from nested Vectors to a Vector.")
@show flat_vals = collect( Base.Iterators.flatten(vals) )
@show vals
println()

#flat_vals .+= 1

# I am here. write a type stable, mutating recursive iterator that loops through a nested vector given a flat vector.
# - a versoin that mutates flat vector from st_ind, reading in from nested.
# - a version that mutates nested vector, reading in from flat vector from st_ind.
## not type safe.
# generic read with nested vectors.
# assumes:
#   - each parameter is a 1D array.
#   - every name entry has the same N_parameters_per_entry number categories of parameters.
# x[j,n] is the parameter for the j-th parameter category for the n-th name entry.
function convertnestedvectors(c::JSON3.Array, val_type::DataType)

    return collect( nestedconvert(c, i, val_type) for i in eachindex(c) )
end

# this is not type stable..
# This is for nested Vector{Any} but really have base Julia concrete types at the inner-most level.
# position takes value in {1, 2, ..., number of elements of x}.
#   - e.g., position is an element of the eachindex(x) iterator.
function nestedconvert(x::Vector{Any}, position, val_type::DataType)

    if isempty(x)
        ret = Vector{val_type}(undef, 0)
        return convert(Vector{val_type}, ret)
    end

    y = x[begin+position-1]
    current_type = typeof(y)

    while current_type <: JSON3.Array
        
        return convertnestedvectors(y, val_type) # returns a Vector{} or Vctor{Vector{T}} or higher-order nested Vector{T}'s.
    end

    return convert(val_type, y)
end

function nestedconvert(x::JSON3.Array, position, val_type::DataType)

    if isempty(x)
        ret = Vector{val_type}(undef, 0)
        return convert(Vector{val_type}, ret)
    end

    y = x[begin+position-1]
    current_type = typeof(y)

    while current_type <: JSON3.Array

        return convertnestedvectors(y, val_type) # returns a Vector{} or Vctor{Vector{T}} or higher-order nested Vector{T}'s.
    end

    return convert(val_type, y)
end