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


# buffer_d, buffer_κs_d, buffer_κs_β, views_d, views_κs_d,
# views_κs_β = NMRSignalSimulator.extractshiftviews(Bs)

shifts, phases, T2s = MSS.shifts, MSS.phases, MSS.T2s

N_κs_d = NMRSignalSimulator.getNvars(shifts)

N_κs_β = NMRSignalSimulator.getNvars(phases)

N_κs_λ = NMRSignalSimulator.getNvars(T2s)
@assert sum(length( As[n].αs ) for n in eachindex(As)) == N_κs_λ

# for coherence shift.

mapping = NMRSignalSimulator.getParamsMapping(shifts, phases, T2s)


# test loading parameters.
#= p = randn(N_κs_d+N_κs_β+N_κs_λ)
#p = collect(1:(N_κs_d+N_κs_β+N_κs_λ))

NMRSignalSimulator.updateparameters!(T2s, p, mapping.T2)
NMRSignalSimulator.updateparameters!(shifts, p, mapping.shift)
NMRSignalSimulator.updateparameters!(phases, p, mapping.phase)

p_test, maping_test = NMRSignalSimulator.testupdateparameters1!(shifts, phases, T2s)
discrepancy = NMRSignalSimulator.testupdateparameters2(p_test, shifts, phases, T2s)
@show discrepancy =#

# # forward model eval.
# get q(u, κ_d, κ_β, κ_λ) and verify against q from Bs, As, from lorentzian.jl exmplae.

q_U = q.(U_rad)

h = uu->NMRSignalSimulator.evalforwardmodel(uu, MSS, MS)

h_U = h.(U_rad)

@show norm(h_U-q_U)

using BenchmarkTools
@btime q.(U_rad)
@btime h.(U_rad)

@assert 1==2

# suppose cos_β and sin_β are updated.

### test derivatives.
n = 1
i = 1
k = 1
sr = MSS.srs[n][i][k]
si = MSS.sis[n][i][k]
qr = (rr, ξξ, bb)->real(evalq(sr, si, rr, ξξ, MSS.phases[n][i][k].cos_β, MSS.phases[n][i][k].sin_β))
qi = (xrr, ξξ, bb)->imag(evalq(sr, si, r, ξ, MSS.phases[n][i][k].cos_β, MSS.phases[n][i][k].sin_β))

# gradient without β.
∇sr = MSS.∇srs![n][i][k]
∇si = MSS.∇sis![n][i][k]

u_rad_test = ppm2hzfunc(3.3984) *2*π
κ_λ_test = 1.12
β_test = 2.31


# use FiniteDifference.jl instead for accuracy control.
# https://github.com/JuliaDiff/FiniteDifferences.jl
buf_r = zeros(2)
buf_i = zeros(2)

phases[n].cos_β[i][k] = cos(β_test)
phases[n].sin_β[i][k] = sin(β_test)
∂qr_∂x, ∂qr_∂κ_λ, ∂qr_∂β = NMRSignalSimulator.derivative!(
    buf_r,
    buf_i,
    u_rad_test,
    κ_λ_test,
    Bs[n].λ0,
    ∇sr,
    ∇si,
    sr,
    si,
    phases,
    n,
    i,
    k,
)

∂qr_eval = [∂qr_∂x; ∂qr_∂κ_λ; ∂qr_∂β]


# verify the derivatives 
import FiniteDiff
import FiniteDifferences


x = [u_rad_test; κ_λ_test; β_test]

#df_x_buffer = similar(x)
#df_x = FiniteDiff.finite_difference_gradient!(df_x_buffer, xx->qr(xx...), x)

dqr_x_AN = ∂qr_eval

ND_accuracy_order = 8
dqr_x_ND = FiniteDifferences.grad(FiniteDifferences.central_fdm(ND_accuracy_order, 1), xx->qr(xx...), x)[1]

dqi_x_ND = FiniteDifferences.grad(FiniteDifferences.central_fdm(ND_accuracy_order, 1), xx->qi(xx...), x)[1]

#@show norm(df_x_buffer - df_x_AN)
@show norm(df_x_ND - df_x_AN)
println("df_x_ND")
display(df_x_ND)

println("df_x_AN")
display(df_x_AN)


println("[df_x_ND df_x_AN]")
display([df_x_ND df_x_AN])

println("df_x_ND ./ df_x_AN")
display(df_x_ND ./ df_x_AN)

@assert 1==2

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