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

### custom parameters start.
# test loading parameters.
p = randn(N_κs_d+N_κs_β+N_κs_λ)
#p = collect(1:(N_κs_d+N_κs_β+N_κs_λ))

NMRSignalSimulator.updateparameters!(T2s, p, mapping.T2)
NMRSignalSimulator.updateparameters!(shifts, p, mapping.shift)
NMRSignalSimulator.updateparameters!(phases, p, mapping.phase)

p_test, maping_test = NMRSignalSimulator.testupdateparameters1!(shifts, phases, T2s)
discrepancy = NMRSignalSimulator.testupdateparameters2(p_test, shifts, phases, T2s)
@show discrepancy
### custom parameters end.

# # forward model eval.
# get q(u, κ_d, κ_β, κ_λ) and verify against q from Bs, As, from lorentzian.jl exmplae.

q_U0 = copy(q_U)
q_U = q.(U_rad)

h = uu->NMRSignalSimulator.evalforwardmodel(uu, MSS, MS)

h_U = h.(U_rad)

@show norm(h_U-q_U)


using BenchmarkTools
println("Timing: q")
@btime q.(U_rad)
println("Timing: h")
@btime h.(U_rad)
println("Timing: g")
@btime g.(U_rad)
#h and g are about the same speed. Both faster than q.


# try flat version.
FM = NMRSignalSimulator.FlatMixtureModel(As, MSS, MS)

g = uu->NMRSignalSimulator.evalforwardmodel(uu, FM, MS)
g_U = g.(U_rad)

@show norm(g_U-q_U)


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_display, real.(g_U), label = "g")
PyPlot.plot(P_display, real.(q_U), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("g vs q")


@assert 1==2

# using BenchmarkTools
# println("Timing: q")
# @btime q.(U_rad)
# println("Timing: h")
# @btime h.(U_rad)
# println("Timing: g")
# @btime g.(U_rad)
# #h and g are about the same speed. Both faster than q.

### why is it so slow?
l = 1
u_rad_test = ppm2hzfunc(3.3984) *2*π

out = NMRSignalSimulator.evalq(
    FM.srs[l],
    FM.sis[l],
    u_rad_test - FM.ds[l],
    FM.ξs[l],
    FM.cos_βs[l],
    FM.sin_βs[l],
)

@btime NMRSignalSimulator.evalq(
    FM.srs[l],
    FM.sis[l],
    u_rad_test - FM.ds[l],
    FM.ξs[l],
    FM.cos_βs[l],
    FM.sin_βs[l],
)

@btime FM.srs[l](u_rad_test, FM.ξs[l])


real_setp, imag_setp, _ = loadalphaglucose()

out2r = real_setp(u_rad_test - FM.ds[l], FM.ξs[l])
out2i = imag_setp(u_rad_test - FM.ds[l], FM.ξs[l])

@btime real_setp(u_rad_test - FM.ds[l], FM.ξs[l])


ds = FM.ds
@btime u_rad_test - ds[3]


x = randn()
y = randn()
z = u_rad_test - ds[3]
@btime real_setp(z, x)
@btime x - y 



@btime x-y setup=(x=rand(),y=randn()) # 17ns
@btime real_setp(x,y) setup=(x=ppm2hzfunc(rand()+3)*2*π,y=randn()) # 61 ns
@btime real_setp(x,y) setup=(x=rand(),y=randn()) # 62ns.
@btime real_setp(u_rad_test - FM.ds[l], FM.ξs[l]) setup=(l=rand(1:length(FM.ds))) # 227ns.

@btime real_setp(FM.ds[l], FM.ξs[l]) setup=(l=rand(1:length(FM.ds))) # 194.

# idea: f.(flat) to speed up.

