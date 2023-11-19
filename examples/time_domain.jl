
Random.seed!(25)

PLT.close("all")
fig_num = 1



N_t = 2^15
offset_ppm = convert(T, 0.3)

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

t_range = DSU.gettimerange(N_t, fs)

offset_Hz = ν_0ppm - (ppm2hzfunc(offset_ppm)-ppm2hzfunc(zero(T)))
U_DFT, U_y, DFT2y_inds = DSU.getwraparoundDFTfreqs(N_t, fs, offset_Hz)
y2DFT_inds = DSU.getinversemap(DFT2y_inds)


FID_proxy_config = SIG.FIDSurrogateConfig{T}(
    Δt = convert(T, 1e-5),
    t_lb = zero(T),
    t_ub = convert(T, 3.0),
)

Cs, MSS_FID, itp_samps_FID = SIG.fitfidproxies(As, λ0, FID_proxy_config)

model_params_FID = SIG.MixtureModelParameters(MSS_FID, copy(w_oracle))
SIG.importmodel!(model_params_FID, x_oracle)


### plot.

# reference.
U_y_rad = U_y .* convert(T, 2*π)
q_U_y = q.(U_y_rad)

q_U_DFT = q_U_y[y2DFT_inds]

# The discrete-time Fourier transform (DTFT) is an approximation of the Fourier Transform that is scaled by 1/fs.
# Need to boost it by fs to match the scale of the FID model.
ifft_q_t = ifft(q_U_DFT) .* fs

# FID surrogates: check spline approximation.
w_fid = copy(w_oracle)
f_fid = uu->SIG.evalfidmixture(uu, As, Cs; w = w_fid)
g = uu->SIG.evalfidproxymixture(uu, Cs; w = w_fid)

g_t = g.(t_range)

# verify g is an approximation of f.
f_fid_t = f_fid.(t_range)
discrepancies = abs.(f_fid_t-g_t)
max_val, ind = findmax(discrepancies)
println("Spline approximation for surrogate model:")
println("relative l-2 discrepancy = ", norm(f_fid_t-g_t)/norm(f_fid_t))
println("max l-2 discrepancy: ", max_val, ", at time $(t_range[ind]).")
println()

##
PLT.figure(fig_num)
fig_num += 1

PLT.plot(t_range, real.(ifft_q_t), linewidth = "2.5", label = "ifft(q)")
PLT.plot(t_range, real.(g_t), "--", label = "g")

PLT.legend()
PLT.xlabel("t")
PLT.ylabel("real")
PLT.title("Free-induction decay. q is complex Lorentzian surrogate, g is FID surrogate")

##
PLT.figure(fig_num)
fig_num += 1

PLT.plot(t_range, imag.(ifft_q_t), linewidth = "2.5", label = "ifft(q)")
PLT.plot(t_range, imag.(g_t), "--", label = "g")

PLT.legend()
PLT.xlabel("t")
PLT.ylabel("imag")
PLT.title("Free-induction decay. q is complex Lorentzian surrogate, g is FID surrogate")

##
PLT.figure(fig_num)
fig_num += 1

PLT.plot(t_range, abs.(ifft_q_t), linewidth = "2.5", label = "ifft(q)")
PLT.plot(t_range, abs.(g_t), "--", label = "g")

PLT.legend()
PLT.xlabel("t")
PLT.ylabel("abs")
PLT.title("Free-induction decay. q is complex Lorentzian surrogate, g is FID surrogate")


nothing