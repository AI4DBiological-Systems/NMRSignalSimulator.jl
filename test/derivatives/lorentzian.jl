
# ### verify lorentzian derivatves
r_test = ppm2hzfunc(0.001) * 2*π
x_test = [r_test; λ_test]

p_test = [ # 0 ppm peak of DSS (molecule entry 3).
    As[3].αs[2][1];
    As[3].Ωs[2][1];
    Bs[3].ss_params.phase.cos_β[2][1];
    Bs[3].ss_params.phase.sin_β[2][1];
]
fr = xx->NMRSignalSimulator.evalabsorptionlorentzian(xx[1], xx[2], p_test...)
fi = xx->NMRSignalSimulator.evaldispersionlorentzian(xx[1], xx[2], p_test...)
fr = xx->MSS.srs[3][2][1](xx[1],xx[2])
fi = xx->MSS.sis[3][2][1](xx[1],xx[2])

df_AN_eval = ones(2) .* NaN
dfr = xx->NMRSignalSimulator.evalabsorptionlorentzianderivatives!(df_AN_eval, xx[1], xx[2], p_test...)
dfi = xx->NMRSignalSimulator.evaldispersionlorentzianderivatives!(df_AN_eval, xx[1], xx[2], p_test...)
dfr = xx->MSS.∇srs![3][2][1](df_AN_eval, xx[1],xx[2])
dfi = xx->MSS.∇sis![3][2][1](df_AN_eval, xx[1],xx[2])

ND_accuracy_order = 8
dfr_ND = pp->FiniteDifferences.grad(FiniteDifferences.central_fdm(ND_accuracy_order, 1), fr, pp)[1]
dfi_ND = pp->FiniteDifferences.grad(FiniteDifferences.central_fdm(ND_accuracy_order, 1), fi, pp)[1]

# real.
dfr_ND_eval = dfr_ND(x_test)
dfr(x_test)

table_re = Table(AN = df_AN_eval, ND = dfr_ND_eval)
display(table_re)

# imag.
dfi_ND_eval = dfi_ND(x_test)
dfi(x_test)

table_im = Table(AN = df_AN_eval, ND = dfi_ND_eval)
display(table_im)
