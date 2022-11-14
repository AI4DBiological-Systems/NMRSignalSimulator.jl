

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i in eachindex(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i in eachindex(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end



###### forward model and derivatives.


# This should be the same as evalclproxymixture()
function forwardmodelmolecule(
    u_rad::T,
    As,
    Bs,
    ) where T <: AbstractFloat

    out = zero(Complex{T})
    for n in eachindex(Bs)
        A = As[n]
        B = Bs[n]

        # d = x.d
        # κs_λ = x.κs_λ

        # for i in eachindex(Bs[n].qs)
            
        #     r = u_rad - d[i]
            
        #     for k in eachindex(Bs[n].qs[i])
        #         out += qs[i][k](r, κs_λ[i])
        #     end
        # end

        out_sys = NMRSignalSimulator.evalclproxysys(B.qs, u_rad, B.ss_params)

        out_singlets = NMRSignalSimulator.evalclsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
        B.β_singlets, B.λ0, B.κs_λ_singlets)

        out += out_sys + out_singlets
    end

    return out
end


##### debug derivatives.

function loadalphaglucose()

    λ0 = 3.4
    α = [0.014473202716171404, 0.014432448018050659, 0.01579127786786877, 0.014472514849190824, 0.01531313643938015, 0.014419856076392434, 0.014471791589693256, 0.015839333899974762, 0.014431729118944636, 0.015261950218007846, 0.014379194597565995, 0.01579064133843473, 0.016751940468504804, 0.015733196128264886, 0.010100854477667377, 0.004394087079550436, 0.004792753853127242, 0.010576465257294589, 0.01538704558664968, 0.01443396034708299, 0.0157926161216909, 0.014471102326732594, 0.01531125077577049, 0.01441844796278814, 0.01583873329241979, 0.01681169733577139, 0.01578112016287714, 0.010071817434038661, 0.004382096746769656, 0.004776065212253049, 0.010541880185742819, 0.015335545043430338, 0.011040232404787645, 0.004775840844333532, 0.005229233541496932, 0.011583880244536644, 0.01683339972548277, 0.015391346498934396, 0.01583785239730858, 0.01443324285521636, 0.015263948975065915, 0.014380703424819072, 0.015791981061200118, 0.016753768222180207, 0.01573453155978096, 0.010099840343126746, 0.004393677978786782, 0.004792131179742463, 0.010575199009286542, 0.015385146645125953, 0.011073567921199698, 0.004790867494886511, 0.0052495892571826405, 0.011623418092208232, 0.01689354241442127, 0.015339817988118296, 0.016838314946364445, 0.01583725044888322, 0.016809714886427123, 0.015779642015283257, 0.010072899981490453, 0.0043825389304573345, 0.00477672662682473, 0.010543221285965777, 0.015337557773855247, 0.01104119515792097, 0.004776228800015373, 0.005229841299140748, 0.011585103316711068, 0.01683524082067617, 0.015389446201314603, 0.016898493320779317, 0.011072497807161895, 0.004790443451139437, 0.0052489405339315445, 0.011622080811992145, 0.016891545903554923, 0.015341832142907803, 0.016840157426352444, 0.016896495360376185]
    Ω = [81848.47957380593, 81848.59881351044, 81909.66578558966, 81848.7366115856, 81906.53233950376, 81848.86350727621, 81848.47961032146, 81909.78502085072, 81848.85586139376, 81906.65433442319, 81848.98279120357, 81909.92325162838, 81967.71818214105, 81910.05050211935, 81855.98212039226, 81836.54149598608, 81919.22636661968, 81899.78574221351, 81907.05819914796, 81848.60444592091, 81909.6602230994, 81848.7366481428, 81906.53238884089, 81848.8635439818, 81910.04250361935, 81967.84019844834, 81910.16974978999, 81856.10127037951, 81836.66061804123, 81919.34853896323, 81899.90788662495, 81907.18024256936, 81917.14860372478, 81897.70413289768, 81980.43747714229, 81960.9930063152, 81968.24475788718, 81907.17317738908, 81909.78505425296, 81848.86149391747, 81906.65997798476, 81848.98842372393, 81909.91768915378, 81967.71263297307, 81910.04493972709, 81855.98215722728, 81836.54153282104, 81919.22641587292, 81899.78579146668, 81907.05824870011, 81917.268985439, 81897.82472915802, 81980.558169673, 81961.11391339201, 81968.36679514011, 81907.29522598434, 81968.36021782777, 81910.04253710088, 81967.84024349123, 81910.16978322496, 81856.1069058871, 81836.66625401804, 81919.35417926768, 81899.91352739862, 81907.18588627686, 81917.14304445099, 81897.69857419084, 81980.43192426773, 81960.98745400757, 81968.23920884836, 81907.17322694202, 81968.48226288757, 81917.26902425931, 81897.8247689235, 81980.55820850337, 81961.11395316756, 81968.36684025562, 81907.30086977014, 81968.35466874571, 81968.48230803385]
    A_r = 78757.19441619715:1.0:93018.19441619715
    A_ξ = 0.5:0.05:2.5

    Interpolations = NMRSignalSimulator.Interpolations

    f = (rr,ξξ)->NMRSignalSimulator.evalclpart(rr, α, Ω, ξξ*λ0)
    #f = (rr,ξξ)->evalclpart(rr, α, Ω, λ[ind])
    A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]

    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_ξ)
    real_setp = Interpolations.extrapolate(real_sitp, 0.0) # zero outside interp range.

    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_ξ)
    imag_setp = Interpolations.extrapolate(imag_sitp, 0.0) # zero outside interp range.

    return real_setp, imag_setp, f
end


