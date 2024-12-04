using DrWatson
@quickactivate "masters1semester_project"

using Plots
using NonlinearSolve

include(srcdir("NEC1S1763_5segm_approx.jl"))
using .NEC1S1763_5segm_approx

#########################################################################################

R = 5 # Ohm
C = 2e-9 # F
L = 5e-6 # H
V_B = 0.2 # V

R₀ = NEC1S1763_5segm_approx.V_1/NEC1S1763_5segm_approx.I_1
E = V_B/NEC1S1763_5segm_approx.I_1/R
δ = R₀/R
ε = C*R₀*R/L

#########################################################################################

u_start, u_end, u_N = 0.0, 6.5, 200
u_range = range(u_start, u_end, u_N)

v_curve_1 = NEC1S1763_5segm_approx.BAX_5segm_approx2.(u_range)
v_curve_2 = NEC1S1763_5segm_approx.load_line_function.(u_range, E, δ)

t_start, t_end, t_N = 0.0, 0.3, 200
t_range = range(t_start, t_end, t_N)

#########################################################################################

v₀_ls = -3.8065

v₀_spans = [
    (-0.45, NEC1S1763_5segm_approx.v_1),
    (NEC1S1763_5segm_approx.v_1, 1.65,),
    (NEC1S1763_5segm_approx.v_2, 1.65,),
    (NEC1S1763_5segm_approx.v_3, 1.65,),
    (NEC1S1763_5segm_approx.v_4, 1.65,),
    (-0.45, NEC1S1763_5segm_approx.v_4),
    (-0.45, NEC1S1763_5segm_approx.v_3),
    (-0.45, NEC1S1763_5segm_approx.v_2),
]
u₀s = [
    NEC1S1763_5segm_approx.u_1,
    NEC1S1763_5segm_approx.u_1,
    NEC1S1763_5segm_approx.u_2,
    NEC1S1763_5segm_approx.u_3,
    NEC1S1763_5segm_approx.u_4,
    NEC1S1763_5segm_approx.u_4,
    NEC1S1763_5segm_approx.u_3,
    NEC1S1763_5segm_approx.u_2
]
u_L_ends = circshift(u₀s, -1)

obl1 = NEC1S1763_5segm_approx.create_obl_from_obl_N(1, δ, E, ε)
obl2 = NEC1S1763_5segm_approx.create_obl_from_obl_N(2, δ, E, ε)
obl3 = NEC1S1763_5segm_approx.create_obl_from_obl_N(3, δ, E, ε)
obl4 = NEC1S1763_5segm_approx.create_obl_from_obl_N(4, δ, E, ε)
obl5 = NEC1S1763_5segm_approx.create_obl_from_obl_N(5, δ, E, ε)
obl6 = NEC1S1763_5segm_approx.create_obl_from_obl_N(6, δ, E, ε)
obl7 = NEC1S1763_5segm_approx.create_obl_from_obl_N(7, δ, E, ε)
obl8 = NEC1S1763_5segm_approx.create_obl_from_obl_N(8, δ, E, ε)

obls = [obl1, obl2, obl3, obl4, obl5, obl6, obl7, obl8]

t_spans = [
    (0.0001, 0.6),
    (0.0, 0.3),
    (0.0, 0.26),
    (0.0, 0.15),
    (0.0001, 0.3),
    (0.0, 0.3),
    (0.0, 0.3),
    (0.0, 0.3)
]

#########################################################################################

function L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    u(t) = obl.sol_u(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    v(t) = obl.sol_v(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    f(τ,p) = u(τ) - u_L_end
    sol = solve(IntervalNonlinearProblem(f, t_span))
    τ₁ = sol.u
    return τ₁
end

function L_map_by_formula(obl, u₀, v₀, u_L_end, t_span)
    u(t) = obl.sol_u(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    v(t) = obl.sol_v(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    return v(τ₁)
end

function L_sol_by_formula(obl, u₀, v₀, u_L_end, t_span)
    u(t) = obl.sol_u(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    v(t) = obl.sol_v(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    t_range = range(0.0, τ₁, 100)
    return (u.(t_range), v.(t_range))
end

#########################################################################################

p1 = plot(title="",
    titlefontsize=12,
    xlabel="u", ylabel="v", 
    label=nothing,
    xlims=(u_start, u_end),
    ylims=(-0.5, 1.7),
    minorgrid=true
)

# Координатные оси
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# ВАХ и нагрузочная прямая
plot!(p1, u_range, v_curve_1, ls=:dash, label="")
plot!(p1, u_range, v_curve_2, ls=:dash, label="")

# Границы областей
borders = [NEC1S1763_5segm_approx.u_1, NEC1S1763_5segm_approx.u_2, NEC1S1763_5segm_approx.u_3, NEC1S1763_5segm_approx.u_4]
vline!(p1, borders, ls=:dot, label=nothing)

for i in eachindex(obls)
    global v₀_ls
    obl = obls[i]
    v₀_span = v₀_spans[i]
    v₀_range = range(v₀_span[1], v₀_span[2], 10)
    u₀, u_L_end = u₀s[i], u_L_ends[i]
    t_span = t_spans[i]

    L_calc_τ_(v₀) = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    L_sol_by_formula_(v₀) = L_sol_by_formula(obl, u₀, v₀, u_L_end, t_span)
    L_map_by_formula_(v₀) = L_map_by_formula(obl, u₀, v₀, u_L_end, t_span)
    
    # τs = L_calc_τ_.(v₀_range)
    Ls = L_sol_by_formula_.(v₀_range)

    for L in Ls
        Lx, Ly = L
        plot!(p1, Lx, Ly, ls=:solid, color=:red, label="", alpha=0.5)
    end

    L_ls_u, L_ls_v = L_sol_by_formula_(v₀_ls)
    v₀_ls = L_ls_v[end]
    plot!(p1, L_ls_u, L_ls_v, ls=:solid, color=:black, label="")
end

# for i in eachindex(L_ls)
#     L_lsx, L_lsy = L_ls[i]
#     plot!(p1, L_lsx, L_lsy, ls=:solid, color=:black, label="", alpha=1.0)
# end
    

# Траектории путем интегрирования
# for sol in sols_int
#     plot!(p1, sol[1,:], sol[2,:], ls=:solid, color=:black, label="", alpha=0.3)
# end

# Траектории из формулы
# for L in Ls
#     Lx, Ly = L
#     plot!(p1, Lx, Ly, ls=:solid, color=:red, label="", alpha=0.3)
# end

# Концы траекторий
# plot!(p1, u_L_end*ones(length(ends)), ends, label="", seriestype=:scatter, ms=2)

# состояние равновесия
# plot!(p1, [obl8.u_st], [obl8.v_st], seriestype=:scatter, label=nothing)


p = plot(p1, legend=false, size=(500,400), dpi=300)

savedir = joinpath("tmp", "22-5_segm_S0-S0_map")
savepath = plotsdir(savedir, "22-5_segm_S0-S0_map_$(time_ns()).pdf")
savefig(p, savepath)