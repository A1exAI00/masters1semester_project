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

u₀ = NEC1S1763_5segm_approx.u_1
v₀_span = (-0.5, NEC1S1763_5segm_approx.v_1)
v₀_range = range(v₀_span[1], v₀_span[2], 200)

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
    (0.001, 0.6),
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
    if isnan(v₀) return (NaN, NaN) end
    u(t) = obl.sol_u(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    v(t) = obl.sol_v(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    if τ₁ == 0.0
        return (NaN, NaN)
    end
    return (u(τ₁), v(τ₁))
end

function L_sol_by_formula(obl, u₀, v₀, u_L_end, t_span)
    u(t) = obl.sol_u(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    v(t) = obl.sol_v(obl.C₁(u₀,v₀), obl.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
    t_range = range(0.0, τ₁, 100)
    return (u.(t_range), v.(t_range))
end

#########################################################################################

function map_interate(v₀)
    u₁, v₁ = L_map_by_formula(obl1, u₀, v₀, obl1.u_max_bound, t_spans[1])
    u₂, v₂ = L_map_by_formula(obl2, u₁, v₁, obl2.u_max_bound, t_spans[2])
    u₃, v₃ = L_map_by_formula(obl3, u₂, v₂, obl3.u_max_bound, t_spans[3])
    u₄, v₄ = L_map_by_formula(obl4, u₃, v₃, obl4.u_max_bound, t_spans[4])
    u₅, v₅ = L_map_by_formula(obl5, u₄, v₄, obl4.u_max_bound, t_spans[5])
    u₆, v₆ = L_map_by_formula(obl4, u₅, v₅, obl4.u_min_bound, t_spans[6])
    u₇, v₇ = L_map_by_formula(obl3, u₆, v₆, obl3.u_min_bound, t_spans[7])
    u₈, v₈ = L_map_by_formula(obl2, u₇, v₇, obl2.u_min_bound, t_spans[8])
    return v₈
end

v_new = map_interate.(v₀_range)

valid_v₀ = []
valid_v_new = []
for i in eachindex(v₀_range)
    if !isnan(v_new[i])
        push!(valid_v₀, v₀_range[i])
        push!(valid_v_new, v_new[i])
    end
end


#########################################################################################

p1 = plot(title="",
    titlefontsize=12,
    xlabel="v", ylabel="v_new", 
    label=nothing,
    xlims=(minimum(valid_v₀), maximum(valid_v₀)),
    ylims=(minimum(valid_v_new), maximum(valid_v_new)),
    minorgrid=true
)

# Координатные оси
# hline!(p1, [0.0], color=:black, label=nothing)
# vline!(p1, [0.0], color=:black, label=nothing)

plot!(p1, valid_v₀, valid_v_new, ls=:solid, color=:blue, label="", alpha=1.0)
# plot!(p1, v₀_range, v_new.-v₀_range, ls=:solid, color=:blue, label="", alpha=1.0)
plot!(p1, [-1, 1], [-1, 1], ls=:dot, color=:black, label="", alpha=1.0)

# ВАХ и нагрузочная прямая
# plot!(p1, u_range, v_curve_1, ls=:dash, label="")
# plot!(p1, u_range, v_curve_2, ls=:dash, label="")

# Границы областей
# borders = [NEC1S1763_5segm_approx.u_1, NEC1S1763_5segm_approx.u_2, NEC1S1763_5segm_approx.u_3, NEC1S1763_5segm_approx.u_4]
# vline!(p1, borders, ls=:dot, label=nothing)

# for i in eachindex(obls)
#     obl = obls[i]
#     v₀_span = v₀_spans[i]
#     v₀_range = range(v₀_span[1], v₀_span[2], 10)
#     u₀, u_L_end = u₀s[i], u_L_ends[i]
#     t_span = t_spans[i]

#     L_calc_τ_(v₀) = L_calc_τ(obl, u₀, v₀, u_L_end, t_span)
#     L_sol_by_formula_(v₀) = L_sol_by_formula(obl, u₀, v₀, u_L_end, t_span)
    
#     τs = L_calc_τ_.(v₀_range); println(i, τs)
#     Ls = L_sol_by_formula_.(v₀_range)
#     # ends = L_map_by_formula.(v₀_range)

#     for L in Ls
#         Lx, Ly = L
#         plot!(p1, Lx, Ly, ls=:solid, color=:red, label="", alpha=0.5)
#     end
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

savedir = joinpath("tmp", "23-5_segm_grand_map")
savepath = plotsdir(savedir, "23-5_segm_grand_map_$(time_ns()).pdf")
savefig(p, savepath)