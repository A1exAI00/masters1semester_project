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
t_span = (t_start, t_end)
t_range = range(t_start, t_end, t_N)

v₀_start, v₀_end, v₀_N = NEC1S1763_5segm_approx.I_2/NEC1S1763_5segm_approx.I_1, 1.7, 10
v₀_range = range(v₀_start, v₀_end, v₀_N)

#########################################################################################

obl4 = NEC1S1763_5segm_approx.create_obl_from_obl_N(4, δ, E, ε)

NEC1S1763_5segm_approx.display_params(obl4)
# println(obl4.is_λ_real)
# println((obl4.u_st, obl4.v_st))
# println(obl4.C₁(1.0,0.5))
# println(obl4.C₂(1.0,0.5))

u₀ = obl4.u_min_bound
u_L_end = obl4.u_max_bound

#########################################################################################

function L_calc_τ(v₀)
    u(t) = obl4.sol_u(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    v(t) = obl4.sol_v(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    f(τ,p) = u(τ) - u_L_end
    sol = solve(IntervalNonlinearProblem(f, (0.0, 0.1)))
    τ₁ = sol.u
    return τ₁
end

function L_map_by_formula(v₀)
    u(t) = obl4.sol_u(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    v(t) = obl4.sol_v(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(v₀)
    return v(τ₁)
end

function L_sol_by_formula(v₀)
    u(t) = obl4.sol_u(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    v(t) = obl4.sol_v(obl4.C₁(u₀,v₀), obl4.C₂(u₀,v₀), t)
    τ₁ = L_calc_τ(v₀)
    t_range = range(0.0, τ₁, 100)
    return (u.(t_range), v.(t_range))
end

τs = L_calc_τ.(v₀_range)
Ls = L_sol_by_formula.(v₀_range)
ends = L_map_by_formula.(v₀_range)

#########################################################################################

p1 = plot(title="Траектории L₄ для разных НУ",
    titlefontsize=12,
    xlabel="u", ylabel="v", 
    label=nothing,
    xlims=(u_start, u_end),
    ylims=(-0.7, 1.7),
    minorgrid=true
)

# ВАХ и нагрузочная прямая
plot!(p1, u_range, v_curve_1, ls=:dash, label="")
plot!(p1, u_range, v_curve_2, ls=:dash, label="")

# Траектории путем интегрирования
# for sol in sols_int
#     plot!(p1, sol[1,:], sol[2,:], ls=:solid, color=:black, label="", alpha=0.3)
# end

# Траектории из формулы
for L in Ls
    Lx, Ly = L
    plot!(p1, Lx, Ly, ls=:solid, color=:red, label="", alpha=0.3)
end

# Концы траекторий
plot!(p1, u_L_end*ones(length(ends)), ends, label="", seriestype=:scatter, ms=2)

# Координатные оси
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# Границы областей и состояние равновесия
border1 = 1.0
border2 = NEC1S1763_5segm_approx.V_2/NEC1S1763_5segm_approx.V_1
border3 = NEC1S1763_5segm_approx.V_3/NEC1S1763_5segm_approx.V_1
border4 = NEC1S1763_5segm_approx.V_4/NEC1S1763_5segm_approx.V_1
borders = [border1, border2, border3, border4]
vline!(p1, borders, ls=:dot, label=nothing)
plot!(p1, [obl4.u_st], [obl4.v_st], seriestype=:scatter, label=nothing)


p = plot(p1, legend=false, size=(500,400), dpi=300)

savedir = joinpath("tmp", "17-5_segm_S3-S4_map")
savepath = plotsdir(savedir, "17-5_segm_S3-S4_map_$(time_ns()).pdf")
savefig(p, savepath)