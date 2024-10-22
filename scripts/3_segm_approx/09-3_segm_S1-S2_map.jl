using DrWatson
@quickactivate "masters1semester_project"

using Plots
using NonlinearSolve

include(srcdir("NEC1S1763_3segm_approx.jl"))
using .NEC1S1763_3segm_approx

#########################################################################################

u_start, u_end, u_N = 0.0, 5.0, 200
u_range = range(u_start, u_end, u_N)

v_curve_1 = NEC1S1763_3segm_approx.BAX_3segm_approx2.(u_range)
v_curve_2 = NEC1S1763_3segm_approx.load_line_function.(u_range)

t_start, t_end, t_N = 0.0, 0.3, 100
t_span = (t_start, t_end)
t_range = range(t_start, t_end, t_N)

v₀_start, v₀_end, v₀_N = 1.0, 1.5, 10
v₀_range = range(v₀_start, v₀_end, v₀_N)

u₀ = 1.0

#########################################################################################

function L₂_sol_integration(v₀)
    sol = NEC1S1763_3segm_approx.all_obl_integrate([u₀, v₀], t_span, 2; saveat=t_range)
    return sol
end

sols_int = L₂_sol_integration.(v₀_range)

#########################################################################################

function L₂_sol_by_formula(v₀)
    u(t) = NEC1S1763_3segm_approx.obl2.u(t,v₀)
    v(t) = NEC1S1763_3segm_approx.obl2.v(t,v₀)
    return (u.(t_range), v.(t_range))
end

sols_form = L₂_sol_by_formula.(v₀_range)

#########################################################################################

function L₂_map_by_formula(v₀)
    u(t) = NEC1S1763_3segm_approx.obl2.u(t,v₀)
    v(t) = NEC1S1763_3segm_approx.obl2.v(t,v₀)

    f(τ,p) = u(τ) - NEC1S1763_3segm_approx.obl2.u_max
    sol = solve(IntervalNonlinearProblem(f, t_span))
    τ₁ = sol.u
    return v(τ₁)
end

ends = L₂_map_by_formula.(v₀_range)

#########################################################################################

p1 = plot(title="Траектории L₂ для разных НУ",
    titlefontsize=12,
    xlabel="u", ylabel="v", 
    label=nothing,
    xlims=(0.0, 4.8),
    ylims=(-0.7, 1.5),
    minorgrid=true
)
plot!(p1, u_range, v_curve_1, ls=:dash, label="")
plot!(p1, u_range, v_curve_2, ls=:dash, label="")

for sol in sols_int
    plot!(p1, sol[1,:], sol[2,:], ls=:solid, color=:black, label="", alpha=0.3)
end
# for sol in sols_form
#     plot!(p1, sol[1], sol[2], ls=:solid, color=:red, label="", alpha=0.3)
# end

tmp_a = NEC1S1763_3segm_approx.V_v/NEC1S1763_3segm_approx.V_p
plot!(p1, tmp_a.*ones(length(ends)), ends, label="", seriestype=:scatter, ms=2)

hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

tmp_a = NEC1S1763_3segm_approx.V_v/NEC1S1763_3segm_approx.V_p
tmp_b = NEC1S1763_3segm_approx.obl2.u₀_st
tmp_c = NEC1S1763_3segm_approx.obl2.v₀_st
vline!(p1, [1.0, tmp_a], ls=:dot, label=nothing)
plot!(p1, [tmp_b], [tmp_c], label="", seriestype=:scatter, ms=2)


p = plot(p1, legend=false, size=(500,400), dpi=300)

savedir = joinpath("tmp", "09-3_segm_S1-S2_map")
savepath = plotsdir(savedir, "09-3_segm_S1-S2_map_$(time_ns()).pdf")
savefig(p, savepath)