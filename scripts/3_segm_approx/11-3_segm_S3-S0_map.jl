using DrWatson
@quickactivate "masters_julia_project"

using Plots
using NonlinearSolve

include(srcdir("segm_3_model.jl"))
using .segm_3_model

#########################################################################################

u_start, u_end, u_N = 0.0, 5.0, 200
u_range = range(u_start, u_end, u_N)

v_curve_1 = segm_3_model.f.(u_range)
v_curve_2 = segm_3_model.g.(u_range)

t_start, t_end, t_N = 0.0, 0.3, 200
t_span = (t_start, t_end)
t_range = range(t_start, t_end, t_N)

v₀_start, v₀_end, v₀_N = -0.25, segm_3_model.I_v./segm_3_model.I_p, 10
v₀_range = range(v₀_start, v₀_end, v₀_N)

u₀ = segm_3_model.V_v/segm_3_model.V_p

#########################################################################################

function L₄_sol_integration(v₀)
    sol = segm_3_model.all_obl_integrate([u₀, v₀], t_span, 4; saveat=t_range)
    return sol
end
sols_int = L₄_sol_integration.(v₀_range)


function L₄_sol_by_formula(v₀)
    u(t) = segm_3_model.u₄(t,v₀)
    v(t) = segm_3_model.v₄(t,v₀)
    return (u.(t_range), v.(t_range))
end
sols_form = L₄_sol_by_formula.(v₀_range)


function L₄_map_by_formula(v₀)
    u(t) = segm_3_model.u₄(t,v₀)
    v(t) = segm_3_model.v₄(t,v₀)

    f(τ,p) = u(τ) - 1
    sol = solve(IntervalNonlinearProblem(f, (0.0, 0.3)))
    τ₁ = sol.u
    return v(τ₁)
end
ends = L₄_map_by_formula.(v₀_range)

#########################################################################################

p1 = plot(title="Траектории L₄ для разных НУ",
    titlefontsize=12,
    xlabel="u", ylabel="v", 
    label=nothing,
    xlims=(0.0, 5.0),
    ylims=(-1.0, 1.5),
    minorgrid=true
)
plot!(p1, u_range, v_curve_1, ls=:dash, label="")
plot!(p1, u_range, v_curve_2, ls=:dash, label="")

for sol in sols_int
    plot!(p1, sol[1,:], sol[2,:], ls=:solid, color=:black, label="", alpha=0.3)
end
for sol in sols_form
    plot!(p1, sol[1], sol[2], ls=:solid, color=:red, label="", alpha=0.3)
end

plot!(p1, ones(length(ends)), ends, label="", seriestype=:scatter, ms=2)

hline!(p1, [0.0], color=:black, label="")
vline!(p1, [0.0], color=:black, label="")

vline!(p1, [1.0, segm_3_model.V_v/segm_3_model.V_p], ls=:dot, label="")
plot!(p1, [segm_3_model.u₀_st_2], [segm_3_model.v₀_st_2], label="", seriestype=:scatter, ms=2)

p = plot(p1, legend=false, size=(500,400), dpi=300, 
    left_margin=(20,:px)
)

savedir = joinpath("tmp", "11-3_segm_S3-S0_map")
savepath = plotsdir(savedir, "11-3_segm_S3-S0_map_$(time_ns()).pdf")
savefig(p, savepath)