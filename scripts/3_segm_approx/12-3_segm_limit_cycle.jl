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

t_start, t_end, t_N = 0.0, 0.3, 200
t_span = (t_start, t_end)
t_range = range(t_start, t_end, t_N)

u_S₀_S₁ = 1
u_S₂_S₃ = NEC1S1763_3segm_approx.V_v/NEC1S1763_3segm_approx.V_p

u₁(C₁_1, C₂_1, t) = NEC1S1763_3segm_approx.obl1.u(C₁_1, C₂_1, t)
v₁(C₁_1, C₂_1, t) = NEC1S1763_3segm_approx.obl1.v(C₁_1, C₂_1, t)
u₂(C₁_2, C₂_2, t) = NEC1S1763_3segm_approx.obl2.u(C₁_2, C₂_2, t)
v₂(C₁_2, C₂_2, t) = NEC1S1763_3segm_approx.obl2.v(C₁_2, C₂_2, t)
u₃(C₁_3, C₂_3, t) = NEC1S1763_3segm_approx.obl3.u(C₁_3, C₂_3, t)
v₃(C₁_3, C₂_3, t) = NEC1S1763_3segm_approx.obl3.v(C₁_3, C₂_3, t)
u₄(C₁_4, C₂_4, t) = NEC1S1763_3segm_approx.obl4.u(C₁_4, C₂_4, t)
v₄(C₁_4, C₂_4, t) = NEC1S1763_3segm_approx.obl4.v(C₁_4, C₂_4, t)

#########################################################################################

function f_err(U, p)
    C₁_1, C₂_1, C₁_2, C₂_2, C₁_3, C₂_3, C₁_4, C₂_4, τ₁, τ₂, τ₃, τ₄ = U

    err01 = u₁(C₁_1, C₂_1, 0) - u₄(C₁_4, C₂_4, τ₄)
    err02 = u₄(C₁_4, C₂_4, τ₄) - u_S₀_S₁
    err03 = v₁(C₁_1, C₂_1, 0) - v₄(C₁_4, C₂_4, τ₄)

    err04 = u₁(C₁_1, C₂_1, τ₁) - u₂(C₁_2, C₂_2, 0)
    err05 = u₂(C₁_2, C₂_2, 0) - u_S₀_S₁
    err06 = v₁(C₁_1, C₂_1, τ₁) - v₂(C₁_2, C₂_2, 0)

    err07 = u₂(C₁_2, C₂_2, τ₂) - u₃(C₁_3, C₂_3, 0)
    err08 = u₃(C₁_3, C₂_3, 0) - u_S₂_S₃
    err09 = v₂(C₁_2, C₂_2, τ₂) - v₃(C₁_3, C₂_3, 0)

    err10 = u₃(C₁_3, C₂_3, τ₃) - u₄(C₁_4, C₂_4, 0)
    err11 = u₄(C₁_4, C₂_4, 0) - u_S₂_S₃
    err12 = v₃(C₁_3, C₂_3, τ₃) - v₄(C₁_4, C₂_4, 0)

    return [err01, err02, err03, err04, err05, err06, err07, err08, err09, err10, err11, err12]
end

function segm_3_solve_for_limit_cycle()
    U₀ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1]
    U₀[1] = NEC1S1763_3segm_approx.obl1.C₁(-0.5)
    U₀[2] = NEC1S1763_3segm_approx.obl1.C₂(-0.5)
    U₀[3] = NEC1S1763_3segm_approx.obl2.C₁(1.1)
    U₀[4] = NEC1S1763_3segm_approx.obl2.C₂(1.1)
    U₀[5] = NEC1S1763_3segm_approx.obl3.C₁(1.0)
    U₀[6] = NEC1S1763_3segm_approx.obl3.C₂(1.0)
    U₀[7] = NEC1S1763_3segm_approx.obl4.C₁(0.0)
    U₀[8] = NEC1S1763_3segm_approx.obl4.C₂(0.0)
    prob = NonlinearProblem(f_err, U₀, nothing)
    sol = solve(prob)
    return sol
end

U_sol = segm_3_solve_for_limit_cycle()

#########################################################################################

println(U_sol)
println(U_sol.retcode)

C₁_1, C₂_1, C₁_2, C₂_2, C₁_3, C₂_3, C₁_4, C₂_4, τ₁, τ₂, τ₃, τ₄ = U_sol

u₁_sol_func(t) = u₁(C₁_1, C₂_1, t)
v₁_sol_func(t) = v₁(C₁_1, C₂_1, t)
u₁_sol = u₁_sol_func.(range(0.0, τ₁, 100))
v₁_sol = v₁_sol_func.(range(0.0, τ₁, 100))
println(u₁_sol_func(0.0), v₁_sol_func(0.0))

u₂_sol_func(t) = u₂(C₁_2, C₂_2, t)
v₂_sol_func(t) = v₂(C₁_2, C₂_2, t)
u₂_sol = u₂_sol_func.(range(0.0, τ₂, 100))
v₂_sol = v₂_sol_func.(range(0.0, τ₂, 100))

u₃_sol_func(t) = u₃(C₁_3, C₂_3, t)
v₃_sol_func(t) = v₃(C₁_3, C₂_3, t)
u₃_sol = u₃_sol_func.(range(0.0, τ₃, 100))
v₃_sol = v₃_sol_func.(range(0.0, τ₃, 100))

u₄_sol_func(t) = u₄(C₁_4, C₂_4, t)
v₄_sol_func(t) = v₄(C₁_4, C₂_4, t)
u₄_sol = u₄_sol_func.(range(0.0, τ₄, 100))
v₄_sol = v₄_sol_func.(range(0.0, τ₄, 100))

#########################################################################################

p1 = plot(title="Предельный цикл в системе с кусочно-линейной аппроксимацией",
    titlefontsize=8,
    xlabel="u", ylabel="v", 
    label=nothing,
    xlims=(0.0, 5.0),
    ylims=(-1.0, 1.5),
    minorgrid=true
)
plot!(p1, u_range, v_curve_1, ls=:dash, label="")
plot!(p1, u_range, v_curve_2, ls=:dash, label="")

plot!(p1, u₁_sol, v₁_sol, ls=:solid, color=:black, label="")
plot!(p1, u₂_sol, v₂_sol, ls=:solid, color=:black, label="")
plot!(p1, u₃_sol, v₃_sol, ls=:solid, color=:black, label="")
plot!(p1, u₄_sol, v₄_sol, ls=:solid, color=:black, label="")

hline!(p1, [0.0], color=:black, label="")
vline!(p1, [0.0], color=:black, label="")

vline!(p1, [1.0, NEC1S1763_3segm_approx.V_v/NEC1S1763_3segm_approx.V_p], ls=:dot, label="")

p = plot(p1, legend=false, size=(500,400), dpi=300, 
    left_margin=(20,:px)
)

savedir = joinpath("tmp", "12-3_segm_limit_cycle")
savepath = plotsdir(savedir, "12-3_segm_limit_cycle_$(time_ns()).pdf")
savefig(p, savepath)