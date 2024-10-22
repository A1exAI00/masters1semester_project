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

u_S₀_S₁ = 1
u_S₂_S₃ = segm_3_model.V_v/segm_3_model.V_p

#########################################################################################

function f_err(U, p)
    C₁_1, C₂_1, C₁_2, C₂_2, C₁_3, C₂_3, C₁_4, C₂_4, τ₁, τ₂, τ₃, τ₄ = U

    err01 = segm_3_model.u₁(C₁_1, C₂_1, 0) - segm_3_model.u₄(C₁_4, C₂_4, τ₄)
    err02 = segm_3_model.u₄(C₁_4, C₂_4, τ₄) - u_S₀_S₁
    err03 = segm_3_model.v₁(C₁_1, C₂_1, 0) - segm_3_model.v₄(C₁_4, C₂_4, τ₄)

    err04 = segm_3_model.u₁(C₁_1, C₂_1, τ₁) - segm_3_model.u₂(C₁_2, C₂_2, 0)
    err05 = segm_3_model.u₂(C₁_2, C₂_2, 0) - u_S₀_S₁
    err06 = segm_3_model.v₁(C₁_1, C₂_1, τ₁) - segm_3_model.v₂(C₁_2, C₂_2, 0)

    err07 = segm_3_model.u₂(C₁_2, C₂_2, τ₂) - segm_3_model.u₃(C₁_3, C₂_3, 0)
    err08 = segm_3_model.u₃(C₁_3, C₂_3, 0) - u_S₂_S₃
    err09 = segm_3_model.v₂(C₁_2, C₂_2, τ₂) - segm_3_model.v₃(C₁_3, C₂_3, 0)

    err10 = segm_3_model.u₃(C₁_3, C₂_3, τ₃) - segm_3_model.u₄(C₁_4, C₂_4, 0)
    err11 = segm_3_model.u₄(C₁_4, C₂_4, 0) - u_S₂_S₃
    err12 = segm_3_model.v₃(C₁_3, C₂_3, τ₃) - segm_3_model.v₄(C₁_4, C₂_4, 0)

    # err_vec = [err01, err02, err03, err04, err05, err06, err07, err08, err09, err10, err11, err12]
    # fin_err = sum(err_vec)
    # return fin_err

    return [err01, err02, err03, err04, err05, err06, err07, err08, err09, err10, err11, err12]
end

function segm_3_solve_for_limit_cycle()
    U₀ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1]
    U₀[1] = segm_3_model.C₁_1(-0.5)
    U₀[2] = segm_3_model.C₂_1(-0.5)
    U₀[3] = segm_3_model.C₁_2(1.1)
    U₀[4] = segm_3_model.C₂_2(1.1)
    U₀[5] = segm_3_model.C₁_3(1.0)
    U₀[6] = segm_3_model.C₂_3(1.0)
    U₀[7] = segm_3_model.C₁_4(0.0)
    U₀[8] = segm_3_model.C₂_4(0.0)
    prob = NonlinearProblem(f_err, U₀, nothing)
    sol = solve(prob)
    return sol
end

U_sol = segm_3_solve_for_limit_cycle()


function segm_3_solve_for_limit_cycle_start_random()
    C_min, C_max = -2.0, 2.0
    success_sol = nothing
    i = 1
    while true
        println("Try #", i)
        C₀_random = rand(8).*(C_max-C_min) .- C_min
        U₀_random = [C₀_random..., 0.2, 0.4, 0.6, 0.8]
        prob = NonlinearProblem(f_err, U₀_random, nothing)
        sol = solve(prob)
        # println(sol.retcode)
        if SciMLBase.successful_retcode(sol)
            success_sol = sol
            break
        end
        i += 1
    end
    return success_sol
end

# U_sol = segm_3_solve_for_limit_cycle_start_random()

#########################################################################################

println(U_sol)
println(U_sol.retcode)

C₁_1, C₂_1, C₁_2, C₂_2, C₁_3, C₂_3, C₁_4, C₂_4, τ₁, τ₂, τ₃, τ₄ = U_sol

u₁_sol_func(t) = segm_3_model.u₁(C₁_1, C₂_1, t)
v₁_sol_func(t) = segm_3_model.v₁(C₁_1, C₂_1, t)
u₁_sol = u₁_sol_func.(range(0.0, τ₁, 100))
v₁_sol = v₁_sol_func.(range(0.0, τ₁, 100))
println(u₁_sol_func(0.0), v₁_sol_func(0.0))

u₂_sol_func(t) = segm_3_model.u₂(C₁_2, C₂_2, t)
v₂_sol_func(t) = segm_3_model.v₂(C₁_2, C₂_2, t)
u₂_sol = u₂_sol_func.(range(0.0, τ₂, 100))
v₂_sol = v₂_sol_func.(range(0.0, τ₂, 100))

u₃_sol_func(t) = segm_3_model.u₃(C₁_3, C₂_3, t)
v₃_sol_func(t) = segm_3_model.v₃(C₁_3, C₂_3, t)
u₃_sol = u₃_sol_func.(range(0.0, τ₃, 100))
v₃_sol = v₃_sol_func.(range(0.0, τ₃, 100))

u₄_sol_func(t) = segm_3_model.u₄(C₁_4, C₂_4, t)
v₄_sol_func(t) = segm_3_model.v₄(C₁_4, C₂_4, t)
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

vline!(p1, [1.0, segm_3_model.V_v/segm_3_model.V_p], ls=:dot, label="")

p = plot(p1, legend=false, size=(500,400), dpi=300, 
    left_margin=(20,:px)
)

savedir = joinpath("tmp", "12-3_segm_limit_cycle")
savepath = plotsdir(savedir, "12-3_segm_limit_cycle_$(time_ns()).pdf")
savefig(p, savepath)