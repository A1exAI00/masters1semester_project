using DrWatson
@quickactivate "masters1semester_project"

using Plots
# using NonlinearSolve
using Optim
using DifferentialEquations

include(srcdir("narahara_papers.jl"))
using .narahara_papers

#########################################################################################

# as = [a₀, ..., aₙ]
function general_poly(x, as)
    res = 0.0
    for i in 0:length(as)-1
        res += as[i+1]*x^i
    end
    return res 
end

function general_poly_approx_loss(params, x_data, y_data, func)
    func_(x) = func(x, params)
    sq_err = (y_data .- func_.(x_data)).^2
    return sum(sq_err)/length(sq_err)
end

function general_poly_approx_optimize(loss, params₀)
    loss_(params) = loss(params, narahara_papers.IV_curve_V_data, narahara_papers.IV_curve_I_data)
    options = Optim.Options(g_tol=0.0, iterations=10000,)
    optim_result = optimize(loss_, params₀, LBFGS(), options)
    display(optim_result)
    return Optim.minimizer(optim_result)
end

#########################################################################################

function general_model(approx, dU, U, p, t)
    as, δ, E, ε = p; u, v = U
    dU[1] = 1/ε * (v - approx(u, as))
    dU[2] = E - v - δ*u
    return nothing
end

function general_integrate(model, U₀, t_span, as;
    abstol=nothing, reltol=nothing, alg=nothing, saveat=nothing)
    ABSTOL = 1e-5; abstol = if isnothing(abstol) ABSTOL end
    RELTOL = 1e-5; reltol = if isnothing(reltol) RELTOL end
    ALG = Tsit5(); alg = if isnothing(alg) ALG end

    prob = ODEProblem(model, U₀, t_span, (as, δ, E, ε))
    sol = if isnothing(saveat) 
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL) else
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, saveat=saveat) end
    return sol
end

#########################################################################################

for degree in 5:13

    loss(params, x_data, y_data) = general_poly_approx_loss(params, x_data, y_data, general_poly)
    approx_optimize(params₀) = general_poly_approx_optimize(loss, params₀)
    optim_params = approx_optimize(zeros(degree))
    poly3approx_optim(x) = general_poly(x, optim_params)


    p1 = plot(title="NEC1S1763 I-V curve",
        titlefontsize=12,
        xlabel="V, В", ylabel="I, A", 
        label=nothing,
        xlims=(minimum(narahara_papers.IV_curve_V_data), maximum(narahara_papers.IV_curve_V_data.*1.1)),
        ylims=(minimum(narahara_papers.IV_curve_I_data), maximum(narahara_papers.IV_curve_I_data.*1.1)),
        minorgrid=true
    )

    # Координатные оси
    hline!(p1, [0.0], color=:black, label=nothing)
    vline!(p1, [0.0], color=:black, label=nothing)

    # Экспериментальные данные
    plot!(p1, narahara_papers.IV_curve_V_data, narahara_papers.IV_curve_I_data, ls=:solid, label="")

    # Аппроксимации
    plot!(p1, poly3approx_optim, ls=:solid, label="$(degree)-я степень")

    p = plot(p1, legend=false, size=(500,500), dpi=300, left_margin=(20,:px))

    savedir = joinpath("tmp", "24-poly_approx")
    savepath = plotsdir(savedir, "24_$(lpad(degree,3,"0"))-poly_approx_$(time_ns()).pdf")
    savefig(p, savepath)
end


