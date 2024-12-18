using DrWatson
@quickactivate "masters1semester_project"

using Plots
# using NonlinearSolve
# using Optim
# using DifferentialEquations

include(srcdir("narahara_papers.jl"))
using .narahara_papers

include(srcdir("NEC1S1763_poly_approx.jl"))
using .NEC1S1763_poly_approx

include(srcdir("tools.jl"))

#########################################################################################

target_Vs_nonuniform = narahara_papers.IV_curve_V_data
target_Is_nonuniform = narahara_papers.IV_curve_I_data
curr_segm_approx(x) = general_segm_approx_symmetrical(x, target_Vs_nonuniform, target_Is_nonuniform)

target_Vs_uniform = range(0.0, 1.0*maximum(target_Vs_nonuniform), 200)
target_Is_uniform = curr_segm_approx.(target_Vs_uniform)

#########################################################################################

for degree in 5:13
    # Calc oprimized parameters and optimized function
    weights = ones(length(target_Vs_uniform))
    optimized_params = NEC1S1763_poly_approx.easy_poly_optimize(degree, target_Vs_uniform, target_Is_uniform, weights)
    poly_approx_optim(x) = NEC1S1763_poly_approx.general_poly(x, optimized_params)

    p1 = plot(title="NEC1S1763 I-V curve",
        titlefontsize=12,
        xlabel="V, В", ylabel="I, A", 
        label=nothing,
        xlims=(minimum(target_Vs_uniform), maximum(target_Vs_uniform)),
        ylims=(minimum(target_Is_uniform), maximum(target_Is_uniform)),
        minorgrid=true
    )

    # Координатные оси
    hline!(p1, [0.0], color=:black, label=nothing)
    vline!(p1, [0.0], color=:black, label=nothing)

    # Экспериментальные данные
    plot!(p1, target_Vs_uniform, target_Is_uniform, ls=:solid, label="")
    plot!(p1, target_Vs_nonuniform, target_Is_nonuniform, ls=:solid, label="")

    # Аппроксимации
    plot!(p1, target_Vs_uniform, poly_approx_optim.(target_Vs_uniform), 
        ls=:solid, label="$(degree)-я степень"
    )

    p = plot(p1, legend=false, size=(500,500), dpi=300, left_margin=(20,:px))

    savedir = joinpath("tmp", "24-poly_approx")
    savepath = plotsdir(savedir, "24_$(lpad(degree,3,"0"))-poly_approx_$(time_ns()).pdf")
    savefig(p, savepath)
end


