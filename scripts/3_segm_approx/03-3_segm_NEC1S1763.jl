using DrWatson
@quickactivate "masters1semester_project"

using Plots
using Optim

include(srcdir("narahara_papers.jl"))
using .narahara_papers

include(srcdir("NEC1S1763_3segm_approx.jl"))
using .NEC1S1763_3segm_approx

#########################################################################################

V_start, V_end, V_N = 0.0, 0.58, 200
target_V_data = range(V_start, V_end, V_N)

target_I_data = narahara_papers.I_D_1.(target_V_data)

#########################################################################################

piecewise_linear_3_func = NEC1S1763_3segm_approx.general_3segm_approx

function piecewise_linear_3_loss(params)
    func(x) = piecewise_linear_3_func(x, params)
    err = target_I_data .- func.(target_V_data)
    sq_err = err .* err
    return sum(sq_err)/length(sq_err)
end

#########################################################################################

# xa, xb, ya, yb, k
piecewise_linear_3_params_0 = [65e-3, 0.3, 6.3e-3, 1.575e-3, 0.1]

# optim_result = optimize(piecewise_linear_3_loss, piecewise_linear_3_params_0, BFGS())
optim_result = optimize(piecewise_linear_3_loss, piecewise_linear_3_params_0, LBFGS())
optim_params = Optim.minimizer(optim_result)
println(optim_result)
println(optim_params)

optim_params = [0.1, 0.4, 0.0058, 0.001, 0.025]

func(x) = piecewise_linear_3_func(x, optim_params)
I_approx = func.(target_V_data)

#########################################################################################

p1 = plot(title="Piecewise-linear func. approx. for NEC 1S1763",
    titlefontsize=10,
    xlabel="V, В", ylabel="I, A", 
    label=nothing,
    # xlims=,
    # ylims=,
    minorgrid=true
)
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

plot!(p1, target_V_data, target_I_data, ls=:dash, label="Original")
plot!(p1, target_V_data, I_approx, ls=:solid, label="Approx.")

plot!(p1, 
    xlims=(V_start, V_end), 
    ylims=(0.0, max(maximum(target_I_data), maximum(I_approx)))
)


p = plot(p1, legend=true, size=(500,500), dpi=300, 
    # left_margin=(20,:px)
)

savedir = joinpath("tmp", "03-3_segm_NEC1S1763")
savepath = plotsdir(savedir, "03-3_segm_NEC1S1763_$(time_ns()).pdf")
savefig(p, savepath)