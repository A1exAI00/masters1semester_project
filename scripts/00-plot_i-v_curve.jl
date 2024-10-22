using DrWatson
@quickactivate "masters_julia_project"

using Plots

include(srcdir("narahara_papers.jl"))
using .narahara_papers

#########################################################################################

V_min, V_max, V_N =0.0, 1.2, 200
V_range = range(V_min, V_max, V_N)

# I_V_curve = narahara_papers.I_D_1.(V_range)
I_V_curve = narahara_papers.I_D_2.(V_range)

#########################################################################################

p1 = plot(title="Diode I-V curve",
    titlefontsize=12,
    xlabel="V, В", ylabel="I, A", 
    label=nothing,
    xlims=(V_min, V_max),
    ylims=(minimum(I_V_curve), maximum(I_V_curve)),
    minorgrid=true
)
plot!(p1, V_range, I_V_curve, ls=:solid, label="")
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

p = plot(p1, legend=false, size=(500,500), dpi=300, left_margin=(20,:px))

savedir = joinpath("tmp", "00-plot_i-v_curve")
savepath = plotsdir(savedir, "00-plot_i-v_curve_$(time_ns()).pdf")
savefig(p, savepath)