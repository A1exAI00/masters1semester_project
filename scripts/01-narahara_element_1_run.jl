using DrWatson
@quickactivate "masters1semester_project"

using Plots

include(srcdir("narahara_papers.jl"))
using .narahara_papers

#########################################################################################

V_offset = 0.04
I_offset = 0.0005
N_V = 200

#########################################################################################

V_B = 0.2 # V

t_start, t_end = 0.0, 1e-5/4
t_span = (t_start, t_end)

u₀ = [0.0, 0.0]

#########################################################################################

sol = narahara_papers.integrate_narahara_element_1(u₀, t_span, V_B)

sol_V = sol[1,:]
sol_I = sol[2,:]
sol_t = sol.t

max_V, min_V = maximum(sol_V), minimum(sol_V)
max_I, min_I = maximum(sol_I), minimum(sol_I)

V_range = range(min_V-V_offset, max_V+V_offset, N_V)
V_nullcline = narahara_papers.V_nullcline_1.(V_range)
I_nullcline = narahara_papers.I_nullcline_1.(V_range, V_B)

#########################################################################################

p1 = plot(title="V_B=$(V_B)В",
    titlefontsize=10,
    xlabel="V, В", ylabel="I, A", 
    label=nothing,
    xlims=(min_V-V_offset, max_V+V_offset),
    ylims=(min_I-I_offset, max_I+I_offset),
    minorgrid=true
)
plot!(p1, V_range, V_nullcline, ls=:dash, label="V-nullcline")
plot!(p1, V_range, I_nullcline, ls=:dot, label="I-nullcline")
plot!(p1, sol_V, sol_I, label=nothing)
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

p2 = plot(sol_t, sol_V, 
    title="", 
    titlefontsize=10,
    xlabel="t, s", ylabel="V, В", 
    label=nothing,
    minorgrid=true
)
hline!(p2, [0.0], color=:black, label=nothing)
vline!(p2, [0.0], color=:black, label=nothing)

glob_layout = @layout [a b]
p = plot(p1, p2, legend=true, size=(800,500), dpi=300, layout=glob_layout, left_margin=(20,:px))

savedir = joinpath("tmp", "01-narahara_element_1_run")
savepath = plotsdir(savedir, "01-narahara_element_1_run_$(time_ns()).pdf")
savefig(p, savepath)