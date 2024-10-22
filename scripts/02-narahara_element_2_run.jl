using DrWatson
@quickactivate "masters_julia_project"

using Plots

include(srcdir("narahara_papers.jl"))
using .narahara_papers

#########################################################################################

V_offset = 0.08
I_offset = 0.001
N_V = 200

#########################################################################################

V_B = 0.6 # V

t_start, t_end = 0.0, 1e-5/4
t_span = (t_start, t_end)

u₀ = [0.0, 0.0]

#########################################################################################

sol = narahara_tunnel_diode.integrate_narahara_element_2(u₀, t_span, V_B)

sol_V = sol[1,:]
sol_I = sol[2,:]
sol_t = sol.t

max_V, min_V = maximum(sol_V), minimum(sol_V)
max_I, min_I = maximum(sol_I), minimum(sol_I)

V_range = range(min_V-V_offset, max_V+V_offset, N_V)
V_nullcline = narahara_tunnel_diode.V_nullcline_2.(V_range)
I_nullcline = narahara_tunnel_diode.I_nullcline_2.(V_range, V_B)

#########################################################################################

p1 = plot(title="Фазовое пространство одного элемента: V_B=$(V_B)",
    titlefontsize=12,
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
    title="Осциллограмма V(t)", 
    titlefontsize=12,
    xlabel="t, s", ylabel="V, В", 
    label=nothing,
    minorgrid=true
)
hline!(p2, [0.0], color=:black, label=nothing)
vline!(p2, [0.0], color=:black, label=nothing)

glob_layout = @layout [a{0.9w}; b]
p = plot(p1, p2, legend=true, size=(500,800), dpi=300, layout=glob_layout, left_margin=(20,:px))

savedir = joinpath("tmp", "02-narahara_element_2_run")
savepath = plotsdir(savedir, "02-narahara_element_2_run$(time_ns()).pdf")
savefig(p, savepath)