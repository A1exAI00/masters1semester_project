# using DrWatson
# @quickactivate "masters1semester_project"

# using Plots
# using Optim

# include(srcdir("narahara_papers.jl"))
# using .narahara_papers

# include(srcdir("NEC1S1763_5segm_approx.jl"))
# using .NEC1S1763_5segm_approx

# #########################################################################################

# R = 5 # Ohm
# C = 2e-9 # F
# L = 5e-6 # H
# V_B = 0.2 # V

# R₀ = NEC1S1763_5segm_approx.V_1/NEC1S1763_5segm_approx.I_1
# E = NEC1S1763_5segm_approx.V_B/NEC1S1763_5segm_approx.I_1/R
# δ = R₀/R
# ε = C*R₀*R/L

# #########################################################################################

# u_start, u_end, u_N = 0.0, 5.0, 200
# u_range = range(u_start, u_end, u_N)

# v_curve_1 = NEC1S1763_5segm_approx.BAX_5segm_approx2.(u_range)
# v_curve_2 = NEC1S1763_5segm_approx.load_line_function.(u_range, E, δ)

# #########################################################################################

# p1 = plot(title="Diode I-V curve, NEC 1S1763",
#     titlefontsize=12,
#     xlabel="u", ylabel="v", 
#     label=nothing,
#     xlims=(V_min, V_max),
#     ylims=(minimum(I_V_curve), maximum(I_V_curve)),
#     minorgrid=true
# )
# plot!(p1, V_range, I_V_curve, ls=:solid, label="")
# hline!(p1, [0.0], color=:black, label=nothing)
# vline!(p1, [0.0], color=:black, label=nothing)

# p = plot(p1, legend=false, size=(500,500), dpi=300, left_margin=(20,:px))

# savedir = joinpath("tmp", "00-plot_i-v_curve")
# savepath = plotsdir(savedir, "00-plot_i-v_curve_$(time_ns()).pdf")
# savefig(p, savepath)