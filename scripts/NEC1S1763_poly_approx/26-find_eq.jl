using DrWatson
@quickactivate "masters1semester_project"

using Plots
using NonlinearSolve
# using Optim
# using DifferentialEquations

include(srcdir("narahara_papers.jl"))
using .narahara_papers

include(srcdir("NEC1S1763_poly_approx.jl"))
using .NEC1S1763_poly_approx

include(srcdir("tools.jl"))

#########################################################################################

target_u_data = NEC1S1763_poly_approx.convert_V_to_u(NEC1S1763_poly_approx.IV_curve_V_data)
target_v_data = NEC1S1763_poly_approx.convert_I_to_v(NEC1S1763_poly_approx.IV_curve_I_data)

u_start, u_end, u_N = minimum(target_u_data), maximum(target_u_data), 200
u_range = range(u_start, u_end, u_N)

u_start_root, u_end_root = u_start, u_end
u_span_root = (u_start_root, u_end_root)
u_eps_root = 1.1e-2

#########################################################################################

R = 5 # Ohm
C = 2e-9 # F
L = 5e-6 # H
V_B = 0.2 # V

R₀ = NEC1S1763_poly_approx.V_p/NEC1S1763_poly_approx.I_p
E = V_B/NEC1S1763_poly_approx.I_p/R
δ = R₀/R
ε = C*R₀*R/L

#########################################################################################

function find_roots(f, g, x_span, eps)
    x_segm_N = Int(round((x_span[2] - x_span[1])/eps))

    x_roots = []
    for i in 0:x_segm_N-1
        bound_l = x_span[1] + i*eps
        bound_r = x_span[1] + (i+1)*eps
        root = find_root(f, g, (bound_l, bound_r))
        if !isnothing(root)
            push!(x_roots, root)
        end
    end
    return unique!(x_roots)
end

function find_root(f, g, x_span)
    # delta(x, p) = f(x) - g(x)
    delta(x, p) = f(x) - g(x)
    if delta(x_span[1], 0) * delta(x_span[2], 0) > 0
        return nothing
    end
    prob = IntervalNonlinearProblem(delta, x_span)
    sol = solve(prob)
    return sol.u
end

#########################################################################################

a₀ = -0.0014491777067221945
a₁ = 0.4047903610420292
a₂ = -8.814345521066247
a₃ = 94.05841173048864
a₄ = -593.8071165028581
a₅ = 2362.6810193881056
a₆ = -5964.114570853057
a₇ = 9209.238555278202
a₈ = -7889.603353342613
a₉ = 2861.212170878973

as = [a₀,a₁,a₂,a₃,a₄,a₅,a₆,a₇,a₈,a₉]
bs = NEC1S1763_poly_approx.convert_dim_to_dimless_coeff(as)


f(x) = NEC1S1763_poly_approx.general_poly(x, bs)
g(x, E, δ) = E - δ*x
g(x) = g(x, E, δ)

u_roots = find_roots(f, g, u_span_root, u_eps_root)
v_at_roots = f.(u_roots)

#########################################################################################

p1 = plot(title="NEC1S1763 Eq",
    titlefontsize=12,
    xlabel="V, В", ylabel="I, A", 
    label=nothing,
    xlims=(minimum(target_u_data), maximum(target_u_data)),
    ylims=(minimum(target_v_data), maximum(target_v_data)),
    minorgrid=true
)

# Координатные оси
hline!(p1, [0.0], color=:black, label=nothing)
vline!(p1, [0.0], color=:black, label=nothing)

# Экспериментальные данные ВАХ
plot!(p1, target_u_data, target_v_data, ls=:solid, label="")

# Аппроксимация ВАХ
plot!(p1, u_range, f.(u_range), ls=:solid, label="")

# Нагрузочная прямая
plot!(p1, u_range, g.(u_range), ls=:solid, label="")

# Состоярия равновесия
scatter!(p1, u_roots, v_at_roots, label="")

p = plot(p1, legend=false, size=(500,500), dpi=300, left_margin=(20,:px))

savedir = joinpath("tmp", "26-find_eq")
savepath = plotsdir(savedir, "26-find_eq_$(time_ns()).pdf")
savefig(p, savepath)