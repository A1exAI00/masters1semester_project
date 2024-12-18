using DrWatson
@quickactivate "masters1semester_project"

# using Plots
# using NonlinearSolve
# using Optim
# using DifferentialEquations

include(srcdir("NEC1S1763_poly_approx.jl"))
using .NEC1S1763_poly_approx

include(srcdir("tools.jl"))

#########################################################################################

I_p = 0.0058
V_p = 0.07039393

a₀ = -0.0003295956895430551
a₁ = 0.2511481502100568
a₂ = -3.8863598765249816
a₃ = 25.870115639283853
a₄ = -87.90416287262619
a₅ = 157.12901171913
a₆ = -138.74257222597265
a₇ = 46.96385554336421

as = [a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇]

#########################################################################################

# u = V/Vp
# v = I/Ip
# I_D(V) = ∑(a_i * x^i)
# f(u) = 1/I_p * I_D(V_p*u) = 1/I_p * ∑(a_i * V_p^i * u^i) = ∑(a_i*V_p^i/I_p * u^i) = ∑(b_i * u^i)
# b_i = a_i*V_p^i/I_p

# for i in eachindex(as)
#     a = as[i]
#     b = a * V_p^(i-1) / I_p
#     println("b$(sub_string(string(i-1))) = $(b)")
# end

bs = NEC1S1763_poly_approx.convert_dim_to_dimless_coeff(as)

println(bs)