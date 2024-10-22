# Circut elements parameters
R = 5 # Ohm
C = 2e-9 # F
L = 5e-6 # H
V_B = 0.2

# I-V curve parameters
V_p = 0.1 # V
I_p = 0.0058 # A
V_v = 0.4 # V
I_v = 0.001 # A

k₁ = I_p/V_p
k₂ = (I_v-I_p)/(V_v-V_p)
k₃ = 0.025

# Dimless parameters
R₀ = V_p/I_p
E = V_B/I_p/R
δ = R₀/R
ε = C*R₀*R/L

# lambda parameters
α₁ = (1+ε)/ε
β₁ = (1+δ)/ε
# m₁ = -α₁/2
# n₁ = sqrt(β₁-α₁^2/4)

# α₂ = (R₀*k₂+ε)/ε
# β₂ = ((R₀*k₂)^2+δ*ε)/(ε^2)
# m₂ = -α₂/2
# n₂ = sqrt(β₂-α₂^2/4)

# α₃ = (R₀*k₃+ε)/ε
# β₃ = ((R₀*k₃)^2+δ*ε)/(ε^2)
# m₃ = -α₃/2
# n₃ = sqrt(β₃-α₃^2/4)

#########################################################################################

println("R = ", R, " Ohm")
println("C = ", C, " F")
println("L = ", L, " H")
println("V_B = ", V_B, " V")
println("")
println("V_p = ", V_p)
println("I_p = ", I_p)
println("V_v = ", V_v)
println("I_v = ", I_v)
println("")
println("k₁ = ", k₁)
println("k₂ = ", k₂)
println("k₃ = ", k₃)
println("")
println("R₀ = ", R₀)
println("E = ", E)
println("δ = ", δ)
println("ε = ", ε)
println("")
println("α₁ = ", α₁)
println("β₁ = ", β₁)
# println("m₁ = ", m₁)
# println("n₁ = ", n₁)
# println("")
# println("α₂ = ", α₂)
# println("β₂ = ", β₂)
# println("m₂ = ", m₂)
# println("n₂ = ", n₂)
# println("")
# println("α₃ = ", α₃)
# println("β₃ = ", β₃)
# println("m₃ = ", m₃)
# println("n₃ = ", n₃)
# println("")