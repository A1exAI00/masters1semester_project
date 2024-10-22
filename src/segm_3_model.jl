module segm_3_model
    
using DrWatson
@quickactivate "masters_julia_project"

using DifferentialEquations

#########################################################################################

# Circut elements parameters
R = 5 # Ohm
C = 2e-9 # F
L = 5e-6 # H
V_B = 0.2 # V

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

#########################################################################################

# mutable struct obl
#     u_min
#     u_max
#     v_min
#     v_max
#     u₀
#     v₀
#     α
#     β
#     ss_type
#     λ₁
#     λ₂
#     m
#     n
# end


#########################################################################################

# L₁ parameters
u₀_st_1 = v₀_st_1 = E/(1+δ)
α₁ = (1+ε)/ε
β₁ = (1+δ)/ε
λ₁_I = -α₁/2 - sqrt(α₁^2/4 - β₁)
λ₂_I = -α₁/2 + sqrt(α₁^2/4 - β₁)
C₂_1(v₀) = (v₀+v₀_st_1*ε*λ₁_I-(ε*λ₁_I+1))/ε/(λ₂_I-λ₁_I)
C₁_1(v₀) = 1 - u₀_st_1 - C₂_1(v₀)

u₁(C₁, C₂, t) = C₁*exp(λ₁_I*t) + C₂*exp(λ₂_I*t) + u₀_st_1
v₁(C₁, C₂, t) = C₁*(ε*λ₁_I+1)*exp(λ₁_I*t) + C₂*(ε*λ₂_I+1)*exp(λ₂_I*t) + v₀_st_1
u₁(t, v₀) = u₁(C₁_1(v₀), C₂_1(v₀), t)
v₁(t, v₀) = v₁(C₁_1(v₀), C₂_1(v₀), t)

#########################################################################################

# L₂ parameters
u₀_st_2 = (E - (1-R₀*k₂))/(R₀*k₂+δ)
v₀_st_2 = (R₀*k₂*E + δ*(1-R₀*k₂))/(R₀*k₂+δ)
α₂ = (R₀*k₂+ε)/ε
β₂ = (R₀*k₂+δ)/ε
m₂ = -α₂/2
n₂ = sqrt(β₂-α₂^2/4)
C₁_2(v₀) = 1 - u₀_st_2
C₂_2(v₀) = (v₀ - v₀_st_2 - C₁_2(v₀)*(ε*m₂+R₀*k₂))/(ε*n₂)

u₂(C₁, C₂, t) = exp(m₂*t)*(C₁*cos(n₂*t) + C₂*sin(n₂*t)) + u₀_st_2
v₂(C₁, C₂, t) = exp(m₂*t)*(
    (ε*m₂*C₁ + ε*n₂*C₂ + R₀*k₂*C₁)*cos(n₂*t) + 
    (ε*m₂*C₂ - ε*n₂*C₁ + R₀*k₂*C₂)*sin(n₂*t)
) + v₀_st_2
u₂(t, v₀) = u₂(C₁_2(v₀), C₂_2(v₀), t)
v₂(t, v₀) = v₂(C₁_2(v₀), C₂_2(v₀), t)

#########################################################################################

# L₃ parameters
u₀_st_3 = (E - (I_v-k₃*V_v)/I_p)/(R₀*k₃+δ)
v₀_st_3 = (R₀*k₃*E + δ*(I_v-k₃*V_v)/I_p)/(R₀*k₃+δ)
α₃ = (R₀*k₃+ε)/ε
β₃ = (R₀*k₃+δ)/ε
m₃ = -α₃/2
n₃ = sqrt(β₃-α₃^2/4)
C₁_3(v₀) = V_v/V_p - u₀_st_3
C₂_3(v₀) = (v₀ - v₀_st_3 - C₁_3(v₀)*(ε*m₃+R₀*k₃))/(ε*n₃)

u₃(C₁, C₂, t) = exp(m₃*t)*(C₁*cos(n₃*t) + C₂*sin(n₃*t)) + u₀_st_3
v₃(C₁, C₂, t) = exp(m₃*t)*(
    (ε*m₃*C₁ + ε*n₃*C₂ + R₀*k₃*C₁)*cos(n₃*t) + 
    (ε*m₃*C₂ - ε*n₃*C₁ + R₀*k₃*C₂)*sin(n₃*t)
) + v₀_st_3
u₃(t, v₀) = u₃(C₁_3(v₀), C₂_3(v₀), t)
v₃(t, v₀) = v₃(C₁_3(v₀), C₂_3(v₀), t)

#########################################################################################

# L₄ parameters
m₄ = m₂
n₄ = n₂
C₁_4(v₀) = V_v/V_p - u₀_st_2
C₂_4(v₀) = (v₀ - v₀_st_2 - C₁_4(v₀)*(ε*m₄+R₀*k₂))/(ε*n₄)

u₄(C₁, C₂, t) = u₂(C₁, C₂, t) 
v₄(C₁, C₂, t) = v₂(C₁, C₂, t) 
u₄(t, v₀) = u₄(C₁_4(v₀), C₂_4(v₀), t)
v₄(t, v₀) = v₄(C₁_4(v₀), C₂_4(v₀), t)

#########################################################################################

function print_parameters()
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
    println("λ₁_I = ", λ₁_I)
    println("λ₂_I = ", λ₂_I)
    println("")
    println("α₂ = ", α₂)
    println("β₂ = ", β₂)
    println("m₂ = ", m₂)
    println("n₂ = ", n₂)
    println("")
    println("α₃ = ", α₃)
    println("β₃ = ", β₃)
    println("m₃ = ", m₃)
    println("n₃ = ", n₃)
end

print_parameters()

#########################################################################################

function f(u)
    if 0 ≤ u ≤ 1 
        return u
    elseif 1 ≤ u ≤ V_v/V_p
        return R₀*k₂*(u-1) + 1
    else
        return R₀*k₃*u - (k₃*V_v - I_v)/I_p
    end
end

function g(u) return E-δ*u end

#########################################################################################

function obl1_model(dU, U, p, t)
    u,v = U
    ε,E,δ = p
    dU[1] = 1/ε * (v - u)
    dU[2] = E - v - δ*u
    return nothing
end

function obl2_model(dU, U, p, t)
    u,v = U
    ε,E,δ,R₀,k₂ = p
    dU[1] = 1/ε * (v - (R₀*k₂*u + (1-R₀*k₂)))
    dU[2] = E - v - δ*u
    return nothing
end

function obl3_model(dU, U, p, t)
    u,v = U
    ε,E,δ,R₀,k₃ = p
    dU[1] = 1/ε * (v - (R₀*k₃*u + (I_v-k₃*V_v)/I_p))
    dU[2] = E - v - δ*u
    return nothing
end

#########################################################################################

function obl1_condition(U, t, int)
    u,v = U
    if v ≤ u
        return 1.0
    else
        return 1.0 - u
    end
end

function obl2_condition(U, t, int)
    u,v = U
    return V_v/V_p - u 
end

function obl3_condition(U, t, int)
    u,v = U
    ε,E,δ,R₀,k₃ = int.p
    if v ≥ (R₀*k₃*u + (I_v-k₃*V_v)/I_p)
        return 1.0
    else
        return u - V_v/V_p
    end
end

function obl4_condition(U, t, int)
    u,v = U
    return 1.0 - u 
end

#########################################################################################

function all_obl_affect!(int)
    terminate!(int)
end

#########################################################################################

function all_obl_integrate(U₀, t_span, L_N;
    w_callback=true, abstol=nothing, reltol=nothing, alg=nothing, saveat=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    ALG = Tsit5()

    reltol = if isnothing(reltol) RELTOL end
    abstol = if isnothing(abstol) ABSTOL end
    alg = if isnothing(alg) ALG end

    if L_N == 1
        cb = if w_callback ContinuousCallback(obl1_condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl1_model, U₀, t_span, (ε,E,δ))
    elseif L_N == 2
        cb = if w_callback ContinuousCallback(obl2_condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl2_model, U₀, t_span, (ε,E,δ,R₀,k₂))
    elseif L_N == 3
        cb = if w_callback ContinuousCallback(obl3_condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl3_model, U₀, t_span, (ε,E,δ,R₀,k₃))
    elseif L_N == 4
        cb = if w_callback ContinuousCallback(obl4_condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl2_model, U₀, t_span, (ε,E,δ,R₀,k₂))
    end
    sol = if isnothing(saveat) 
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, callback=cb) else
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, callback=cb, saveat=saveat) end
    return sol
end

end
