module NEC1S1763_3segm_approx
    
using DrWatson
@quickactivate "masters1semester_project"

#########################################################################################

module NEC1S1763_3segm_params
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

    function display_params()
        println("General parameters")
        println("   R = ", R, " Ohm")
        println("   C = ", C, " F")
        println("   L = ", L, " H")
        println("   V_B = ", V_B, " V")
        println("")
        println("   V_p = ", V_p)
        println("   I_p = ", I_p)
        println("   V_v = ", V_v)
        println("   I_v = ", I_v)
        println("")
        println("   k₁ = ", k₁)
        println("   k₂ = ", k₂)
        println("   k₃ = ", k₃)
        println("")
    end

    export R, C, L, V_B, V_p, I_p, V_v, I_v, k₁, k₂, k₃, R₀, E, δ, ε
end

#########################################################################################

module obl1
    using ..NEC1S1763_3segm_params
    u_min = 0.0
    u_max = 1.0

    u₀_st = v₀_st = E/(1+δ)
    α = (1+ε)/ε
    β = (1+δ)/ε
    λ₁ = -α/2 - sqrt(α^2/4 - β)
    λ₂ = -α/2 + sqrt(α^2/4 - β)
    C₂(v₀) = (v₀+v₀_st*ε*λ₁-(ε*λ₁+1))/ε/(λ₂-λ₁)
    C₁(v₀) = 1 - u₀_st - C₂(v₀)

    u(C₁, C₂, t) = C₁         *exp(λ₁*t) + C₂         *exp(λ₂*t) + u₀_st
    v(C₁, C₂, t) = C₁*(ε*λ₁+1)*exp(λ₁*t) + C₂*(ε*λ₂+1)*exp(λ₂*t) + v₀_st
    u(t, v₀) = u(C₁(v₀), C₂(v₀), t)
    v(t, v₀) = v(C₁(v₀), C₂(v₀), t)

    function display_params()
        println("I obl:")
        println("   α = ", α)
        println("   β = ", β)
        println("   λ₁ = ", λ₁)
        println("   λ₂ = ", λ₂)
        println("")
    end

    function model(dU, U, p, t)
        u,v = U
        ε,E,δ = p
        dU[1] = 1/ε * (v - u)
        dU[2] = E - v - δ*u
        return nothing
    end

    function condition(U, t, int)
        u,v = U
        if v ≤ u
            return 1.0
        else
            return u_max - u
        end
    end
end

module obl2
    using ..NEC1S1763_3segm_params
    u_min = 1.0
    u_max = V_v/V_p

    u₀_st = (E - (1-R₀*k₂))/(R₀*k₂+δ)
    v₀_st = (R₀*k₂*E + δ*(1-R₀*k₂))/(R₀*k₂+δ)
    α = (R₀*k₂+ε)/ε
    β = (R₀*k₂+δ)/ε
    m = -α/2
    n = sqrt(β-α^2/4)
    C₁(v₀) = 1 - u₀_st
    C₂(v₀) = (v₀ - v₀_st - C₁(v₀)*(ε*m+R₀*k₂))/(ε*n)

    u(C₁, C₂, t) = exp(m*t)*(C₁*cos(n*t) + C₂*sin(n*t)) + u₀_st
    v(C₁, C₂, t) = exp(m*t)*(
        (ε*m*C₁ + ε*n*C₂ + R₀*k₂*C₁)*cos(n*t) + 
        (ε*m*C₂ - ε*n*C₁ + R₀*k₂*C₂)*sin(n*t)
    ) + v₀_st
    u(t, v₀) = u(C₁(v₀), C₂(v₀), t)
    v(t, v₀) = v(C₁(v₀), C₂(v₀), t)

    function display_params()
        println("II obl:")
        println("   α = ", α)
        println("   β = ", β)
        println("   m = ", m)
        println("   n = ", n)
        println("")
    end

    function model(dU, U, p, t)
        u,v = U
        ε,E,δ,R₀,k₂ = p
        dU[1] = 1/ε * (v - (R₀*k₂*u + (1-R₀*k₂)))
        dU[2] = E - v - δ*u
        return nothing
    end

    function condition(U, t, int)
        u,v = U
        return u_max - u 
    end
end

module obl3
    using ..NEC1S1763_3segm_params
    u_min = V_v/V_p
    u_max = Inf

    u₀_st = (E - (I_v-k₃*V_v)/I_p)/(R₀*k₃+δ)
    v₀_st = (R₀*k₃*E + δ*(I_v-k₃*V_v)/I_p)/(R₀*k₃+δ)
    α = (R₀*k₃+ε)/ε
    β = (R₀*k₃+δ)/ε
    m = -α/2
    n = sqrt(β-α^2/4)
    C₁(v₀) = V_v/V_p - u₀_st
    C₂(v₀) = (v₀ - v₀_st - C₁(v₀)*(ε*m+R₀*k₃))/(ε*n)

    u(C₁, C₂, t) = exp(m*t)*(C₁*cos(n*t) + C₂*sin(n*t)) + u₀_st
    v(C₁, C₂, t) = exp(m*t)*(
        (ε*m*C₁ + ε*n*C₂ + R₀*k₃*C₁)*cos(n*t) + 
        (ε*m*C₂ - ε*n*C₁ + R₀*k₃*C₂)*sin(n*t)
    ) + v₀_st
    u(t, v₀) = u(C₁(v₀), C₂(v₀), t)
    v(t, v₀) = v(C₁(v₀), C₂(v₀), t)

    function display_params()
        println("III obl:")
        println("   α = ", α)
        println("   β = ", β)
        println("   m = ", m)
        println("   n = ", n)
        println("")
    end

    function model(dU, U, p, t)
        u,v = U
        ε,E,δ,R₀,k₃ = p
        dU[1] = 1/ε * (v - (R₀*k₃*u + (I_v-k₃*V_v)/I_p))
        dU[2] = E - v - δ*u
        return nothing
    end

    function condition(U, t, int)
        u,v = U
        ε,E,δ,R₀,k₃ = int.p
        if v ≥ (R₀*k₃*u + (I_v-k₃*V_v)/I_p)
            return 1.0
        else
            return u - u_min
        end
    end
end

module obl4
    using ..NEC1S1763_3segm_params
    using ..obl2
    u_min = 1.0
    u_max = V_v/V_p

    u₀_st = obl2.u₀_st
    v₀_st = obl2.v₀_st
    α = obl2.α
    β = obl2.β
    m = obl2.m
    n = obl2.n
    C₁(v₀) = V_v/V_p - obl2.u₀_st
    C₂(v₀) = (v₀ - obl2.v₀_st - C₁(v₀)*(ε*m+R₀*k₂))/(ε*n)

    u(C₁, C₂, t) = obl2.u(C₁, C₂, t)
    v(C₁, C₂, t) = obl2.v(C₁, C₂, t)
    u(t, v₀) = u(C₁(v₀), C₂(v₀), t)
    v(t, v₀) = v(C₁(v₀), C₂(v₀), t)

    function display_params()
        println("IIII obl:")
        println("   α = ", α)
        println("   β = ", β)
        println("   m = ", m)
        println("   n = ", n)
        println("")
    end

    model = obl2.model

    function condition(U, t, int)
        u,v = U
        return u_min - u 
    end
end

using .NEC1S1763_3segm_params
using .obl1
using .obl2
using .obl3
using .obl4

NEC1S1763_3segm_params.display_params()
obl1.display_params()
obl2.display_params()
obl3.display_params()
obl4.display_params()

#########################################################################################

# general 3 segment approximations
function general_3segm_approx(x, params)
    xa, xb, ya, yb, k = params
    y = NaN
    if 0.0 ≤ x ≤ xa
        y = ya/xa * x
    elseif xa ≤ x ≤ xb
        y = (yb-ya)/(xb-xa) * (x-xa) + ya
    elseif x ≥ xb
        y = k*(x-xb) + yb
    end
    return y
end

# 3 segment I-V curve approximations
BAX_3segm_approx1(V) = general_3segm_approx(V, (V_p, V_v, I_p, I_v, k₃))
BAX_3segm_approx2(u) = general_3segm_approx(u, (1.0, V_v/V_p, 1.0, I_v/I_p, k₃*R₀))

# Нагрузочная прямая (load line)
function load_line_function(u) return E-δ*u end

#########################################################################################

using DifferentialEquations

function all_obl_affect!(int)
    terminate!(int)
end

function all_obl_integrate(U₀, t_span, L_N;
    w_callback=true, abstol=nothing, reltol=nothing, alg=nothing, saveat=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    ALG = Tsit5()

    reltol = if isnothing(reltol) RELTOL end
    abstol = if isnothing(abstol) ABSTOL end
    alg = if isnothing(alg) ALG end

    if L_N == 1
        cb = if w_callback ContinuousCallback(obl1.condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl1.model, U₀, t_span, (ε,E,δ))
    elseif L_N == 2
        cb = if w_callback ContinuousCallback(obl2.condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl2.model, U₀, t_span, (ε,E,δ,R₀,k₂))
    elseif L_N == 3
        cb = if w_callback ContinuousCallback(obl3.condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl3.model, U₀, t_span, (ε,E,δ,R₀,k₃))
    elseif L_N == 4
        cb = if w_callback ContinuousCallback(obl4.condition, all_obl_affect!) else nothing end
        prob = ODEProblem(obl4.model, U₀, t_span, (ε,E,δ,R₀,k₂))
    end
    sol = if isnothing(saveat) 
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, callback=cb) else
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, callback=cb, saveat=saveat) end
    return sol
end

end