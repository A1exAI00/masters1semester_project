module narahara_papers

using DrWatson
@quickactivate "masters1semester_project"

using DifferentialEquations

#########################################################################################

# Element parameters
# paper: Trigger waves in an oscillator loop and their applications in voltage-controlled oscillation with high tuning gain
L_se = L_sh = 5e-6 # H
C = 2e-9 # F
R_se = R_sh = 5 # Ohm

# I-V curve data for NEC 1S1763
# paper: Self-injection Locking of Rotary Traveling Pulses in Resonant-Tunneling-Diode Transmission-Line Loop
IV_curve_V_data = [0.0, 0.01028524, 0.02062652, 0.02902881, 0.03937009, 0.04971137, 0.06005265, 0.07039393, 0.07879622, 0.0891375, 0.09947878, 0.12016134, 0.13955124, 0.1602338, 0.17897737, 0.19965993, 0.24942734, 0.29984108, 0.34960849, 0.40002223, 0.44978964, 0.49955705]
IV_curve_I_data = [0.0, 0.0016478873, 0.00299999994, 0.00396126752, 0.0047112675, 0.00519718298, 0.00548239424, 0.00557746466, 0.00353873232, 0.00329577458, 0.00320070416, 0.0030739436, 0.0030739436, 0.00292605628, 0.00288380276, 0.002862676, 0.00278873234, 0.00261971826, 0.0022288732, 0.00072887324, 0.00102464788, 0.00297887318]

#########################################################################################

# I-V curve approximation for NEC 1S1763
# paper: Trigger waves in an oscillator loop and their applications in voltage-controlled oscillation with high tuning gain
function I_D_1(x)
    a₁ = 0.005;  a₂ = 0.175
    a₃ = 1.4;    a₄ = 1.8535
    a₅ = 19.014; a₆ = 0.612
    a₇ = 54;     a₈ = 1.067
    a₉ = 0.0018; a₁₀ = 11.6
    a = (a₃+atan(a₄-a₅*x))
    d = log((1+exp(-a₆+a₇*x))/(1+exp(-a₆-a₇*x)))
    # I = a₁ * (a₂*(a*d-a₈) + a₉*(exp(a₁₀*x)-1))
    I = a₁ * (a₂*(a*d) + a₉*(exp(a₁₀*x)-1))
    return I
end

function V_nullcline_1(V)
    return I_D_1(V)
end

function I_nullcline_1(V, V_B)
    return 1/R_sh * (V_B-V)
end

#########################################################################################

# I-V curve approximation for RTD
# paper: Self-injection Locking of Rotary Traveling Pulses in Resonant-Tunneling-Diode Transmission-Line Loop
function I_D_2(x)
    S = 1 # μm^2
    jₚ = 1e5 # A/cm^2
    a₁ = 9.5e-10; a₂ = 5.792e-12
    a₃ = -0.59; a₄ = 13.0
    a₅ = 14.0; a₆ = 35.0; a₇ = 0.17
    a = 1 + exp(a₃+a₄*x)
    b = 1 + exp(a₃-a₄*x)
    c = a₁ * log(a/b)
    d = (π/2 + atan(a₅-a₆*x))
    e = a₂ * (exp(x/a₇) - 1)
    I = S * jₚ * (c*d + e)
    return I
end

function V_nullcline_2(V)
    return I_D_2(V)
end

function I_nullcline_2(V, V_B)
    return 1/R_sh * (V_B-V)
end

#########################################################################################

function narahara_element_1(du, u, p, t)
    L_sh, C, R_sh, V_B = p
    V, I = u
    du[1] = 1/C * (I - I_D_1(V))
    du[2] = 1/L_sh * (V_B - I*R_sh - V)
    return nothing
end

function narahara_element_2(du, u, p, t)
    L_sh, C, R_sh, V_B = p
    V, I = u
    du[1] = 1/C * (I - I_D_2(V))
    du[2] = 1/L_sh * (V_B - I*R_sh - V)
    return nothing
end

#########################################################################################

function integrate_narahara_element_1(u₀, t_span, V_B, 
    abstol=nothing, reltol=nothing, alg=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    ALG = Tsit5()

    if isnothing(reltol)
        reltol = RELTOL
    end
    if isnothing(abstol)
        abstol = ABSTOL
    end
    if isnothing(alg)
        alg = ALG
    end

    prob = ODEProblem(narahara_element_1, u₀, t_span, (L_sh, C, R_sh, V_B))
    sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)
    return sol
end

function integrate_narahara_element_2(u₀, t_span, V_B,
    abstol=nothing, reltol=nothing, alg=nothing)
    ABSTOL = 1e-5
    RELTOL = 1e-5
    ALG = Tsit5()

    if isnothing(reltol)
        reltol = RELTOL
    end
    if isnothing(abstol)
        abstol = ABSTOL
    end
    if isnothing(alg)
        alg = ALG
    end

    prob = ODEProblem(narahara_element_2, u₀, t_span, (L_sh, C, R_sh, V_B))
    sol = solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL)
    return sol
end


end