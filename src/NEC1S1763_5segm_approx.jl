module NEC1S1763_5segm_approx

using DrWatson
@quickactivate "masters1semester_project"

#########################################################################################

struct obl
    obl_N::Int8

    u_st::Float64
    v_st::Float64

    u_min_bound::Float64
    u_max_bound::Float64
    is_u_in_bounds::Function

    k::Float64
    b::Float64
    δ::Float64
    E::Float64

    is_λ_real::Bool
    λ₁::Float64
    λ₂::Float64
    m::Float64
    n::Float64

    C₁::Function
    C₂::Function
    sol_u::Function
    sol_v::Function

    model::Function
    display_params::Function
end

obl(obl_N, u_min_bound, u_max_bound, k, b, δ, E, ε) = begin
    u_st = (E-b)/(k+δ)
    v_st = (k*E+b*δ)/(k+δ)

    is_u_in_bounds(u) = if u_min_bound<u<u_max_bound true else false end
    is_λ_real = if (k^2/ε^2 - (2k+4δ)/ε + 1 ≤ 0) false else true end

    α = (k+ε)/ε; β = (k+δ)/ε
    # λ₁, λ₂, m, n = NaN, NaN, NaN, NaN
    # C₁, C₂ = NaN, NaN
    # sol_u, sol_v = NaN, NaN
    if is_λ_real
        m = NaN; n = NaN
        λ₁ = -α/2 - sqrt(α^2/4 - β)
        λ₂ = -α/2 + sqrt(α^2/4 - β)
        C₂_real(u₀,v₀) = (v₀-v_st - (u₀-u_st)*(ε*λ₁+k))/ε/(λ₂-λ₁)
        C₁_real(u₀,v₀) = u₀-u_st - C₂_real(u₀, v₀)
        u_real(C₁,C₂,t) = C₁*exp(λ₁*t) + C₂*exp(λ₂*t) + u_st
        v_real(C₁,C₂,t) = C₁*(ε*λ₁+k)*exp(λ₁*t) + C₂*(ε*λ₂+k)*exp(λ₂*t) + v_st
        sol_u = u_real; sol_v = v_real
        C₁, C₂ = C₁_real, C₂_real
    else
        m = -α/2; n = sqrt(β-α^2/4)
        λ₁ = NaN; λ₂ = NaN
        C₁_imag(u₀,v₀) = u₀-u_st
        C₂_imag(u₀,v₀) = (v₀-v_st - C₁_imag(u₀,v₀)*(ε*m+k))/(ε*n)
        u_imag(C₁,C₂,t) = exp(m*t)*(C₁*cos(n*t) + C₂*sin(n*t)) + u_st
        v_imag(C₁,C₂,t) = exp(m*t)*(
            (ε*m*C₁ + ε*n*C₂ + k*C₁)*cos(n*t) + 
            (ε*m*C₂ - ε*n*C₁ + k*C₂)*sin(n*t)
        ) + v_st
        sol_u = u_imag; sol_v = v_imag
        C₁, C₂ = C₁_imag, C₂_imag
    end

    function model(dU, U, p, t)
        u,v = U
        k, b, δ, E, ε = p
        dU[1] = 1/ε * (v - (k*u + b))
        dU[2] = E - v - δ*u
        return nothing
    end

    obl(obl_N, u_st, v_st, 
        u_min_bound, u_max_bound, is_u_in_bounds,
        k, b, δ, E, 
        is_λ_real, λ₁, λ₂, m, n,
        C₁, C₂,
        sol_u, sol_v, 
        model,
        display_params
    )
end

function display_params(cur_obl::obl)
    labels = ["obl_N", "u_st", "v_st", "λ₁", "λ₂", "m", "n"]
    vals = [cur_obl.obl_N, cur_obl.u_st, cur_obl.v_st, cur_obl.λ₁, cur_obl.λ₂, cur_obl.m, cur_obl.n]
    println("Obl #", vals[1])
    for i in 2:length(labels)
        print("    $(labels[i]) = $(vals[i])")
    end
end

#########################################################################################

# general segment approximations
function general_segm_approx(x, xs, ys, k)
    i_match = NaN
    for i in 1:length(xs)-1
        if xs[i]≤x≤xs[i+1] i_match = i; break end
    end

    if isnan(i_match) return k*(x-xs[end]) + ys[end] end
    return (ys[i_match+1]-ys[i_match])/(xs[i_match+1]-xs[i_match]) * (x-xs[i_match]) + ys[i_match]
end

general_segm_approx(x, params) = general_segm_approx(x, params[1:3], params[4:6], params[end])

#########################################################################################

V_0, V_1, V_2, V_3, V_4 = 0.0, 0.07, 0.089, 0.3496, 0.42
I_0, I_1, I_2, I_3, I_4 = 0.0, 0.0056, 0.0033, 0.00223, 0.00073
k_dim = 0.0392671690168325

R₀ = NEC1S1763_5segm_approx.V_1/NEC1S1763_5segm_approx.I_1

u_0, u_1, u_2, u_3, u_4 = (V_0, V_1, V_2, V_3, V_4) ./ V_1
v_0, v_1, v_2, v_3, v_4 = (I_0, I_1, I_2, I_3, I_4) ./ I_1
k_dimless = k_dim * V_1/I_1

BAX_5segm_approx2(u) = general_segm_approx(u, 
    (u_0, u_1, u_2, u_3, u_4), (v_0, v_1, v_2, v_3, v_4), 
    k_dimless
)

# Нагрузочная прямая (load line)
function load_line_function(u, E, δ) return E-δ*u end

k₁ = 1.0;                       b₁ = 0.0
k₂ = R₀*(I_2-I_1)/(V_2-V_1);    b₂ = -(k₂ - 1.0)
k₃ = R₀*(I_3-I_2)/(V_3-V_2);    b₃ = -((I_3-I_2)/(V_3-V_2)*V_2 - I_2)/I_1
k₄ = R₀*(I_4-I_3)/(V_4-V_3);    b₄ = -((I_4-I_3)/(V_4-V_3)*V_3 - I_3)/I_1
k₅ = R₀*k_dimless;              b₅ = -(k_dimless*V_4 - I_4)/I_1

#########################################################################################

# TEST obl STRUCT
# obl(obl_N, u_min_bound, u_max_bound, k, b, δ, E, ε)
# obl1 = obl(1, 0.0, 1.0, 1.0, 0.0, -1.1, 0.0, 1.0)
# obl1.display_params()

# TEST general_4segm_approx
# general_4segm_approx(x) = general_4segm_approx(x, (0.0, 1.0, 2.0), (0.0, 1.0, 0.0), 1.0)
# println(general_4segm_approx(0.0))
# println(general_4segm_approx(0.5))
# println(general_4segm_approx(1.5))
# println(general_4segm_approx(3.5))

#########################################################################################

function create_obl_from_obl_N(obl_N, δ, E, ε)
    if obl_N == 1 return obl(1, -1e10, u_1, k₁, b₁, δ, E, ε) end
    if obl_N == 2 return obl(2, u_1, u_2, k₂, b₂, δ, E, ε) end
    if obl_N == 3 return obl(3, u_2, u_3, k₃, b₃, δ, E, ε) end
    if obl_N == 4 return obl(4, u_3, u_4, k₄, b₄, δ, E, ε) end
    if obl_N == 5 return obl(5, u_4, 1e10, k₅, b₅, δ, E, ε) end
    if obl_N == 6 return obl(6, u_3, u_4, k₄, b₄, δ, E, ε) end
    if obl_N == 7 return obl(7, u_2, u_3, k₃, b₃, δ, E, ε) end
    if obl_N == 8 return obl(8, u_1, u_2, k₂, b₂, δ, E, ε) end
end

#########################################################################################
    
end