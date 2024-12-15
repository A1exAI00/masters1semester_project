module NEC1S1763_poly_approx

using DrWatson
@quickactivate "masters1semester_project"

using Optim

#########################################################################################

IV_curve_V_data = [0.0, 0.01028524, 0.02062652, 0.02902881, 0.03937009, 0.04971137, 0.06005265, 0.07039393, 0.07879622, 0.0891375, 0.09947878, 0.12016134, 0.13955124, 0.1602338, 0.17897737, 0.19965993, 0.24942734, 0.29984108, 0.34960849, 0.40002223, 0.44978964, 0.49955705]
IV_curve_I_data = [0.0, 0.0016478873, 0.00299999994, 0.00396126752, 0.0047112675, 0.00519718298, 0.00548239424, 0.00557746466, 0.00353873232, 0.00329577458, 0.00320070416, 0.0030739436, 0.0030739436, 0.00292605628, 0.00288380276, 0.002862676, 0.00278873234, 0.00261971826, 0.0022288732, 0.00072887324, 0.00102464788, 0.00297887318]

V_p = 0.00557746466
I_p = 0.07039393

#########################################################################################

OPTIM_G_TOL = 0.0
OPTIM_ITERATIONS = 100000
OPTIM_ALG = LBFGS()
OPTIM_DISPLAY = false

#########################################################################################


function general_poly(x, as)
    res = 0.0
    for i in 0:length(as)-1
        res += as[i+1]*x^i
    end
    return res 
end

function general_poly_weighted_loss(params, x_data, y_data, weights)
    func_(x) = general_poly(x, params)
    sq_err = weights .* (y_data .- func_.(x_data)).^2
    return sum(sq_err)/length(sq_err)
end

function general_poly_optimize(params₀, x_data, y_data, weights)
    loss_(params) = general_poly_weighted_loss(params, x_data, y_data, weights)
    options = Optim.Options(g_tol=OPTIM_G_TOL, iterations=OPTIM_ITERATIONS)
    optim_result = optimize(loss_, params₀, OPTIM_ALG, options)
    if OPTIM_DISPLAY display(optim_result) end
    return Optim.minimizer(optim_result)
end

function easy_poly_optimize(degree, x_data, y_data, weights)
    params₀ = zeros(degree+1)
    return general_poly_optimize(params₀, x_data, y_data, weights)
end

#########################################################################################

function general_model(approx, dU, U, p, t)
    as, δ, E, ε = p; u, v = U
    dU[1] = 1/ε * (v - approx(u, as))
    dU[2] = E - v - δ*u
    return nothing
end

function general_integrate(model, U₀, t_span, as;
    abstol=nothing, reltol=nothing, alg=nothing, saveat=nothing)
    ABSTOL = 1e-5; abstol = if isnothing(abstol) ABSTOL end
    RELTOL = 1e-5; reltol = if isnothing(reltol) RELTOL end
    ALG = Tsit5(); alg = if isnothing(alg) ALG end

    prob = ODEProblem(model, U₀, t_span, (as, δ, E, ε))
    sol = if isnothing(saveat) 
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL) else
        solve(prob; alg=ALG, reltol=RELTOL, abstol=ABSTOL, saveat=saveat) end
    return sol
end

end