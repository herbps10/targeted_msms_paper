using Distributed

println("Adding ", Distributed.nworkers(), " processes")
addprocs(Distributed.nworkers())

using CSV
using SharedArrays

@everywhere begin
using TargetedMSM
using LinearAlgebra
using Distributions
using Random
using DataFrames
using DataFramesMeta
#using ProgressMeter
using Optim
#using Plots
using GLM
#using LaTeXStrings

include("simulate_data.jl")

function linear_working_model(beta::Vector, X::DataFrameRow)
    return first([1 X[:X4]] * beta)
end

function squared_error_loss(t, m)
    return (t - m)^2
end

function dL(t, beta, X)
    -2 * (t - linear_working_model(beta, X)) .* [1 X[:X4]]'
end

function dL_dβ(t, beta, X)
   2 * [1 X[:X4]]' * [1 X[:X4]] 
end
end

rng = MersenneTwister(1234)

large_data = simulate_data(100000, rng, linear = false)
optim = optimize(beta -> sum((large_data[:, :diff] .- beta[1] - beta[2] .* large_data[:, :X4]).^2), repeat([0.0], 2), LBFGS(), autodiff = :forward)
beta0 = Optim.minimizer(optim)

# TMLE simulation study
g_correct = [true, false]
Q_correct = [true, false]
Ns = [50, 100, 250, 500, 750, 1000];
num_simulations = 500;

simulations = DataFrame(
    collect(Iterators.product(g_correct, Q_correct, Ns, 1:num_simulations))
);

rename!(simulations, ["g_correct", "Q_correct", "N", "index"])

bayes = true

println("Starting simulations...")
@sync @distributed for index = 1:size(simulations, 1)
    println(index)
    rng = MersenneTwister(index)
    data = simulate_data(simulations[index, :N], rng, linear = false)

    g_formula = @formula(A ~ X1 + X2 + X3 + X4)
    Q_formula = @formula(Y ~ X1 + X2 + X3 + X4 + A + X1*A + X2*A + X3*A + X4*A)

    g_formula_misspec = @formula(A ~ X1 + X4)
    Q_formula_misspec = @formula(Y ~ X1 + X4 + A + X1*A + X4*A)

    result = TargetedMSM.treatment_effect_modification(
        data[:, :Y], 
        data[:, :A], 
        data[:, [:X1, :X2, :X3, :X4]], 
        (simulations[index, :g_correct] == true ? g_formula : g_formula_misspec),
        (simulations[index, :Q_correct] == true ? Q_formula : Q_formula_misspec),
        2, 
        squared_error_loss,
        linear_working_model;
        var_formula = @formula(Y2 ~ X1 + X2 + X3 + X4 + A + X1*A + X2*A + X3*A + X4*A),
        seed = 1234,
        proposal_sd = nothing,
				linear = false,
        iterations = 10_000,
        bayes = bayes,
        include_prior = true,
        prior = Normal(0, 1),
        dL = dL,
        dL_dβ = dL_dβ
    )

    output = string("worker-results/results-", Distributed.myid(), ".csv")

    CSV.write(output, DataFrame(
       seed = index,
       N = simulations[index, :N],
       g_correct = simulations[index, :g_correct],
       Q_correct = simulations[index, :Q_correct],
       beta1 = result[:beta][1], 
       beta2 = result[:beta][2], 
       beta1_lower = result[:beta_lower][1],
       beta1_upper = result[:beta_upper][1],
       beta2_lower = result[:beta_lower][2],
       beta2_upper = result[:beta_upper][2],

       bayes_beta1 = bayes == true ? result[:beta_bayes][1] : 0, 
       bayes_beta2 = bayes == true ? result[:beta_bayes][2] : 0, 
       bayes_beta1_lower = bayes == true ? result[:beta_lower_bayes][1] : 0,
       bayes_beta1_upper = bayes == true ? result[:beta_upper_bayes][1] : 0,
       bayes_beta2_lower = bayes == true ? result[:beta_lower_bayes][2] : 0,
       bayes_beta2_upper = bayes == true ? result[:beta_upper_bayes][2] : 0
    ), append = true)

    GC.gc()
end
