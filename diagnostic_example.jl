using TargetedMSM
using LinearAlgebra
using Distributions
using Random
using DataFrames
using DataFramesMeta
using ProgressMeter
using Optim
using CSV
using Plots
using StatsPlots
using StatsFuns
using GLM

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

rng = MersenneTwister(1253)
data = simulate_data(200, rng, linear = false)

results1 = TargetedMSM.treatment_effect_modification(
    data[:, :Y], 
    data[:, :A], 
    data[:, [:X1, :X2, :X3, :X4]], 
    @formula(A ~ X1 + X2 + X3 + X4),
    @formula(Y ~ X1 + X2 + X3 + X4 + A + X1*A + X2*A + X3*A + X4*A),
    2, 
    squared_error_loss,
    linear_working_model; 
    linear = false,
    bayes = true,
    iterations = 10_000,
    seed = 1234,
    proposal_sd = nothing,
    include_prior = true,
    prior = Normal(0, 1),
    dL = dL,
    dL_dβ = dL_dβ
)

density(results1[:beta_post][:, 2])
scatter(results1[:epsilon_post][:, 2], results1[:beta_post][:, 2])

# 1254
rng = MersenneTwister(1254)
data = simulate_data(20, rng, linear = false)

results2 = TargetedMSM.treatment_effect_modification(
    data[:, :Y], 
    data[:, :A], 
    data[:, [:X1, :X2, :X3, :X4]], 
    @formula(A ~ X1 + X2 + X3 + X4),
    @formula(Y ~ X1 + X2 + X3 + X4 + A + X1*A + X2*A + X3*A + X4*A),
    2, 
    squared_error_loss,
    linear_working_model; 
    linear = false,
    bayes = true,
    iterations = 10_000,
    seed = 1234,
    proposal_sd = nothing,
    include_prior = true,
    prior = Normal(0, 1),
    dL = dL,
    dL_dβ = dL_dβ
)

density(results2[:beta_post][:, 2])
scatter(results2[:epsilon_post][:, 2], results2[:beta_post][:, 2])

CSV.write(
    "results/diagnostic.csv",
    DataFrame(
        example1_beta_post = results1[:beta_post][:, 2],
        example1_epsilon_post = results1[:epsilon_post][:, 2],

        example2_beta_post = results2[:beta_post][:, 2],
        example2_epsilon_post = results2[:epsilon_post][:, 2]
    )
)
