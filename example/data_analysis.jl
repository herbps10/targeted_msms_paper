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

data = DataFrame(CSV.File("./pnas_data.csv"))

function linear_working_model(beta::Vector, X::DataFrameRow)
    return first([1 X[:total_alive_bl]] * beta)
    #return beta[1]
end

function squared_error_loss(t, m)
    return (t - m)^2
end

Y = data[:, :curr_use_inj_adj]
A = data[:, :treatment]
X = data[:, [:total_alive_bl]]
    
result = TargetedMSM.treatment_effect_modification(
    Y,
    A,
    X,
    nothing,
    nothing,
    2, 
    squared_error_loss,
    linear_working_model; 
    Qbar = data[:, :Qbar],
    Qbar0 = data[:, :Qbar0],
    Qbar1 = data[:, :Qbar1],
    g = data[:, :g],
    linear = false,
    bayes = true,
    iterations = 50_000,
    seed = 1234,
    #proposal_sd = 0.09375,
    include_prior = true,
    prior = Normal(0, 1)
);

total_alive = [0:9;]
pred = result[:beta][1] .+ result[:beta][2] * total_alive

h = scatter(data[:, :total_alive_bl], result[:Psi_star])
h = scatter!(total_alive, pred)

density(result[:beta_post][:, 1])
plot!(Normal(result[:beta][1], result[:beta_se][1] / sqrt(size(data)[1])))

density(result[:beta_post][:, 2])
plot!(Normal(result[:beta][2], result[:beta_se][2] / sqrt(size(data)[1])))

h = scatter(result[:epsilon_post][:, 1], result[:beta_post][:, 1])
h = scatter(result[:epsilon_post][:, 2], result[:beta_post][:, 2])

CSV.write(
    "results/pnas_frequentist.csv",
    DataFrame(
        name  = ["beta1", "beta2"],
        n     = repeat([size(data)[1]], 2),
        beta  = result[:beta],
        lower = result[:beta_lower][:, 1],
        upper = result[:beta_upper][:, 1],
        se    = result[:beta_se][:, 1],
    )
)

CSV.write(
    "results/pnas_posterior.csv",
    DataFrame(
        beta1 = result[:beta_post][:, 1],
        beta2 = result[:beta_post][:, 2],
        epsilon1 = result[:epsilon_post][:, 1],
        epsilon2 = result[:epsilon_post][:, 2]
    )
)
