using StatsFuns, Random, Distributions

function simulate_data(N, rng; linear = true)
    data = DataFrame(
        p1 = fill(0., N),
        p2 = fill(0., N),
        X1 = fill(0., N),
        X2 = fill(0., N),
        X3  = fill(0., N),
        X4  = fill(0., N),
        A  = fill(0, N),
        Y0 = fill(0.0, N),
        Qbar0 = fill(0.0, N),
        Y1 = fill(0.0, N),
        Qbar1 = fill(0.0, N),
        Y  = fill(0.0, N),
        diff = fill(0.0, N)
    )

    data[:, :X1] = rand(rng, Normal(0, 1), N)
    data[:, :X2] = rand(rng, Normal(0, 1), N)
    data[:, :X3] = rand(rng, Normal(0, 1), N)
    data[:, :X4] = rand(rng, Normal(0, 1), N)

    for i = 1:N
        p1 = 0.5 * data[i, :X1] - 0.5 * data[i, :X2] + 0.2 * data[i, :X3] - 0.1 * data[i, :X4]
        p2 = data[i, :X2] + data[i, :X3]

        data[i, :p1] = p1
        data[i, :p2] = p2

        data[i, :A]  = rand(rng, Binomial(1, logistic(p1)), 1)[1]


        if linear == true
            data[i, :Qbar0] = p2
            data[i, :Qbar1] = p2 + 3.0 + 1.5 * data[i, :X4]

            data[i, :Y0] = rand(rng, Normal(data[i, :Qbar0], 0.1))
            data[i, :Y1] = rand(rng, Normal(data[i, :Qbar1], 0.1))

            data[i, :Y]  = data[i, :A] == 1 ? data[i, :Y1] : data[i, :Y0]
        else
            data[i, :Qbar0] = p2
            data[i, :Qbar1] = p2 + 2 + data[i, :X4]

            data[i, :Y0] = rand(rng, Bernoulli(logistic(data[i, :Qbar0])))
            data[i, :Y1] = rand(rng, Bernoulli(logistic(data[i, :Qbar1])))

            data[i, :Y]  = data[i, :A] == 1 ? data[i, :Y1] : data[i, :Y0]

        end

        data[i, :diff] = data[i, :Y1] - data[i, :Y0]
    end

    return data
end
