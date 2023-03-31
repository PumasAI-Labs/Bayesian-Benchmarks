using Pumas
using DataFrames
using CSV

depot_1cmt_prop = @model begin

    @param begin
        TVCL ~ LogNormal(log(0.4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVKA ~ LogNormal(log(1), 1)
        #σ_p ~ Constrained(Normal(0, 0.5), lower = 0, upper = Inf)
        σ_p ~ truncated(Normal(0, 0.5), 0, Inf)
        C ~ LKJCholesky(3, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(3), Diagonal([0.4, 0.4, 0.4].^2)),
            lower = zeros(3),
            upper = fill(Inf, 3),
            init = ones(3)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(3)) # Z in the Stan code
    end

    @covariates lloq

    @pre begin

        _lloq = lloq

        # compute the η from the ηstd
        # using lower Cholesky triangular matrix
        η = ω .* (getchol(C).L * ηstd)

        # PK parameters
        CL = TVCL * exp(η[1])
        Vc = TVVC * exp(η[2])
        Ka = TVKA * exp(η[3])
    end

    @dynamics begin
        Depot' = -Ka * Depot
        Central' = Ka * Depot - (CL / Vc) * Central
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Censored(truncated(Normal(cp, cp*σ_p), 0, Inf), _lloq, Inf)
    end
end

df = CSV.read("02-depot_1cmt_linear/data/single_dose.csv", DataFrame, 
              missingstring = ".")
rename!(lowercase, df)

pop = read_pumas(df, 
                 covariates = [:lloq])

iparams = (;
    TVCL = 3.7,
    TVVC = 80,
    TVKA = 1.2,
    σ_p = 0.3,
    C = float.(Matrix(I(3))),
    ω = [0.2, 0.3, 0.4]
)

pumas_fit = fit(
    depot_1cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 100,
        nadapts = 50,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true)
    )

Pumas.truncate(pumas_fit; burnin = 50)