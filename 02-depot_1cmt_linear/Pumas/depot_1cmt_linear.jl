using Pumas
using DataFrames
using CSV
using Serialization

depot_1cmt_prop = @model begin

    @param begin
        TVCL ~ LogNormal(log(4), 1)
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
df_multi = CSV.read("02-depot_1cmt_linear/data/multiple_dose.csv", DataFrame, 
              missingstring = ".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df, 
                 covariates = [:lloq])
pop_multi = read_pumas(df_multi, 
                 covariates = [:lloq])

iparams = (;
    TVCL = exp(1.29700),
    TVVC = exp(3.82100),
    TVKA = exp(-0.06489),
    σ_p = 0.18470,
    C = float.(Matrix(I(3))),
    ω = [0.29950, 0.23180, 0.16000]
)

pumas_fit = fit(
    depot_1cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true)
    )

my_fit = Pumas.truncate(pumas_fit; burnin = 500)
serialize("02-depot_1cmt_linear/Pumas/fit_single_dose.jls", my_fit)

pumas_fit_multi = fit(
    depot_1cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true)
    )

my_fit_multi = Pumas.truncate(pumas_fit_multi; burnin = 500)
serialize("02-depot_1cmt_linear/Pumas/fit_multi_dose.jls", my_fit_multi)