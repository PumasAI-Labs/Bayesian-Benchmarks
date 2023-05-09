using Pumas
using DataFrames
using CSV

depot_2cmt_prop = @model begin

    @param begin
        TVCL ~ LogNormal(log(4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVQ ~ LogNormal(log(4), 1)
        TVVP ~ LogNormal(log(40), 1)
        TVKA ~ LogNormal(log(1), 1)
        #σ_p ~ Constrained(Normal(0, 0.5), lower = 0, upper = Inf)
        σ_p ~ truncated(Normal(0, 0.5), 0, Inf)
        C ~ LKJCholesky(5, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(5), Diagonal([0.4, 0.4, 0.4, 0.4, 0.4].^2)),
            lower = zeros(5),
            upper = fill(Inf, 5),
            init = ones(5)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(5)) # Z in the Stan code
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
        Q = TVCL * exp(η[3])
        Vp = TVVC * exp(η[4])
        Ka = TVKA * exp(η[5])

        k_cp = Q/Vc
        k_pc = Q/Vp
    end

    @dynamics begin
        Depot' = -Ka * Depot
        Central' = Ka * Depot - (CL/Vc + k_cp) * Central + k_pc*Peripheral
        Peripheral' = k_cp * Central - k_pc*Peripheral
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Censored(truncated(Normal(cp, cp*σ_p), 0, Inf), _lloq, Inf)
    end
end

df = CSV.read("03-depot_2cmt_linear/data/single_dose.csv", DataFrame, 
              missingstring = ".")
df_multi = CSV.read("03-depot_2cmt_linear/data/multiple_dose.csv", DataFrame, 
              missingstring = ".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df, 
                 covariates = [:lloq])
pop_multi = read_pumas(df_multi, 
                 covariates = [:lloq])

iparams = (;
    TVCL = exp(1.2970),
    TVVC = exp(3.8210),
    TVQ = exp(1.3210),
    TVVP = exp(3.6870),
    TVKA = exp(-0.2577),
    σ_p = 0.1957,
    C = float.(Matrix(I(5))),
    ω = [0.1600, 0.2770, 0.2909, 0.2360, 0.3500]
)

pumas_fit = fit(
    depot_2cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true,
        max_chunk_size=16,
        )
    )

my_fit = Pumas.truncate(pumas_fit; burnin = 500)
serialize("03-depot_2cmt_linear/Pumas/fit_single_dose.jls", my_fit)

pumas_fit_multi = fit(
    depot_2cmt_prop,
    pop_multi,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true,
        max_chunk_size=16,
        )
    )

my_fit_multi = Pumas.truncate(pumas_fit_multi; burnin = 500)
serialize("03-depot_2cmt_linear/Pumas/fit_multi_dose.jls", my_fit_multi)