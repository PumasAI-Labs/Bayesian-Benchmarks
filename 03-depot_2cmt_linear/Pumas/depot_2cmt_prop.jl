using Pumas
using DataFrames
using CSV

depot_2cmt_prop = @model begin

    @param begin
        TVCL ~ LogNormal(log(0.4), 1)
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
rename!(lowercase, df)

pop = read_pumas(df, 
                 covariates = [:lloq])

iparams = (;
    TVCL = 3.7,
    TVVC = 80,
    TVQ = 4.3,
    TVVP = 35,
    TVKA = 1.2,
    σ_p = 0.3,
    C = float.(Matrix(I(5))),
    ω = [0.2, 0.3, 0.4, 0.3, 0.2]
)

pumas_fit = fit(
    depot_2cmt_prop,
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