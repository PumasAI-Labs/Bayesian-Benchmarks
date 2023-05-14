using Pumas
using DataFrames
using CSV

depot_1cmt_mm_prop = @model begin
    @options begin
        inplace = false
    end

    @param begin
        # TVCL ~ LogNormal(log(4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVVMAX ~ LogNormal(log(1), 1)
        TVKM ~ LogNormal(log(0.25), 1)
        TVKA ~ LogNormal(log(1), 1)
        #σ_p ~ Constrained(Normal(0, 0.5), lower = 0, upper = Inf)
        σ_p ~ truncated(Normal(0, 0.5), 0.0, Inf)
        C ~ LKJCholesky(4, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(4), Diagonal([0.4, 0.4, 0.4, 0.4].^2)),
            lower = zeros(4),
            upper = fill(Inf, 4),
            init = ones(4)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(4)) # Z in the Stan code
    end

    @covariates lloq

    @pre begin

        _lloq = lloq

        # compute the η from the ηstd
        # using lower Cholesky triangular matrix
        η = ω .* (getchol(C).L * ηstd)

        # PK parameters
        Vc = TVVC * exp(η[1])
        VMAX = TVVMAX * exp(η[2])
        KM = TVKM * exp(η[3])
        Ka = TVKA * exp(η[4])

    end

    @vars
        conc = Central/Vc
    end

    @dynamics begin
        Depot' = -Ka*Depot
        Central' = Ka*Depot -  (VMAX*conc)/(KM + conc)
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Censored(
            truncated(
                Normal(cp, cp*σ_p + 1e-10),
                0.0,
                Inf,
            ),
            _lloq,
            Inf,
        )
    end
end

df = CSV.read("04-depot_1cmt_mm/data/single_dose.csv", DataFrame, 
              missingstring = ".")
rename!(lowercase, df)

pop = read_pumas(df, 
                 covariates = [:lloq])

iparams = (;
    TVVC = 80,
    TVVMAX = 4.3,
    TVKM = 35,
    TVKA = 1.2,
    σ_p = 0.3,
    C = float.(Matrix(I(4))),
    ω = [0.2, 0.3, 0.4, 0.3]
)

pumas_fit = fit(
    depot_1cmt_mm_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true)
    )

Pumas.truncate(pumas_fit; burnin = 500)

serialize("04-depot_1cmt_mm/Pumas/fit_single_dose", my_fit)