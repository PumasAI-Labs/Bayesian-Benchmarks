using Pumas
using DataFrames
using CSV

depot_2cmt_friberg_prop = @model begin
    @options begin
        inplace = false
    end

    @param begin
        TVCL ~ LogNormal(log(0.4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVQ ~ LogNormal(log(4), 1)
        TVVP ~ LogNormal(log(40), 1)
        TVKA ~ LogNormal(log(1), 1)
        TVMTT ~ LogNormal(log(125), 1)
        TVCIRC0 ~ LogNormal(log(5), 1)
        TVGAMMA ~ LogNormal(log(0.17), 1)
        TVALPHA ~ LogNormal(log(3e-4), 1)
        σ_p ~ Constrained(Normal(0, 0.5); lower=0.0)
        σ_p_pd ~ Constrained(Normal(0, 0.5); lower=0.0)
        C ~ LKJCholesky(9, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(9), Diagonal([0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4].^2)),
            lower = zeros(9),
            init = ones(9)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(9)) # Z in the Stan code
    end

    @pre begin
        # compute the η from the ηstd
        # using lower Cholesky triangular matrix
        η = ω .* (getchol(C).L * ηstd)

        # PK parameters
        CL = TVCL * exp(η[1])
        VC = TVVC * exp(η[2])
        Q = TVQ * exp(η[3])
        VP = TVVP * exp(η[4])
        KA = TVKA * exp(η[5])
        MTT = TVMTT * exp(η[6])
        CIRC0 = TVCIRC0 * exp(η[7])
        γ = TVGAMMA * exp(η[8])
        α = TVALPHA * exp(η[9])

        k_cp = Q/VC
        k_pc = Q/VP
        k_tr = 4/MTT
    end

    @init begin
        Prol     = CIRC0
        Transit1 = CIRC0
        Transit2 = CIRC0
        Transit3 = CIRC0
        Circ     = CIRC0
    end

    @vars begin
        conc  = Central/VC
        e_drug = min(α*conc, 1)
    end

    @dynamics begin
        Depot' = -KA * Depot
        Central' = KA * Depot - (CL/VC + k_cp) * Central + k_pc*Peripheral
        Peripheral' = k_cp * Central - k_pc*Peripheral

        Prol'       = k_tr*Prol*(1 - e_drug)*abs(CIRC0/Circ)^γ - k_tr*Prol
        Transit1'   = k_tr*Prol     - k_tr*Transit1
        Transit2'   = k_tr*Transit1 - k_tr*Transit2
        Transit3'   = k_tr*Transit2 - k_tr*Transit3
        Circ'       = k_tr*Transit3 - k_tr*Circ
    end

    @derived begin
        cp := @. Central / VC
        dv ~ @. Normal(cp, cp*σ_p)
        e  ~ @. Normal(Circ, Circ*σ_p_pd)
    end
end

df = CSV.read("05-friberg/data/multiple_dose_pumas.csv", DataFrame, 
              missingstring = ".")
rename!(lowercase, df)

pop = read_pumas(
    df;
    observations=[:dv, :e],
)

iparams = (;
    TVCL = 3.7,
    TVVC = 80,
    TVQ = 4.3,
    TVVP = 35,
    TVKA = 1.2,
    TVMTT = 130,
    TVCIRC0 = 6,
    TVGAMMA = 0.15,
    TVALPHA = 2.8e-4,
    σ_p = 0.3,
    σ_p_pd = 0.3,
    C = float.(Matrix(I(9))),
    ω = [0.2, 0.3, 0.4, 0.3, 0.2, 0.3, 0.4, 0.3, 0.2]
)

pumas_fit = fit(
    depot_2cmt_friberg_prop,
    pop,
    iparams,
    BayesMCMC(
        nsamples = 1500,
        nadapts = 500,
        nchains = 4,
        parallel_chains = true,
        parallel_subjects = true)
    )

my_fit = discard(pumas_fit; burnin=500)

serialize("05-friberg/Pumas/fit_multiple_dose", my_fit)