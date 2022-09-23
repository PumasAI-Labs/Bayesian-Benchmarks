using Pumas
using CSV
using DataFrames

#########################################################################################
# HCV model from "Methods and software tools for design evaluation in population        #
# pharmacokinetics– pharmacodynamics studies" in BJCP written by Nyberg et al           #
#########################################################################################

hcv_model = @model begin
    # The "@param" block specifies the parameters
    @param begin
        # fixed effects with lower and upper bounds
        logθKa ~ Normal(log(0.8), 1)
        logθKe ~ Normal(log(0.15), 1)
        logθVd ~ Normal(log(100), 1)
        logθn ~ Normal(log(2.0), 1)
        logθδ ~ Normal(log(0.20), 1)
        logθc ~ Normal(log(7.0), 1)
        logθEC50 ~ Normal(log(0.12), 1)
        # random effects variance parameters, must be posisitive
        ω²Ka ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²Ke ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²Vd ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²n ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²δ ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²c ~ TruncatedNormal(0.25, 1, 0, Inf)
        ω²EC50 ~ TruncatedNormal(0.25, 1, 0, Inf)
        # variance parameter in proportional error model
        σ²PK ~ TruncatedNormal(0.04, 1, 0, Inf)
        σ²PD ~ TruncatedNormal(0.04, 1, 0, Inf)
    end

    # The random block allows us to specify variances for, and covariances
    # between, the random effects
    @random begin
        ηKa ~ Normal(0.0, sqrt(ω²Ka))
        ηKe ~ Normal(0.0, sqrt(ω²Ke))
        ηVd ~ Normal(0.0, sqrt(ω²Vd))
        ηn ~ Normal(0.0, sqrt(ω²n))
        ηδ ~ Normal(0.0, sqrt(ω²δ))
        ηc ~ Normal(0.0, sqrt(ω²c))
        ηEC50 ~ Normal(0.0, sqrt(ω²EC50))
    end

    @pre begin
        # constants
        p = 100.0
        d = 0.001
        e = 1e-7
        s = 20000.0

        logKa = logθKa + ηKa
        logKe = logθKe + ηKe
        logVd = logθVd + ηVd
        logn = logθn + ηn
        logδ = logθδ + ηδ
        logc = logθc + ηc
        logEC50 = logθEC50 + ηEC50
    end

    @init begin
        T = exp(logc + logδ) / (p * e)
        I = (s * e * p - d * exp(logc + logδ)) / (p * exp(logδ) * e)
        W = (s * e * p - d * exp(logc + logδ)) / (exp(logc + logδ) * e)
    end

    # The dynamics block is used to describe the evolution of our variables.
    @dynamics begin
        X' = -exp(logKa) * X
        A' = exp(logKa) * X - exp(logKe) * A
        T' = s - T * (e * W + d)
        I' = e * W * T - exp(logδ) * I
        W' =
            p / ((A / exp(logVd) / exp(logEC50) + 1e-100)^exp(logn) + 1) * I - exp(logc) * W
    end

    # The derived block is used to model the dependent variables. Both will
    # be available in our simulated data, but only `dv` has a distribution
    # here (~ read "ditributed as").
    @derived begin
        conc := @. A / exp(logVd)
        log10W := @. log10(W)
        yPK ~ @. TruncatedNormal(A / exp(logVd), sqrt(σ²PK), 0, Inf)
        yPD ~ @. TruncatedNormal(log10W, sqrt(σ²PD), 0, Inf)
    end
end

# data
df = CSV.read(joinpath(pwd(), "data", "hcv.csv"), DataFrame)
select!(df, Not(:route))

# population
pop = read_pumas(df; observations=[:yPK, :yPD])

# params
parms = init_params(hcv_model)

# Fit
hcv_fit = fit(
    hcv_model,
    pop,
    parms,
    Pumas.BayesMCMC(; nsamples=2_000, nadapts=1_000, target_accept=0.8, nchains=4, progress = false);
    diffeq_options=(; alg=Rodas5())
)

Pumas.truncate(hcv_fit; burnin=1_000)
# Wall duration = 12666.02 seconds
