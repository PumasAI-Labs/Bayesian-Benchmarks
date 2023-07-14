using Pumas
using DataFrames
using CSV
using Serialization

iv_2cmt_prop = @model begin
    @param begin
        TVCL ~ LogNormal(log(4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVQ ~ LogNormal(log(4), 1)
        TVVP ~ LogNormal(log(50), 1)
        σ_p ~ Constrained(Normal(0, 0.5); lower=0.0)
        C ~ LKJCholesky(4, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(4), Diagonal([0.4, 0.4, 0.4, 0.4] .^ 2)),
            lower=zeros(4),
            init=ones(4)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(4)) # Z in the Stan code
    end

    @pre begin

        # compute the η from the ηstd
        # using lower Cholesky triangular matrix
        η = ω .* (getchol(C).L * ηstd)

        # PK parameters
        CL = TVCL * exp(η[1])
        Vc = TVVC * exp(η[2])
        Q = TVQ * exp(η[3])
        Vp = TVVP * exp(η[4])

        k_cp = Q / Vc
        k_pc = Q / Vp
    end

    @dynamics begin
        Central' = -(CL / Vc + k_cp) * Central + k_pc * Peripheral
        Peripheral' = k_cp * Central - k_pc * Peripheral
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. truncated(Normal(cp, cp * σ_p); lower=0.0)
    end
end

df = CSV.read("01-iv_2cmt_linear/data/single_dose.csv", DataFrame,
    missingstring=".")
df_multi = CSV.read("01-iv_2cmt_linear/data/multiple_dose.csv", DataFrame,
    missingstring=".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df,
    covariates=[:lloq])
pop_multi = read_pumas(df_multi,
    covariates=[:lloq])

iparams = (;
    TVCL=exp(1.2970),
    TVVC=exp(3.8210),
    TVQ=exp(1.3210),
    TVVP=exp(3.9100),
    σ_p=0.1573,
    C=float.(Matrix(I(4))),
    ω=[0.2318, 0.1600, 0.2770, 0.2909]
)

pumas_fit = fit(
    iv_2cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples=1500,
        nadapts=500,
        nchains=4,
        parallel_chains=true,
        parallel_subjects=true,
        max_chunk_size=16,
    )
)

my_fit = Pumas.discard(pumas_fit; burnin=500)
serialize("01-iv_2cmt_linear/Pumas/fit_single_dose.jls", my_fit)

pumas_fit_multi = fit(
    iv_2cmt_prop,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples=1500,
        nadapts=500,
        nchains=4,
        parallel_chains=true,
        parallel_subjects=true,
        max_chunk_size=16,
    )
)

my_fit_multi = Pumas.discard(pumas_fit_multi; burnin=500)
serialize("01-iv_2cmt_linear/Pumas/fit_multi_dose.jls", my_fit_multi)

# io = IOBuffer()
# versioninfo(io)
# write("01-iv_2cmt_linear/Pumas/juliahub_versioninfo.txt", String(take!(io)))
