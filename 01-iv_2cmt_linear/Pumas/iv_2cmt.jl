using Pumas
using DataFrames
using CSV
using Serialization
using JSON3

iv_2cmt_exp = @model begin
    @param begin
        TVCL ~ LogNormal(log(4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVQ ~ LogNormal(log(4), 1)
        TVVP ~ LogNormal(log(50), 1)
        σ ~ Constrained(Normal(0, 0.5); lower=0.0)
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
    end

    @dynamics Central1Periph1

    @derived begin
        cp := @. Central / Vc
        dv ~ @. LogNormal(log(cp), σ)
    end
end

df = CSV.read("01-iv_2cmt_linear/data/single_dose.csv", DataFrame,
    missingstring=".")
df_multi = CSV.read("01-iv_2cmt_linear/data/multiple_dose.csv", DataFrame,
    missingstring=".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df)
pop_multi = read_pumas(df_multi)

json_inits = filter(
    x -> endswith(x, ".json"),
    readdir("01-iv_2cmt_linear/data/inits/"; join=true)
)
# TODO: run for all 4 inits and <=10 runs (now it is 5 runs)
# for now we'll just take the first chain the other will be random
filter!(x -> contains(x, r"inits_[1|2|3|4|5]_1"), json_inits)

function parse_json(json_path; n=4)
    iparams = JSON3.read(json_path, Dict{Symbol,Any})
    delete!(iparams, :Z)
    delete!(iparams, :L)
    iparams[:omega] = Float64.(iparams[:omega])
    iparams[:ω] = iparams[:omega]
    delete!(iparams, :omega)
    iparams[:σ] = iparams[:sigma]
    delete!(iparams, :sigma)
    iparams = (;
        iparams...,
        C=I(n)
    )
    return iparams
end

iparams = map(parse_json, json_inits)

pumas_fits = map(
    p -> fit(
        iv_2cmt_exp,
        pop,
        p,
        BayesMCMC(
            nsamples=1500,
            nadapts=500,
            nchains=4,
            parallel_chains=true,
            parallel_subjects=true,
        )
    ),
    iparams
)

my_fits = map(x -> discard(x; burnin=500), pumas_fits)
map(
    (i, f) -> serialize("01-iv_2cmt_linear/Pumas/fit_single_dose_$i.jls", f),
    1:length(my_fits),
    my_fits
)

pumas_fits_multi = map(
    p -> fit(
        iv_2cmt_exp,
        pop_multi,
        p,
        BayesMCMC(
            nsamples=1500,
            nadapts=500,
            nchains=4,
            parallel_chains=true,
            parallel_subjects=true,
        )
    ),
    iparams
)

my_fits_multi = map(x -> discard(x; burnin=500), pumas_fits_multi)
map(
    (i, f) -> serialize("01-iv_2cmt_linear/Pumas/fit_multi_dose_$i.jls", f),
    1:length(my_fits_multi),
    my_fits_multi
)
