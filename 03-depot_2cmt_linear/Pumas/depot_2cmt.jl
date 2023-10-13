using Distributed
if nprocs() < 5
    addprocs(
        5 - nprocs(),
        exeflags=["--threads=$(Threads.nthreads())", "--project=$(Base.active_project())"],
    )
end
@everywhere using Pumas
using DataFrames
using CSV
using Serialization
using JSON3

depot_2cmt_exp = @model begin
    @param begin
        TVCL ~ LogNormal(log(4), 1)
        TVVC ~ LogNormal(log(70), 1)
        TVQ ~ LogNormal(log(4), 1)
        TVVP ~ LogNormal(log(40), 1)
        TVKA ~ LogNormal(log(1), 1)
        σ ~ Constrained(Normal(0, 0.5), lower=0.0)
        C ~ LKJCholesky(5, 2) # L in the Stan code is the lower triangular part of the Cholesky decomposition
        ω ∈ Constrained(
            MvNormal(zeros(5), Diagonal([0.4, 0.4, 0.4, 0.4, 0.4] .^ 2)),
            lower=zeros(5),
            init=ones(5)
        )
    end

    @random begin
        ηstd ~ MvNormal(I(5)) # Z in the Stan code
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
        Ka = TVKA * exp(η[5])
    end

    @dynamics Depots1Central1Periph1

    @derived begin
        cp := @. Central / Vc
        dv ~ @. LogNormal(log(cp), σ)
    end
end

df = CSV.read("03-depot_2cmt_linear/data/single_dose.csv", DataFrame,
    missingstring=".")
df_multi = CSV.read("03-depot_2cmt_linear/data/multiple_dose.csv", DataFrame,
    missingstring=".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df)
pop_multi = read_pumas(df_multi)

json_inits = filter(
    x -> endswith(x, ".json"),
    readdir("03-depot_2cmt_linear/data/inits/"; join=true)
)
# TODO: run for all 4 inits and <=10 runs (now it is 5 runs)
# for now we'll just take the first chain the other will be random
filter!(x -> contains(x, r"inits_[1|2|3|4|5]_1"), json_inits)

function parse_json(json_path; n=5)
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

# dummy fit to trigger precompilation
fit(
    depot_2cmt_exp,
    pop,
    iparams[1],
    BayesMCMC(
        nsamples=10,
        nadapts=5,
        nchains=4,
        parallel_chains=true,
        parallel_subjects=true,
        ensemblealg=EnsembleSplitThreads(),
    )
)

pumas_fits = map(
    p -> fit(
        depot_2cmt_exp,
        pop,
        p,
        BayesMCMC(
            nsamples=1500,
            nadapts=500,
            nchains=4,
            parallel_chains=true,
            parallel_subjects=true,
            ensemblealg=EnsembleSplitThreads(),
        )
    ),
    iparams
)

my_fits = map(x -> discard(x; burnin=500), pumas_fits)
map(
    (i, f) -> serialize("03-depot_2cmt_linear/Pumas/fit_single_dose_$i.jls", f),
    1:length(my_fits),
    my_fits
)

# dummy fit to trigger precompilation
fit(
    depot_2cmt_exp,
    pop_multi,
    iparams[1],
    BayesMCMC(
        nsamples=10,
        nadapts=5,
        nchains=4,
        parallel_chains=true,
        parallel_subjects=true,
        ensemblealg=EnsembleSplitThreads(),
    )
)

pumas_fits_multi = map(
    p -> fit(
        depot_2cmt_exp,
        pop_multi,
        p,
        BayesMCMC(
            nsamples=1500,
            nadapts=500,
            nchains=4,
            parallel_chains=true,
            parallel_subjects=true,
            ensemblealg=EnsembleSplitThreads(),
        )
    ),
    iparams
)

my_fits_multi = map(x -> discard(x; burnin=500), pumas_fits_multi)
map(
    (i, f) -> serialize("03-depot_2cmt_linear/Pumas/fit_multi_dose_$i.jls", f),
    1:length(my_fits_multi),
    my_fits_multi
)
