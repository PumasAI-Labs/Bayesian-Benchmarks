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

depot_1cmt_mm_exp = @model begin
    @param begin
        TVVC ~ LogNormal(log(70), 1)
        TVVMAX ~ LogNormal(log(1), 1)
        TVKM ~ LogNormal(log(0.25), 1)
        TVKA ~ LogNormal(log(1), 1)
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
        Vc = TVVC * exp(η[1])
        VMAX = TVVMAX * exp(η[2])
        KM = TVKM * exp(η[3])
        Ka = TVKA * exp(η[4])
    end

    @vars begin
        conc = Central / Vc
    end

    @dynamics begin
        Depot' = -Ka * Depot
        Central' = Ka * Depot - (VMAX * conc) / (KM + conc)
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. LogNormal(log(cp), σ)
    end
end

df = CSV.read("04-depot_1cmt_mm/data/single_dose.csv", DataFrame,
    missingstring=".")
df_multi = CSV.read("04-depot_1cmt_mm/data/multiple_dose.csv", DataFrame,
    missingstring=".")
rename!(lowercase, df)
rename!(lowercase, df_multi)

pop = read_pumas(df)
pop_multi = read_pumas(df_multi)

json_inits = filter(
    x -> endswith(x, ".json"),
    readdir("04-depot_1cmt_mm/data/inits/"; join=true)
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

# dummy fit to trigger precompilation
fit(
    depot_1cmt_mm_exp,
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
        depot_1cmt_mm_exp,
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
    (i, f) -> serialize("04-depot_1cmt_mm/Pumas/fit_single_dose_$i.jls", f),
    1:length(my_fits),
    my_fits
)

# dummy fit to trigger precompilation
fit(
    depot_1cmt_mm_exp,
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
        depot_1cmt_mm_exp,
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
    (i, f) -> serialize("04-depot_1cmt_mm/Pumas/fit_multi_dose_$i.jls", f),
    1:length(my_fits_multi),
    my_fits_multi
)
