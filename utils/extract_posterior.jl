using Arrow
using MCMCChains
using DataFrames
using DataFramesMeta
using CSV
using Serialization

function get_chains_stan(
    arrow;
    nchains=4,
    internals=Symbol.(
        ["lp__"]
    ),
)
    df = DataFrame(Arrow.Table(arrow))
    names_df = setdiff(names(df), [".chain", ".iteration", ".draw", "warmup", "sampling", "total"])
    times = combine(
        groupby(
            df,
            ".chain"
        ),
        :total => first
    )
    mat = Array{Float64}(undef, Int(nrow(df)/nchains), length(names_df), nchains)
    for c in 1:nchains
        mat[:, :, c ] .= Matrix(
            select(
                subset(df, Symbol(".chain") => ByRow(==(c))),
            names_df)
        )
    end
    info = (; start_time=fill(0.0, nchains), stop_time=times[:, 2])
    chn = Chains(mat, names_df, Dict(:internals => internals); info)
    return chn
end

function get_chains_nonmen(
        arrow;
        nchains=4,
    )
    df = DataFrame(Arrow.Table(arrow))
    names_df = setdiff(names(df), ["chain", "iteration", "time"])
    times = combine(
        groupby(
            df,
            "chain"
        ),
        :time => first
    )
    mat = Array{Float64}(undef, Int(nrow(df)/nchains), length(names_df), nchains)
    for c in 1:nchains
        mat[:, :, c ] .= Matrix(
            select(
                subset(df, Symbol("chain") => ByRow(==(c))),
            names_df)
        )
    end
    info = (; start_time=fill(0.0, nchains), stop_time=times[:, 2])
    chn = Chains(mat, names_df; info)
    return chn
end

function mean_ess_sec(chn)
    summ_df = DataFrame(summarystats(chn))
    filter!(r -> !(isnan(r.ess_per_sec)), summ_df)
    mean_ess_sec = mean(summ_df[:, :ess_per_sec])
    return mean_ess_sec
end

files_01 = filter(f -> endswith(f, ".arrow"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))
files_02 = filter(f -> endswith(f, ".arrow"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "Stan", "Torsten", "Fits"); join=true))
files_03 = filter(f -> endswith(f, ".arrow"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))

chn_01_mult_dose, chn_01_mult_dose_mat_exp, chn_01_single_dose, chn_01_single_dose_mat_exp = get_chains_stan.(files_01)
chn_02_mult_dose, chn_02_mult_dose_mat_exp, chn_02_single_dose, chn_02_single_dose_mat_exp = get_chains_stan.(files_02)
chn_03_mult_dose, chn_03_mult_dose_mat_exp, chn_03_single_dose, chn_03_single_dose_mat_exp = get_chains_stan.(files_03)

stan_df = DataFrame(;
    model=[
    "01_mult_dose","01_mult_dose_mat_exp","01_single_dose","01_single_dose_mat_exp",
    "02_mult_dose","02_mult_dose_mat_exp","02_single_dose","02_single_dose_mat_exp",
    "03_mult_dose","03_mult_dose_mat_exp","03_single_dose","03_single_dose_mat_exp",
    ],
    mean_ess_sec=mean_ess_sec.(
        [
            chn_01_mult_dose, chn_01_mult_dose_mat_exp, chn_01_single_dose, chn_01_single_dose_mat_exp,
            chn_02_mult_dose, chn_02_mult_dose_mat_exp, chn_02_single_dose, chn_02_single_dose_mat_exp,
            chn_03_mult_dose, chn_03_mult_dose_mat_exp, chn_03_single_dose, chn_03_single_dose_mat_exp,
        ]
    )
)

CSV.write("results/stan.csv", stan_df)

nonmem_01_sd = get_chains_nonmen(joinpath(pwd(), "01-iv_2cmt_linear", "NONMEM", "iv-2cmt-linear", "chains", "chains.arrow"))
nonmem_01_md = get_chains_nonmen(joinpath(pwd(), "01-iv_2cmt_linear", "NONMEM", "iv-2cmt-linear-md", "chains", "chains.arrow"))
nonmem_02_sd = get_chains_nonmen(joinpath(pwd(), "02-depot_1cmt_linear", "NONMEM", "depot-1cmt-linear", "chains", "chains.arrow"))
nonmem_02_md = get_chains_nonmen(joinpath(pwd(), "02-depot_1cmt_linear", "NONMEM", "depot-1cmt-linear-md", "chains", "chains.arrow"))
nonmem_03_sd = get_chains_nonmen(joinpath(pwd(), "03-depot_2cmt_linear", "NONMEM", "depot-2cmt-linear", "chains", "chains.arrow"))
nonmem_03_md = get_chains_nonmen(joinpath(pwd(), "03-depot_2cmt_linear", "NONMEM", "depot-2cmt-linear-md", "chains", "chains.arrow"))

nonmem_df = DataFrame(;
    model=[
    "01_single_dose","01_mult_dose_",
    "02_single_dose","02_mult_dose",
    "03_single_dose","03_mult_dose",
    ],
    mean_ess_sec=mean_ess_sec.(
        [
            nonmem_01_sd, nonmem_01_md, nonmem_02_sd, nonmem_02_md, nonmem_03_sd, nonmem_03_md
        ]
    )
)

CSV.write("results/nonmem.csv", nonmem_df)

pumas_01_sd = deserialize(joinpath(pwd(), "01-iv_2cmt_linear", "Pumas", "fit_single_dose.jls"))
pumas_01_md = deserialize(joinpath(pwd(), "01-iv_2cmt_linear", "Pumas", "fit_multi_dose.jls"))
pumas_02_sd = deserialize(joinpath(pwd(), "02-depot_1cmt_linear", "Pumas", "fit_single_dose.jls"))
pumas_02_md = deserialize(joinpath(pwd(), "02-depot_1cmt_linear", "Pumas", "fit_multi_dose.jls"))
pumas_03_sd = deserialize(joinpath(pwd(), "03-depot_2cmt_linear", "Pumas", "fit_single_dose.jls"))
pumas_03_md = deserialize(joinpath(pwd(), "03-depot_2cmt_linear", "Pumas", "fit_multi_dose.jls"))

pumas_df = DataFrame(;
    model=[
    "01_single_dose","01_mult_dose_",
    "02_single_dose","02_mult_dose",
    "03_single_dose","03_mult_dose",
    ],
    mean_ess_sec=mean_ess_sec.(
        [
            pumas_01_sd, pumas_01_md, pumas_02_sd, pumas_02_md, pumas_03_sd, pumas_03_md
        ]
    )
)

CSV.write("results/pumas.csv", pumas_df)

@rtransform! stan_df :software = "stan"
@rtransform! nonmem_df :software = "nonmem"
@rtransform! pumas_df :software = "pumas"

all_df = vcat(stan_df, nonmem_df, pumas_df)
select!(all_df, :software, All())

CSV.write("results/all.csv", all_df)
