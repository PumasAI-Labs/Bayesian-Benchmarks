using Arrow
using MCMCChains
using DataFrames
using CSV

function get_chains_stan(
    arrow;
    nchains=4,
    internals=Symbol.(
        ["lp__"]
    )
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

function get_chains_nonmen(files, lst)
    dfs = CSV.read.(files, DataFrame; delim='\t', skipto=3, header=2)
    df = vcat(dfs...)
    # TODO get times from lst file using readlines
    # grep("Elapsed estimation  time in seconds: ", lstText)
end

function mean_ess_sec(chn)
    summ_df = summarystats(chn)
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