using Arrow
using MCMCChains
using DataFrames
using DataFramesMeta
using CSV
using Serialization
using Dates
using Pumas
using CairoMakie
using AlgebraOfGraphics

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
    mat = Array{Float64}(undef, Int(nrow(df) / nchains), length(names_df), nchains)
    for c in 1:nchains
        mat[:, :, c] .= Matrix(
            select(
                subset(df, Symbol(".chain") => ByRow(==(c))),
                names_df)
        )
    end
    info = (; start_time=fill(0.0, nchains), stop_time=times[:, 2])
    chn = Chains(mat, names_df, Dict(:internals => internals); info)
    return chn
end
function rename_stan(chn)
    names_chn = names(chn)
    theta_names = filter(s -> startswith(s, "TV"), string.(names_chn))
    omega_names = filter(s -> startswith(s, "omega"), string.(names_chn))
    sigma_names = filter(s -> startswith(s, "sigma"), string.(names_chn))
    new_omega_names = "ω" .* getindex.(split.(theta_names, "TV"), 2)
    max_sigma_idx = length(filter(s -> startswith(s, "sigma"), string.(names_chn)))
    if max_sigma_idx > 1
        new_sigma_names = ["σ", "σ_pd"]
    else
        new_sigma_names = "σ"
    end
    # Now we interleave them with `zip`
    names_to_change = Dict(vcat(
        zip(omega_names, new_omega_names)...,
        zip(sigma_names, new_sigma_names)...,
    ))
    return replacenames(chn, names_to_change)
end

function get_chains_nonmen(
    arrow;
    nchains=4
)
    df = DataFrame(Arrow.Table(arrow))
    names_df = setdiff(names(df), ["chain", "iteration", "time"])
    theta_names = filter(s -> startswith(s, "theta"), names_df)
    sigma_names = filter(s -> startswith(s, "sigma"), names_df)
    omega_names = filter(s -> startswith(s, "omega"), names_df)
    # only retain diagonal sigmas (WTF?!?)
    max_sigma_idx = parse(Int64, maximum(getindex.(split.(sigma_names, '_'), 2)))
    sigmas_to_retain = ["sigma_$(i)_$(i)" for i in 1:max_sigma_idx]
    sigmas_to_filter = setdiff(sigma_names, sigmas_to_retain)
    select!(df, Not(sigmas_to_filter))
    filter!(p -> p ∉ sigmas_to_filter, names_df)
    # only retain diagonal omegas
    max_omega_idx = parse(Int64, maximum(getindex.(split.(omega_names, '_'), 2)))
    omegas_to_retain = ["omega_$(i)_$(i)" for i in 1:max_omega_idx]
    omegas_to_filter = setdiff(omega_names, omegas_to_retain)
    select!(df, Not(omegas_to_filter))
    filter!(p -> p ∉ omegas_to_filter, names_df)
    # exponentiate back thetas
    transform!(df, theta_names .=> ByRow(exp); renamecols=false)
    times = combine(
        groupby(
            df,
            "chain"
        ),
        :time => first
    )
    mat = Array{Float64}(undef, Int(nrow(df) / nchains), length(names_df), nchains)
    for c in 1:nchains
        mat[:, :, c] .= Matrix(
            select(
                subset(df, Symbol("chain") => ByRow(==(c))),
                names_df)
        )
    end
    info = (; start_time=fill(0.0, nchains), stop_time=times[:, 2])
    chn = Chains(mat, names_df; info)
    return chn
end
# 01 has 4 thetas/omegas
# 02 has 3 thetas/omegas
# 03 has 5 thetas/omegas
# TODO: 04
# 05 has 9 thetas/omegas and 2 sigmas
dict_nonmem = Dict(
    4 => ["TVCL", "TVVC", "TVQ", "TVVP"],
    3 => ["TVCL", "TVVC", "TVKA"],
    5 => ["TVCL", "TVVC", "TVQ", "TVVP", "TVKA"],
    9 => [
        "TVCL",
        "TVVC",
        "TVQ",
        "TVVP",
        "TVKA",
        "TVMTT",
        "TVCIRC0",
        "TVGAMMA",
        "TVALPHA",
    ],
)
function rename_nonmem(chn)
    names_chn = names(chn)
    theta_names = filter(s -> startswith(s, "theta"), string.(names_chn))
    omega_names = filter(s -> startswith(s, "omega"), string.(names_chn))
    sigma_names = filter(s -> startswith(s, "sigma"), string.(names_chn))
    max_theta_idx = parse(Int64, maximum(last.(split.(theta_names, "theta"))))
    new_theta_names = dict_nonmem[max_theta_idx]
    new_omega_names = "ω" .* getindex.(split.(new_theta_names, "TV"), 2)
    max_sigma_idx = parse(Int64, maximum(getindex.(split.(sigma_names, '_'), 2)))
    if max_sigma_idx > 1
        new_sigma_names = ["σ", "σ_pd"]
    else
        new_sigma_names = "σ"
    end
    # Now we interleave them with `zip`
    names_to_change = Dict(vcat(
        zip(theta_names, new_theta_names)...,
        zip(omega_names, new_omega_names)...,
        zip(sigma_names, new_sigma_names)...,
    ))
    return replacenames(chn, names_to_change)
end

function rename_pumas(chn)
    names_chn = names(chn)
    theta_names = filter(s -> startswith(s, "TV"), string.(names_chn))
    omega_names = filter(s -> startswith(s, 'ω'), string.(names_chn))
    new_omega_names = "ω" .* getindex.(split.(theta_names, "TV"), 2)
    # Now we interleave them with `zip`
    names_to_change = Dict(vcat(
        zip(omega_names, new_omega_names)...,
    ))
    return replacenames(chn, names_to_change)
end
function filter_pumas(chn)
    names_chn = names(chn)
    names_to_drop = filter(s -> startswith(s, "C"), string.(names_chn))
    names_to_keep = setdiff(string.(names_chn), names_to_drop)
    return mapreduce(p -> group(chn, p), hcat, names_to_keep)
end

function _wall_duration(c::Chains; start=MCMCChains.min_start(c), stop=MCMCChains.max_stop(c))
    return Dates.value(stop - start) / 1000
end

function mean_ess_sec(chn)
    return mean_ess(chn) / _wall_duration(chn)
end
function mean_ess(chn)
    summ_df = DataFrame(summarystats(chn))
    mean_ess = mean(summ_df[:, :ess_tail])
    return mean_ess
end
function mean_rhat(chn)
    summ_df = DataFrame(summarystats(chn))
    mean_ess = mean(summ_df[:, :rhat])
    return mean_ess
end

# Stan
# stan_files_01_sd = filter(f -> contains(f, r"single.*\d.arrow"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))
stan_files_01_md = filter(f -> contains(f, r"multi.*\d.arrow"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))
# stan_files_02_sd = filter(f -> contains(f, r"single.*\d.arrow"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "Stan", "Torsten", "Fits"); join=true))
stan_files_02_md = filter(f -> contains(f, r"multi.*\d.arrow"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "Stan", "Torsten", "Fits"); join=true))
# stan_files_03_sd = filter(f -> contains(f, r"single.*\d.arrow"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))
stan_files_03_md = filter(f -> contains(f, r"multi.*\d.arrow"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "Stan", "Torsten", "Fits"); join=true))
# stan_files_04_sd = filter(f -> contains(f, r"single.*\d.arrow"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "Stan", "Torsten", "Fits"); join=true))
# stan_files_04_md = filter(f -> contains(f, r"multi.*\d.arrow"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "Stan", "Torsten", "Fits"); join=true))
stan_files_05 = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "05-friberg", "Stan", "Torsten", "Fits"); join=true))

# stan_chn_01_sd = rename_stan.(get_chains_stan.(stan_files_01_sd))
stan_chn_01_md = rename_stan.(get_chains_stan.(stan_files_01_md))
# stan_chn_02_sd = rename_stan.(get_chains_stan.(stan_files_02_sd))
stan_chn_02_md = rename_stan.(get_chains_stan.(stan_files_02_md))
# stan_chn_03_sd = rename_stan.(get_chains_stan.(stan_files_03_sd))
stan_chn_03_md = rename_stan.(get_chains_stan.(stan_files_03_md))
# stan_chn_04_sd = rename_stan.(get_chains_stan.(stan_files_04_sd))
# stan_chn_04_md = rename_stan.(get_chains_stan.(stan_files_04_md))
stan_chn_05 = rename_stan.(get_chains_stan.(stan_files_05))

stan_df = DataFrame(;
    model=[
        "01_mult_dose",
        "02_mult_dose",
        "03_mult_dose",
        #  "04_mult_dose",
        "05",
    ],
    mean_ess=map(c -> mean(mean_ess.(c)),
        [
            stan_chn_01_md,
            stan_chn_02_md,
            stan_chn_03_md,
            # stan_chn_04_md,
            stan_chn_05,
        ]
    ),
    mean_ess_sec=map(c -> mean(mean_ess_sec.(c)),
        [
            stan_chn_01_md,
            stan_chn_02_md,
            stan_chn_03_md,
            #  stan_chn_04_md,
            stan_chn_05,
        ]
    ),
    mean_rhat=map(c -> mean(mean_rhat.(c)),
        [
            stan_chn_01_md,
            stan_chn_02_md,
            stan_chn_03_md,
            # stan_chn_04_md,
            stan_chn_05,
        ]
    )
)

CSV.write("results/stan.csv", stan_df)

# NONMEM
# nonmem_files_01_sd = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "NONMEM", "iv-2cmt-linear", "chains"); join=true))
nonmem_files_01_md = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "NONMEM", "iv-2cmt-linear-md", "chains"); join=true))
# nonmem_files_02_sd = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "NONMEM", "depot-1cmt-linear", "chains"); join=true))
nonmem_files_02_md = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "NONMEM", "depot-1cmt-linear-md", "chains"); join=true))
# nonmem_files_03_sd = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "NONMEM", "depot-2cmt-linear", "chains"); join=true))
nonmem_files_03_md = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "NONMEM", "depot-2cmt-linear-md", "chains"); join=true))
# nonmem_files_04_sd = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "NONMEM", "depot-1cmt-mm", "chains"); join=true))
# nonmem_files_04_md = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "NONMEM", "depot-1cmt-mm-md", "chains"); join=true))
nonmem_files_05 = filter(f -> contains(f, r"\d.arrow"), readdir(joinpath(pwd(), "05-friberg", "NONMEM", "friberg", "chains"); join=true))

# nonmem_chn_01_sd = rename_nonmem.(get_chains_nonmen.(nonmem_files_01_sd))
nonmem_chn_01_md = rename_nonmem.(get_chains_nonmen.(nonmem_files_01_md))
# nonmem_chn_02_sd = rename_nonmem.(get_chains_nonmen.(nonmem_files_02_sd))
nonmem_chn_02_md = rename_nonmem.(get_chains_nonmen.(nonmem_files_02_md))
# nonmem_chn_03_sd = rename_nonmem.(get_chains_nonmen.(nonmem_files_03_sd))
nonmem_chn_03_md = rename_nonmem.(get_chains_nonmen.(nonmem_files_03_md))
# nonmem_chn_04_sd = rename_nonmem.(get_chains_nonmen.(nonmem_files_04_sd))
# nonmem_chn_04_md = rename_nonmem.(get_chains_nonmen.(nonmem_files_04_md))
nonmem_chn_05 = rename_nonmem.(get_chains_nonmen.(nonmem_files_05))

nonmem_df = DataFrame(;
    model=[
        "01_mult_dose",
        "02_mult_dose",
        "03_mult_dose",
        #  "04_mult_dose",
        "05",
    ],
    mean_ess=map(c -> mean(mean_ess.(c)),
        [
            nonmem_chn_01_md,
            nonmem_chn_02_md,
            nonmem_chn_03_md,
            # nonmem_chn_04_md,
            nonmem_chn_05,
        ]
    ),
    mean_ess_sec=map(c -> mean(mean_ess_sec.(c)),
        [
            nonmem_chn_01_md,
            nonmem_chn_02_md,
            nonmem_chn_03_md,
            # nonmem_chn_04_md,
            nonmem_chn_05,
        ]
    ),
    mean_rhat=map(c -> mean(mean_rhat.(c)),
        [
            nonmem_chn_01_md,
            nonmem_chn_02_md,
            nonmem_chn_03_md,
            # nonmem_chn_04_md,
            nonmem_chn_05,
        ]
    )
)

CSV.write("results/nonmem.csv", nonmem_df)

# Pumas
# Stan
# pumas_files_01_sd = filter(f -> contains(f, r"single.*\d.jls"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "Pumas"); join=true))
pumas_files_01_md = filter(f -> contains(f, r"multi.*\d.jls"), readdir(joinpath(pwd(), "01-iv_2cmt_linear", "Pumas"); join=true))
# pumas_files_02_sd = filter(f -> contains(f, r"single.*\d.jls"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "Pumas"); join=true))
pumas_files_02_md = filter(f -> contains(f, r"multi.*\d.jls"), readdir(joinpath(pwd(), "02-depot_1cmt_linear", "Pumas"); join=true))
# pumas_files_03_sd = filter(f -> contains(f, r"single.*\d.jls"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "Pumas"); join=true))
pumas_files_03_md = filter(f -> contains(f, r"multi.*\d.jls"), readdir(joinpath(pwd(), "03-depot_2cmt_linear", "Pumas"); join=true))
# pumas_files_04_sd = filter(f -> contains(f, r"single.*\d.jls"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "Pumas"); join=true))
# pumas_files_04_md = filter(f -> contains(f, r"multi.*\d.jls"), readdir(joinpath(pwd(), "04-depot_1cmt_mm", "Pumas"); join=true))
pumas_files_05 = filter(f -> contains(f, r"\d.jls"), readdir(joinpath(pwd(), "05-friberg", "Pumas"); join=true))

# pumas_chn_01_sd = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_01_sd))))
pumas_chn_01_md = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_01_md))))
# pumas_chn_02_sd = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_02_sd))))
pumas_chn_02_md = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_02_md))))
# pumas_chn_03_sd = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_03_sd))))
pumas_chn_03_md = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_03_md))))
# pumas_chn_04_sd = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_04_sd))))
# pumas_chn_04_md = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_04_md))))
pumas_chn_05 = filter_pumas.(rename_pumas.(Chains.(deserialize.(pumas_files_05))))

pumas_df = DataFrame(;
    model=[
        "01_mult_dose",
        "02_mult_dose",
        "03_mult_dose",
        # "04_mult_dose",
        "05",
    ],
    mean_ess=map(c -> mean(mean_ess.(c)),
        [
            pumas_chn_01_md,
            pumas_chn_02_md,
            pumas_chn_03_md,
            # pumas_chn_04_md,
            pumas_chn_05,
        ]
    ),
    mean_ess_sec=map(c -> mean(mean_ess_sec.(c)),
        [
            pumas_chn_01_md,
            pumas_chn_02_md,
            pumas_chn_03_md,
            # pumas_chn_04_md,
            pumas_chn_05,
        ]
    ),
    mean_rhat=map(c -> mean(mean_rhat.(c)),
        [
            pumas_chn_01_md,
            pumas_chn_02_md,
            pumas_chn_03_md,
            # pumas_chn_04_md,
            pumas_chn_05,
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

summ_df = @chain all_df begin
    groupby(:software)
    @combine begin
        :mean_ess = mean(:mean_ess)
        :mean_ess_sec = mean(:mean_ess_sec)
        :mean_rhat = mean(:mean_rhat)
    end
end

CSV.write("results/summary.csv", summ_df)

# Plots
plt = data(all_df) *
      mapping(
          :model => renamer([
              "01_mult_dose" => "1",
              "02_mult_dose" => "2",
              "03_mult_dose" => "3",
              "05" => "5"
          ]),
          [:mean_ess, :mean_ess_sec, :mean_rhat];
          color=:software,
          dodge=:software,
          col=dims(1) => renamer(["Mean ESS", "Mean ESS/sec", "Mean Rhat"])
      ) *
      visual(BarPlot)

fig = draw(
    plt;
    axis=(;
        xticklabelrotation=π / 3,
        ylabel=""
    ),
    facet=(; linkyaxes=:none)
)

save("results/results.png", fig; px_per_unit=3)