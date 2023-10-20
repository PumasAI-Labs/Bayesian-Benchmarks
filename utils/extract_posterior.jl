using Arrow
using MCMCChains
using DataFrames
using DataFramesMeta
using CSV
using Serialization
using Dates
using Pumas

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

# More details per model
function extract_chn_model(chns)
    parameters = DataFrame(summarystats(first(chns)))[:, :parameters]
    df = mapreduce(p -> _summary_params(chns, p), vcat, parameters)
    return df
end
function _get_params(chn::Chains, p::Symbol)
    df = @chain DataFrame(summarystats(chn)) begin
        @rsubset :parameters == p
    end
    return Vector(df[1, :])[2:end]
end
function _summary_params(chns, p::Symbol)
    names_df = names(DataFrame(summarystats(first(chns))))
    flat = Iterators.flatten(_get_params.(chns, p)) |> collect
    summary_params = reshape(flat, :, 5)
    means = reshape(mean(summary_params; dims=2), 1, :)
    summ_df = DataFrame(means, names_df[2:end]; copycols=true)
    @rtransform! summ_df :parameters = p
    select!(summ_df, :parameters, :)
    return summ_df
end
function save_chn_model(chns, model, software)
    CSV.write(joinpath(pwd(), "results", software, "$(model).csv"), extract_chn_model(chns))
end

models_names = [
    "01-md",
    "02-md",
    "03-md",
    "05",
]
# Stan
chains_stan = [
    stan_chn_01_md,
    stan_chn_02_md,
    stan_chn_03_md,
    stan_chn_05,
]
map((c, m) -> save_chn_model(c, m, "stan"), chains_stan, models_names)
# NONMEM
chains_nonmem = [
    nonmem_chn_01_md,
    nonmem_chn_02_md,
    nonmem_chn_03_md,
    nonmem_chn_05,
]
map((c, m) -> save_chn_model(c, m, "nonmem"), chains_nonmem, models_names)
# Pumas
chains_pumas = [
    pumas_chn_01_md,
    pumas_chn_02_md,
    pumas_chn_03_md,
    pumas_chn_05,
]
map((c, m) -> save_chn_model(c, m, "pumas"), chains_pumas, models_names)

# Rhats within across software implementations
#= WIP
function merge_chain(chn)
    df = DataFrame(chn)
    @rtransform! df :chain = 1
    select!(df, Not(:chain, :iteration))
    if "lp__" in names(df)
        select!(df, Not("lp__"))
    end
    return Chains(Matrix(df), names(df))
end
function merge_chains(chns)
    return mapreduce(merge_chain, chainscat, chns)
end
function merge_chains_software(vec_of_chns)
    return mapreduce(merge_chains, chainscat, vec_of_chns)
end
merge_chains(stan_chn_01_md)
merge_chains_software(chains_01)
=#