include("extract_posterior.jl")
using CairoMakie
using AlgebraOfGraphics

colors = [
    colorant"#008579", # nonmem
    colorant"#0035c7", # pumas 
    colorant"#b2001d", # stan
]

# Combined chains
function combine_chains(chns)
    return reduce(chainscat, chns)
end
stan_chn_01_md, stan_chn_02_md, stan_chn_03_md, stan_chn_05 = map(combine_chains, chains_stan)
nonmem_chn_01_md, nonmem_chn_02_md, nonmem_chn_03_md, nonmem_chn_05 = map(combine_chains, chains_nonmem)
pumas_chn_01_md, pumas_chn_02_md, pumas_chn_03_md, pumas_chn_05 = map(combine_chains, chains_pumas)

# Stan
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
CSV.write("results/stan/stan.csv", stan_df)

# nonmem
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
CSV.write("results/nonmem/nonmem.csv", nonmem_df)

# Pumas
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
CSV.write("results/pumas/pumas.csv", pumas_df)

# Results Overall
@rtransform! stan_df :software = "Stan"
@rtransform! nonmem_df :software = "NONMEM"
@rtransform! pumas_df :software = "Pumas"
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
    facet=(; linkyaxes=:none),
    palettes=(; color=colors)
)

save("results/results.png", fig; px_per_unit=3)

# density plots
model_specs_01 = Dict(
    :name => "01-md",
    :TVCL => 4,
    :TVVC => 70,
    :TVQ => 4,
    :TVVP => 50,
    :ωCL => 0.3,
    :ωVC => 0.3,
    :ωQ => 0.3,
    :ωVP => 0.3,
    :σ => 0.2,
)
model_specs_02 = Dict(
    :name => "02-md",
    :TVCL => 4,
    :TVVC => 70,
    :TVKA => 1,
    :ωCL => 0.3,
    :ωVC => 0.3,
    :ωKA => 0.3,
    :σ => 0.2,
)
model_specs_03 = Dict(
    :name => "03-md",
    :TVCL => 4,
    :TVVC => 70,
    :TVQ => 4,
    :TVVP => 50,
    :TVKA => 1,
    :ωCL => 0.3,
    :ωVC => 0.3,
    :ωKA => 0.3,
    :ωQ => 0.3,
    :ωVP => 0.3,
    :σ => 0.2,
)
model_specs_05 = Dict(
    :name => "05",
    :TVCL => 4,
    :TVVC => 70,
    :TVQ => 4,
    :TVVP => 50,
    :TVKA => 1,
    :TVMTT => 125,
    :TVCIRC0 => 5,
    :TVGAMMA => 0.17,
    :TVALPHA => 3e-4,
    :ωCL => 0.3,
    :ωVC => 0.3,
    :ωKA => 0.3,
    :ωQ => 0.3,
    :ωVP => 0.3,
    :ωMTT => 0.3,
    :ωCIRC0 => 0.3,
    :ωGAMMA => 0.3,
    :ωALPHA => 0.3,
    :σ => 0.2,
    :σ_pd => 0.2,
)
function _subset_chains(chn, param::Symbol, software::AbstractString)
    df = DataFrame(chn)
    @rselect! df $(param)
    @rtransform! df :software = software
    return df
end
function _combine_chains(chns, param)
    softwares = ["Stan", "NONMEM", "Pumas"]
    df = mapreduce((c, s) -> _subset_chains(c, param, s), vcat, chns, softwares)
    return df
end
function plot_posterior(chns, param::Symbol, model_specs)
    df = _combine_chains(chns, param)
    plt_data = data(df)
    plt_dens = mapping(param; color=:software) * AlgebraOfGraphics.density()# * visual(; alpha=0.5)
    plt_true_val = mapping([model_specs[param]]) * visual(VLines; color=:black, linestyle = :dash, linewidth=3)
    fig = draw(
        plt_data * plt_dens + plt_true_val;
        axis=(;
            title=model_specs[:name],
            ylabel="PDF"
        ),
        palettes=(; color=colors),
        legend=(position=:top, titleposition=:left, framevisible=false, padding=0)
    )
    save(
        joinpath(pwd(), "results", "posterior_plots", "$(model_specs[:name])-$(string(param)).png"),
        fig;
        px_per_unit=3,
    )
    return fig
end

chains_01 = [
    stan_chn_01_md,
    nonmem_chn_01_md,
    pumas_chn_01_md,
]
chains_02 = [
    stan_chn_02_md,
    nonmem_chn_02_md,
    pumas_chn_02_md,
]
chains_03 = [
    stan_chn_03_md,
    nonmem_chn_03_md,
    pumas_chn_03_md,
]
chains_05 = [
    stan_chn_05,
    nonmem_chn_05,
    pumas_chn_05,
]
params_01 = [
    :TVCL,
    :TVVC,
    :TVQ,
    :TVVP,
    :ωCL,
    :ωVC,
    :ωQ,
    :ωVP,
    :σ,
]
params_02 = [
    :TVCL,
    :TVVC,
    :TVKA,
    :ωCL,
    :ωVC,
    :ωKA,
    :σ,
]
params_03 = [
    :TVCL,
    :TVVC,
    :TVQ,
    :TVVP,
    :TVKA,
    :ωCL,
    :ωVC,
    :ωQ,
    :ωVP,
    :ωKA,
    :σ,
]
params_05 = [
    :TVCL,
    :TVVC,
    :TVQ,
    :TVVP,
    :TVKA,
    :TVMTT,
    :TVCIRC0,
    :TVGAMMA,
    :TVALPHA,
    :ωCL,
    :ωVC,
    :ωQ,
    :ωVP,
    :ωKA,
    :ωMTT,
    :ωCIRC0,
    :ωGAMMA,
    :ωALPHA,
    :σ,
    :σ_pd,
]
plot_posterior_01 = map(p -> plot_posterior(chains_01, p, model_specs_01), params_01)
plot_posterior_02 = map(p -> plot_posterior(chains_02, p, model_specs_02), params_02)
plot_posterior_03 = map(p -> plot_posterior(chains_03, p, model_specs_03), params_03)
plot_posterior_05 = map(p -> plot_posterior(chains_05, p, model_specs_05), params_05)

# QR Code to the repo (ask vijay for approval)