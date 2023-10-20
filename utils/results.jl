using CSV
using DataFrames
using DataFramesMeta
using CairoMakie
using AlgebraOfGraphics

stan_files = filter(f -> contains(f, r"\d.*.csv"), readdir(joinpath(pwd(), "results", "stan"); join=true))
nonmem_files = filter(f -> contains(f, r"\d.*.csv"), readdir(joinpath(pwd(), "results", "nonmem"); join=true))
pumas_files = filter(f -> contains(f, r"\d.*.csv"), readdir(joinpath(pwd(), "results", "pumas"); join=true))

function concat_dfs(files)
    df = @chain files begin
        map(
            _ -> CSV.read(_, DataFrame),
            files
        )
        vcat(_...; source=:model => ["01", "02", "03", "05"])
    end
    return df
end
stan_df = concat_dfs(stan_files)
nonmem_df = concat_dfs(nonmem_files)
pumas_df = concat_dfs(pumas_files)
@rtransform! stan_df :software = "Stan"
@rtransform! nonmem_df :software = "NONMEM"
@rtransform! pumas_df :software = "Pumas"
df = vcat(stan_df, nonmem_df, pumas_df)

ess_sec_df = @chain df begin
    groupby(:software)
    @combine :mean_ess_per_sec = mean(:ess_per_sec)
    @transform :comparison = :mean_ess_per_sec ./ minimum(:mean_ess_per_sec)
end

CSV.write(joinpath(pwd(), "results", "ess_per_sec.csv"), ess_sec_df)

# Plots
order_params = [
    "TVCL",
    "TVVC",
    "TVQ",
    "TVVP",
    "TVKA",
    "ωCL",
    "ωVC",
    "ωQ",
    "ωVP",
    "ωKA",
    "σ",
    "TVALPHA",
    "TVCIRC0",
    "TVGAMMA",
    "TVMTT",
    "ωALPHA",
    "ωCIRC0",
    "ωGAMMA",
    "ωMTT",
    "σ_pd",
]

colors = [
    colorant"#008579", # nonmem
    colorant"#0035c7", # pumas 
    colorant"#b2001d", # stan
]

plt_rhat = data(df) *
      mapping(
          :parameters => sorter(order_params) => "Parameters",
          :rhat => "Rhat";
          color=:software,
          dodge=:software,
          row=:model
      ) *
      visual(BarPlot)

fig_rhat = draw(
    plt_rhat;
    axis=(;
        xticklabelrotation=π / 3,
        limits = (nothing, nothing, 0.99, nothing)
    ),
    facet=(; linkyaxes=:none),
    palettes=(; color=colors),
    legend=(position=:top, titleposition=:left, framevisible=false, padding=0),
)
save(joinpath(pwd(), "results", "rhat.png"), fig_rhat; px_per_unit=3)

# TODO: fix y-axis 
plt_ess = data(df) *
      mapping(
          :parameters => sorter(order_params) => "Parameters",
          :ess_tail => "Effective Sample Size";
          color=:software,
          dodge=:software,
          row=:model
      ) *
      visual(BarPlot)

fig_ess = draw(
    plt_ess;
    axis=(;
        xticklabelrotation=π / 3,
    ),
    palettes=(; color=colors),
    legend=(position=:top, titleposition=:left, framevisible=false, padding=0),
)
save(joinpath(pwd(), "results", "ess_tail.png"), fig_ess; px_per_unit=3)

plt_ess_sec = data(df) *
      mapping(
          :parameters => sorter(order_params) => "Parameters",
          :ess_per_sec => "Effective Sample Size / second";
          color=:software,
          dodge=:software,
          row=:model
      ) *
      visual(BarPlot)

fig_ess_sec = draw(
    plt_ess_sec;
    axis=(;
        xticklabelrotation=π / 3,
    ),
    facet=(; linkyaxes=:none),
    palettes=(; color=colors),
    legend=(position=:top, titleposition=:left, framevisible=false, padding=0),
    # figure=(; resolution=(800, 1200))
)
save(joinpath(pwd(), "results", "ess_per_sec.png"), fig_ess_sec; px_per_unit=3)