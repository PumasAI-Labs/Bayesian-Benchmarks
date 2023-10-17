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

# Plots
plt_rhat = data(df) *
      mapping(
          :parameters,
          :rhat;
          color=:software,
          # dodge=:software,
          row=:model
      ) *
      visual(Scatter; alpha=0.3)

fig_rhat = draw(
    plt_rhat;
    axis=(;
        xticklabelrotation=π / 3,
        ylabel=""
    ),
    facet=(; linkyaxes=:none)
)
save(joinpath(pwd(), "results", "rhat.png"), fig_rhat; px_per_unit=3)
    
plt_ess = data(df) *
      mapping(
          :parameters,
          :rhat;
          color=:software,
          # dodge=:software,
          row=:model
      ) *
      visual(Scatter; alpha=0.3)

fig_ess = draw(
    plt_ess;
    axis=(;
        xticklabelrotation=π / 3,
        ylabel=""
    ),
    facet=(; linkyaxes=:none)
)
save(joinpath(pwd(), "results", "ess.png"), fig_ess; px_per_unit=3)