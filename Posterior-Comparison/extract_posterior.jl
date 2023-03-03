using Arrow
using MCMCChains
using DataFrames

function get_chains(
    df;
    nchains=4,
    internals=Symbol.(
        ["lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"]
    )
)
    time = first(df.time)
    df = select(df, Not(:time))
    mat = Array{Float64}(undef, Int(nrow(df)/nchains), ncol(df), nchains)
    for c in 1:nchains
        idx_min = Int(1 + (1_000 * (c - 1)))
        idx_max = Int(1_000 * c)
        # @show idx_min
        # @show size(df[idx_min:idx_max, :])
        mat[:, :, c ] .= df[idx_min:idx_max, :]
    end
    # info = (; start_time=fill(DateTime(today()), nchains), stop_time=fill(DateTime(today()) + Second(Int(floor(time))), nchains))
    info = (; start_time=fill(0.0, nchains), stop_time=fill(time, nchains))
    chn = Chains(mat, names(df), Dict(:internals => internals); info)
    # chn = setinfo(chn, info)
    return chn
end

files = filter(f -> endswith(f, ".arrow"), readdir(joinpath(pwd(), "Model1", "Stan", "Torsten", "Fits"); join=true))
df = mapreduce(f -> DataFrame(Arrow.Table(f)), vcat, files)

chn = get_chains(df)