ENV["CMDSTAN"] = joinpath(@__DIR__, "..", "Torsten", "cmdstan")
ENV["STAN_NUM_THREADS"] = Threads.nthreads()
using Stan
using StanSample
using DataFramesMeta

# setting up CMDSTAN path to Torsten
Stan.set_cmdstan_home!(joinpath(@__DIR__, "..", "Torsten", "cmdstan"))

m_str = read(joinpath(@__DIR__, "..", "Torsten", "example-models", "poppk2cpt", "depot_2cmt_match_metrum_half_normal_omega_ode.stan"), String)
m = SampleModel("poppk2cpt", m_str)

# data and inits
# 1 => linode (matrix exponential), 2 => general ode (rk45), 3 => general ode (bdf)
ode_solver = 1
include(joinpath(@__DIR__, "poppk2cpt-reduce_sum_ode_data.jl"))

rc = stan_sample(
    m;
    use_cpp_chains=true,
    data=stan_data,
    init=stan_init,
    num_cpp_chains=4,
    num_threads=Threads.nthreads(),
    num_samples=1_000,
    num_warmups=1_000,
    delta=0.8
)

# to see the logs while sampling do in a terminal:
# ;cat "$(joinpath(m.tmpdir, m.name*"_log_1.log"))"

if success(rc)
    summary_df = read_summary(m, false)
end

parameters_to_summarize = [:TVCL, :TVVC, :TVQ, :TVVP, :TVKA]

@rsubset summary_df :parameters ∈ parameters_to_summarize