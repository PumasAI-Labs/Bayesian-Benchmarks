ENV["CMDSTAN"] = joinpath(@__DIR__, "..", "Torsten", "cmdstan")
using Stan
using StanSample
using DataFramesMeta

# setting up CMDSTAN path to Torsten
Stan.set_cmdstan_home!(joinpath(@__DIR__, "..", "Torsten", "cmdstan"))

m_str = read(joinpath(@__DIR__, "..", "Torsten", "example-models", "iv_2cmt_macro_ir1", "iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_bloq_lower_coupled.stan"), String)
m = SampleModel("iv_2cmt_macro_ir1", m_str)

# data and inits
include(joinpath(@__DIR__, "iv_2cmt_macro_ir1_data.jl"))

rc = stan_sample(
    m;
    use_cpp_chains=true,
    data=stan_data,
    init=stan_init,
    num_chains=4,
    num_samples=500,
    num_warmups=100,
    delta=0.98,
    max_depth=15,
    refresh=5
)

if success(rc)
    summary_df = read_summary(m, false)
end

parameters_to_summarize = [:TVCL, :TVVC, :TVQ, :TVVP, :TVKIN, :TVKOUT, :TVIC50]

@rsubset summary_df :parameters ∈ parameters_to_summarize