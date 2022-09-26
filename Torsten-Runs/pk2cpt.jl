ENV["CMDSTAN"] = joinpath(@__DIR__, "..", "Torsten", "cmdstan")
using Stan
using StanSample

# setting up CMDSTAN path to Torsten
Stan.set_cmdstan_home!(joinpath(@__DIR__, "..", "Torsten", "cmdstan"))

m_str = read(joinpath(@__DIR__, "..", "Torsten", "example-models", "pk2cpt", "pk2cpt.stan"), String)
m = SampleModel("pk2cpt", m_str)

# data and inits
include(joinpath(@__DIR__, "pk2cpt_data.jl"))

rc = stan_sample(
    m;
    data=stan_data,
    init=stan_init,
    num_chains=4,
    use_cpp_chains=true,
    check_num_chains=false,
    num_cpp_chains=4,
    num_threads=Threads.nthreads(),
    num_samples=1_000,
    num_warmups=1_000,
    delta=0.8
)

if success(rc)
    summary_df = read_summary(m, false)
end