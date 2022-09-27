using Distributions
using JSON3

json_string = read(joinpath(@__DIR__, "data", "iv_2cmt_macro_ir1_cl_vc_q_vp_kin_kout_ic50_ppa_lognormal_001_bloq.json"), String)

stan_data = JSON3.read(json_string, NamedTuple)

stan_init = (;
    TVCL=rand(LogNormal(log(0.3), 0.3)),
    TVVC=rand(LogNormal(log(5), 0.3)),
    TVQ=rand(LogNormal(log(1), 0.3)),
    TVVP=rand(LogNormal(log(10), 0.3)),
    omega=rand(LogNormal(log(0.3), 0.3), 4),
    sigma=rand(LogNormal(log(0.5), 0.3), 2),
    TVKIN=rand(LogNormal(log(110), 0.3)),
    TVKOUT=rand(LogNormal(log(1), 0.3)),
    TVIC50=rand(LogNormal(log(6), 0.3)),
    omega_pd=rand(LogNormal(log(0.3), 0.3), 3),
    sigma_pd=rand(LogNormal(log(0.1), 0.3))
)
