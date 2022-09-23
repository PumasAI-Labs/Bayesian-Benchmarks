using Pumas
using DataFrames

pk2cpt_ode = @model begin
    @param begin
        tvcl ~ LogNormal(log(10), 0.25) # CL
        tvq ~ LogNormal(log(15), 0.5)   # Q
        tvvc ~ LogNormal(log(35), 0.25) # V1
        tvvp ~ LogNormal(log(105), 0.5) # V2
        tvka ~ LogNormal(log(2.5), 1)   # ka
        σ ~ truncated(Cauchy(), 0, Inf) # sigma
    end

    @pre begin
        CL = tvcl
        Vc = tvvc
        Q = tvq
        Vp = tvvp
        Ka = tvka
    end

    @dynamics begin
        Depot'      = -Ka*Depot
        Central'    =  Ka*Depot -(CL+Q)/Vc*Central + Q/Vp*Peripheral
        Peripheral' =                  Q/Vc*Central - Q/Vp*Peripheral
    end

    # Don't linearize, we want a numerical solver!
    @options checklinear=false

    @derived begin
        # Torsten uses log(dv) ~ Normal(log(cp), sigma)
        cp := @. Central/Vc
        dv ~ @. LogNormal(log(cp), σ)
    end
end

# Torsten Data
addl = [14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
amt = [80000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
cmt = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
cObs = [missing, 359.239725273613, 662.647763577673, 1106.23947626728,
1185.25862586227, 1802.42630338229, 2296.47896937277, 2008.04412122616,
2000.94001020581, 1115.28848387488, 902.769414322048, 445.989900210258,
285.638744806589, 333.509854464029, 664.522440159846, 960.044660816181,
1157.24318831826, 1593.58511079002, 2170.22768767046, 2175.56839105518,
2168.6215871781]
evid = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ii = [12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
iObs = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
nObs = 20
nt = 21
rate = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ss = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
time = [0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12,
12.083, 12.167, 12.25, 12.5, 12.75, 13, 13.5]

df = DataFrame(; time, id=1, ii, evid, dv=cObs, cmt, addl, amt)

# Pumas Population
pop = read_pumas(df)

# Inits from Torsten
# TODO: After fix bug order inits same as Torsten
init_params = (;
    tvcl=7.4367958406427,
    tvq=28.0799996152587,
    tvvc=78.4460632446725,
    tvvp=68.1255965629187,
    tvka=1.0811298754049,
    σ=0.589695154260051,
)

pk2cpt_ode_fit = fit(
    pk2cpt_ode,
    pop,
    init_params,
    Pumas.BayesMCMC(
        nsamples=2_000,
        nadapts=1_000,
        target_accept=0.8,
        nchains=4,
        ensemblealg = EnsembleThreads(),
        parallel_chains = true,
        # Same options as Torsten: pmx_solve_rk45(..., 1e-5, 1e-8, 1e5)
        diffeq_options = (;
            alg=Tsit5(), # 4th/5th order RK Solver
            reltol=1e-5,
            abstol=1e-8,
            maxiters=Int(1e5),
        )
    );
)

Pumas.truncate(pk2cpt_ode_fit; burnin=1_000)
