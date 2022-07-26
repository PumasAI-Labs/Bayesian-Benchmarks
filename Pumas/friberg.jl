using Pumas
using DataFrames
using CSV

friberg = @model begin
    @param begin
        θ ~ MvNormal(
            [2.3, 2.7, 3.56, 4.7, 0.693, 4.83, 1.61, -8.11, -1.78, -1.60, -1.60],
            Diagonal([0.25, 0.25, 0.25, 0.25, 0.25, 0.04, 0.04, 1.0, 0.04, 1.0, 1.0]))
        Ω ~ InverseWishart(11, diagm(fill(0.045, 8)) .* (11 + 8 + 1))
    end

    @random begin
        ηstd ~ MvNormal(I(8))
    end

    @covariates vwt

    @pre begin
        η = cholesky(Ω).L * ηstd
        CL    = exp(θ[1] + 0.75*vwt + η[1]) # CL
        Q     = exp(θ[2] + 0.75*vwt + η[2]) # Vc
        Vc    = exp(θ[3] +      vwt + η[3]) # Q
        Vp    = exp(θ[4] +      vwt + η[4]) # Vp
        Ka    = exp(θ[5]            + η[5]) # Ka
        MTT   = exp(θ[6]            + η[6]) # CIRC0
        CIRC0 = exp(θ[7]            + η[7]) # MTT
        α     = exp(θ[8]            + η[8]) # α
        γ     = exp(θ[9])		    # γ
    end

    @init begin
        Prol     = CIRC0
        Transit1 = CIRC0
        Transit2 = CIRC0
        Transit3 = CIRC0
        Circ     = CIRC0
    end

    @vars begin
        KTR   = 4/MTT
        CONC  = Central/Vc
        EDRUG = min(α*CONC, 1)
    end

    @dynamics begin
        #PK kinetic model
        Depot'      = -Ka*Depot
        Central'    =  Ka*Depot - (CL/Vc + Q/Vc)*Central + (Q/Vp)*Peripheral
        Peripheral' =  Q/Vc*Central - Q/Vp*Peripheral
        # PD
        Prol'       = KTR*Prol*(1 - EDRUG)*abs(CIRC0/Circ)^γ - KTR*Prol
        Transit1'   = KTR*Prol     - KTR*Transit1
        Transit2'   = KTR*Transit1 - KTR*Transit2
        Transit3'   = KTR*Transit2 - KTR*Transit3
        Circ'       = KTR*Transit3 - KTR*Circ
    end

    @derived begin
        # dv ~ @. Normal(log(Central/Vc), exp(σpk))
        # E  ~ @. Normal(log(Circ)      , exp(σpd))
        dv ~ @. Normal(log(abs(Central/Vc)), exp(θ[10]))
        E  ~ @. Normal(log(abs(Circ))      , exp(θ[11]))
    end
end

# Read data
tmpdf = CSV.read("data/friberg.csv", DataFrame; missingstring=".")
df_obs = filter(t -> t.evid == 0, tmpdf)
df_obs_unstacked = unstack(df_obs, [:id, :time, :evid], :cmt, :LNDV; renamecols=x->x==2 ? "dv" : "E")
df_events = filter(t -> t.evid != 0, tmpdf)
df = outerjoin(df_obs_unstacked, df_events; on=[:id, :time, :evid])
df[:,:vwt] = log.(df.weight/70)
sort!(df, [:id, :time, :evid])

# Pumas Population
pop = read_pumas(df;
    observations=[:dv, :E],
    covariates=[:vwt]
)

# Inits from Pumas
init_params = (
    # tvcl        =  2.29293650954706,
    # tvq         =  2.78187002197878,
    # tvvc        =  4.10653043167607,
    # tvvp        =  5.79441188940434,
    # tvka        =  0.758743417852145,
    # tvPOP_CIRC0 =  1.7904012807499,
    # tvPOP_MTT   =  4.78975235013628,
    # tvα         = -8.82748705123275,
    # tvγ         = -1.91147638613364,
    # σpk         = -1.96694893289804,
    # σpd         = -1.53523045323129,
    θ = [
         2.29293650954706,
         2.78187002197878,
         4.10653043167607,
         5.79441188940434,
         0.758743417852145,
         4.78975235013628,
         1.7904012807499,
        -8.82748705123275,
        -1.91147638613364,
        -1.96694893289804,
        -1.53523045323129],
    Ω = Pumas.vechinv([
    0.384756084711769
    -0.00862682205633429
    0.133928083343693
    -0.00611937059477515
    -0.00611937059477515
    0.360525225413911
    -0.00896799230792607
    -0.00896799230792607
    -0.00896799230792607
    0.202699208207848
    -0.00016074572673854
    -0.00016074572673854
    -0.00016074572673854
    -0.00016074572673854
    0.206098155610402
    0.00451724596821987
    0.00451724596821987
    0.00451724596821987
    0.00451724596821987
    0.00451724596821987
    0.42147934456551
    -0.00331225454945625
    -0.00331225454945625
    -0.00331225454945625
    -0.00331225454945625
    -0.00331225454945625
    -0.00331225454945625
    0.158636749633111
    0.0168748417643117
    0.0168748417643117
    0.0168748417643117
    0.0168748417643117
    0.0168748417643117
    0.0168748417643117
    0.0168748417643117
    0.0723822402125386]))

# max_chunk_size is number of fixed effects + number of random effects
friberg_fit = fit(friberg,
                 pop,
                 init_params,
                 Pumas.BayesMCMC(nsamples=2, nadapts=1, target_accept=0.8, nchains=4))

Pumas.truncate(friberg_fit; burnin=1_000)
