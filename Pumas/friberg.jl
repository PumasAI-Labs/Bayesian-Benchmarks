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
        η ~ MvNormal(Ω)
    end

    @covariates vwt

    @pre begin
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
iparams = (
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
        -1.53523045323129,
    ],
    Ω = Pumas.vechinv([
        0.384756084711769,
        -0.00862682205633429,
        0.133928083343693,
        -0.00611937059477515,
        -0.00611937059477515,
        0.360525225413911,
        -0.00896799230792607,
        -0.00896799230792607,
        -0.00896799230792607,
        0.202699208207848,
        -0.00016074572673854,
        -0.00016074572673854,
        -0.00016074572673854,
        -0.00016074572673854,
        0.206098155610402,
        0.00451724596821987,
        0.00451724596821987,
        0.00451724596821987,
        0.00451724596821987,
        0.00451724596821987,
        0.42147934456551,
        -0.00331225454945625,
        -0.00331225454945625,
        -0.00331225454945625,
        -0.00331225454945625,
        -0.00331225454945625,
        -0.00331225454945625,
        0.158636749633111,
        0.0168748417643117,
        0.0168748417643117,
        0.0168748417643117,
        0.0168748417643117,
        0.0168748417643117,
        0.0168748417643117,
        0.0168748417643117,
        0.0723822402125386,
    ]; lower=false),
)

friberg_fit = fit(
    friberg,
    pop,
    iparams,
    Pumas.BayesMCMC(
        nsamples=2_000,
        nadapts=1_000,
        target_accept=0.8,
        nchains=4,
        ensemblealg = EnsembleThreads(),
        parallel_subjects=true,
        parallel_chains=true,
        diffeq_options=(; alg=Rodas5P())
    ),
)

Pumas.truncate(friberg_fit; burnin=1_000)

#Chains MCMC chain (1000×47×4 Array{Float64, 3}):
#
#Iterations        = 1:1:1000
#Number of chains  = 4
#Samples per chain = 1000
#Wall duration     = 62474.32 seconds
#Compute duration  = 238713.61 seconds
#parameters        = θ₁, θ₂, θ₃, θ₄, θ₅, θ₆, θ₇, θ₈, θ₉, θ₁₀, θ₁₁, Ω₁,₁, Ω₂,₁, Ω₃,₁, Ω₄,₁, Ω₅,₁, Ω₆,₁, Ω₇,₁, Ω₈,₁, Ω₂,₂, Ω₃,₂, Ω₄,₂, Ω₅,₂, Ω₆,₂, Ω₇,₂, Ω₈,₂, Ω₃,₃, Ω₄,₃, Ω₅,₃, Ω₆,₃, Ω₇,₃, Ω₈,₃, Ω₄,₄, Ω₅,₄, Ω₆,₄, Ω₇,₄, Ω₈,₄, Ω₅,₅, Ω₆,₅, Ω₇,₅, Ω₈,₅, Ω₆,₆, Ω₇,₆, Ω₈,₆, Ω₇,₇, Ω₈,₇, Ω₈,₈

# Row │ parameters  mean          std        naive_se     mcse         ess       rhat      ess_per_sec
#     │ Symbol      Float64       Float64    Float64      Float64      Float64   Float64   Float64
#─────┼────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ θ₁           2.40096      0.173637   0.00274545   0.00426498   1607.85   1.00231    0.00673549
#   2 │ θ₂           3.09503      0.187216   0.00296014   0.00425341   2139.29   0.99982    0.00896176
#   3 │ θ₃           3.37377      0.210057   0.00332129   0.0039951    2792.62   1.00019    0.0116986
#   4 │ θ₄           4.74644      0.18795    0.00297175   0.00365949   2347.18   1.00117    0.0098326
#   5 │ θ₅           0.522187     0.215115   0.00340127   0.00404207   2646.12   1.00059    0.0110849
#   6 │ θ₆           4.78032      0.129326   0.00204483   0.00243849   3318.96   1.00016    0.0139035
#   7 │ θ₇           1.72019      0.132774   0.00209934   0.00221375   3593.0    0.999634   0.0150515
#   8 │ θ₈          -8.05833      0.183918   0.00290799   0.0037257    2717.1    1.00074    0.0113823
#   9 │ θ₉          -1.8947       0.0808207  0.00127789   0.00118303   5802.66   1.00035    0.024308
#  10 │ θ₁₀         -2.29952      0.0459511  0.00072655   0.000573539  4845.59   0.999345   0.0202988
#  11 │ θ₁₁         -2.2417       0.0621105  0.000982054  0.000781926  6081.57   0.999621   0.0254764
#  12 │ Ω₁,₁         0.197526     0.129103   0.0020413    0.0025463    2318.92   1.00033    0.00971423
#  13 │ Ω₂,₁        -0.0293852    0.0914842  0.00144649   0.00188366   2323.75   0.999793   0.00973449
#  14 │ Ω₃,₁        -0.00106438   0.0890756  0.00140841   0.00203389   1768.53   1.00296    0.0074086
#  15 │ Ω₄,₁         0.0433499    0.0924096  0.00146112   0.00192493   2488.62   1.00001    0.0104251
#  16 │ Ω₅,₁         0.000521161  0.0905874  0.00143231   0.00186225   1962.36   1.0008     0.00822057
#  17 │ Ω₆,₁        -0.0192981    0.0784549  0.00124048   0.00137745   2678.43   1.00061    0.0112203
#  18 │ Ω₇,₁         0.0170693    0.0840183  0.00132845   0.00201927   2020.57   1.00132    0.00846442
#  19 │ Ω₈,₁         0.0400457    0.0878807  0.00138952   0.00178394   2247.82   1.00145    0.00941639
#  20 │ Ω₂,₂         0.22494      0.146996   0.00232421   0.00305321   2154.15   1.00199    0.009024
#  21 │ Ω₃,₂        -0.0092238    0.0906514  0.00143333   0.00185336   2086.51   0.999999   0.00874064
#  22 │ Ω₄,₂        -0.0585148    0.112048   0.00177164   0.00246371   1969.93   1.00069    0.00825227
#  23 │ Ω₅,₂         0.00225752   0.0920258  0.00145506   0.0014498    2576.37   1.00002    0.0107927
#  24 │ Ω₆,₂         0.0112175    0.0793366  0.00125442   0.00166719   2377.82   1.00008    0.00996099
#  25 │ Ω₇,₂         0.0304913    0.087295   0.00138026   0.00193373   1973.73   1.00141    0.00826819
#  26 │ Ω₈,₂        -0.0159139    0.0928395  0.00146792   0.00185986   2296.0    1.00101    0.00961823
#  27 │ Ω₃,₃         0.176635     0.124126   0.00196261   0.00308059   1378.01   1.00368    0.00577264
#  28 │ Ω₄,₃        -0.00454505   0.0966621  0.00152836   0.0022592    1761.28   0.999914   0.00737822
#  29 │ Ω₅,₃         0.0260362    0.106191   0.00167903   0.00291724   1111.52   1.0043     0.0046563
#  30 │ Ω₆,₃         0.0050754    0.0790164  0.00124936   0.00217706   1374.43   1.00238    0.00575764
#  31 │ Ω₇,₃        -0.00114014   0.0848096  0.00134096   0.00184111   2036.13   1.00085    0.00852958
#  32 │ Ω₈,₃         0.00410712   0.102657   0.00162315   0.00314508   1192.65   1.00332    0.00499616
#  33 │ Ω₄,₄         0.248914     0.170564   0.00269685   0.00382296   1688.9    1.00365    0.007075
#  34 │ Ω₅,₄        -0.0157457    0.10269    0.00162367   0.00231129   1745.31   1.00039    0.00731132
#  35 │ Ω₆,₄        -0.00512414   0.0824931  0.00130433   0.00204274   1792.03   1.00049    0.00750705
#  36 │ Ω₇,₄        -0.0190371    0.0912626  0.00144299   0.00201719   1915.61   0.999859   0.00802473
#  37 │ Ω₈,₄         0.0204887    0.101382   0.00160299   0.00243199   2019.69   1.00323    0.00846073
#  38 │ Ω₅,₅         0.185193     0.151543   0.0023961    0.00435455   1219.89   1.00311    0.00511025
#  39 │ Ω₆,₅        -0.00207355   0.090913   0.00143746   0.00272314   1196.0    1.00234    0.00501017
#  40 │ Ω₇,₅         0.0077686    0.0869125  0.00137421   0.00181442   2387.7    0.999772   0.0100024
#  41 │ Ω₈,₅        -0.00227091   0.117569   0.00185893   0.00389757    985.131  1.00342    0.00412683
#  42 │ Ω₆,₆         0.152863     0.102699   0.00162382   0.00259236   1628.64   1.0023     0.00682256
#  43 │ Ω₇,₆        -0.012734     0.071411   0.00112911   0.00169953   2129.98   1.00088    0.00892276
#  44 │ Ω₈,₆        -0.0119375    0.0896106  0.00141687   0.00214578   1458.56   1.00163    0.00611009
#  45 │ Ω₇,₇         0.180456     0.119795   0.00189413   0.00268585   1544.73   1.00507    0.00647107
#  46 │ Ω₈,₇         0.018523     0.089131   0.00140928   0.00189746   1988.35   1.00345    0.00832945
#  47 │ Ω₈,₈         0.191414     0.151054   0.00238837   0.00544954   1164.24   1.00825    0.00487712

# Row │ parameters  2.5%        25.0%        50.0%         75.0%        97.5%
#     │ Symbol      Float64     Float64      Float64       Float64      Float64
#─────┼───────────────────────────────────────────────────────────────────────────
#   1 │ θ₁           2.05354     2.2922       2.4061        2.51244      2.73261
#   2 │ θ₂           2.72008     2.97912      3.09806       3.22053      3.449
#   3 │ θ₃           2.95183     3.23707      3.37944       3.51391      3.77555
#   4 │ θ₄           4.37626     4.62438      4.74556       4.86543      5.12417
#   5 │ θ₅           0.101574    0.379525     0.524523      0.666229     0.935271
#   6 │ θ₆           4.52824     4.69342      4.77726       4.86776      5.03719
#   7 │ θ₇           1.44848     1.63294      1.72423       1.81052      1.97151
#   8 │ θ₈          -8.42254    -8.17647     -8.05503      -7.93817     -7.71041
#   9 │ θ₉          -2.05111    -1.94974     -1.89342      -1.83904     -1.7422
#  10 │ θ₁₀         -2.38736    -2.33075     -2.29986      -2.26933     -2.20693
#  11 │ θ₁₁         -2.3609     -2.28357     -2.24206      -2.20197     -2.1184
#  12 │ Ω₁,₁         0.0694626   0.116359     0.16182       0.234534     0.536496
#  13 │ Ω₂,₁        -0.224692   -0.0671492   -0.0232149     0.0125962    0.147643
#  14 │ Ω₃,₁        -0.179261   -0.0381205   -0.00226713    0.0346698    0.180444
#  15 │ Ω₄,₁        -0.128686   -0.00305199   0.0350477     0.0809429    0.252033
#  16 │ Ω₅,₁        -0.173075   -0.037154    -0.00150195    0.0367284    0.189951
#  17 │ Ω₆,₁        -0.1982     -0.0473779   -0.0144774     0.0163289    0.120668
#  18 │ Ω₇,₁        -0.128286   -0.0195358    0.0131698     0.0487396    0.196559
#  19 │ Ω₈,₁        -0.103156   -0.00148911   0.030373      0.0729931    0.239982
#  20 │ Ω₂,₂         0.0795511   0.13455      0.188714      0.26711      0.59212
#  21 │ Ω₃,₂        -0.197993   -0.0464516   -0.00557494    0.0317196    0.169826
#  22 │ Ω₄,₂        -0.301724   -0.096101    -0.0454214    -0.00330148   0.116743
#  23 │ Ω₅,₂        -0.182898   -0.0377751    0.000868498   0.0422264    0.185755
#  24 │ Ω₆,₂        -0.140363   -0.0242532    0.00956991    0.0451482    0.172218
#  25 │ Ω₇,₂        -0.120434   -0.0131696    0.0231987     0.0648365    0.229341
#  26 │ Ω₈,₂        -0.214463   -0.0543009   -0.0120991     0.025934     0.155268
#  27 │ Ω₃,₃         0.0594608   0.10223      0.143888      0.211865     0.492809
#  28 │ Ω₄,₃        -0.207727   -0.0439609   -0.000112946   0.0405854    0.179368
#  29 │ Ω₅,₃        -0.121216   -0.0162159    0.015906      0.055299     0.218061
#  30 │ Ω₆,₃        -0.134316   -0.0273969    0.00331008    0.0332805    0.153738
#  31 │ Ω₇,₃        -0.160008   -0.0361757   -0.00209941    0.0319583    0.161045
#  32 │ Ω₈,₃        -0.158964   -0.0348759    0.00105157    0.0366645    0.174705
#  33 │ Ω₄,₄         0.0858943   0.148495     0.204488      0.293388     0.667596
#  34 │ Ω₅,₄        -0.22269    -0.0562148   -0.0120068     0.0302786    0.171839
#  35 │ Ω₆,₄        -0.162828   -0.0401735   -0.0050903     0.0301398    0.148723
#  36 │ Ω₇,₄        -0.199257   -0.0575547   -0.0155096     0.0218819    0.148601
#  37 │ Ω₈,₄        -0.156232   -0.0217174    0.0180915     0.0602689    0.229197
#  38 │ Ω₅,₅         0.0615278   0.104688     0.1489        0.21814      0.510689
#  39 │ Ω₆,₅        -0.147116   -0.034007    -0.00226252    0.0285476    0.134939
#  40 │ Ω₇,₅        -0.156836   -0.0272202    0.00646034    0.0413076    0.175993
#  41 │ Ω₈,₅        -0.171125   -0.0416502   -0.00413351    0.0310456    0.168721
#  42 │ Ω₆,₆         0.0545498   0.0906867    0.127568      0.18141      0.39627
#  43 │ Ω₇,₆        -0.159942   -0.0420534   -0.00869688    0.0211975    0.114251
#  44 │ Ω₈,₆        -0.163331   -0.0432285   -0.0102393     0.0198051    0.121135
#  45 │ Ω₇,₇         0.0616341   0.107281     0.149137      0.214767     0.511465
#  46 │ Ω₈,₇        -0.133056   -0.0186529    0.0139481     0.0503333    0.19691
#  47 │ Ω₈,₈         0.0646225   0.111143     0.155932      0.225477     0.515419
