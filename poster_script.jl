using Pumas, PumasPlots

model = @model begin
  @param begin
    θ ~ Constrained(
      MvNormal([1.9, 0.0781, 0.0463, 1.5, 0.4], Diagonal([9025, 15.25, 5.36, 5625, 400])),
      lower = [0.1, 0.008, 0.0004, 0.1, 0.0001],
      upper = [5, 0.5, 0.09, 5, 1.5],
      init = [1.9, 0.0781, 0.0463, 1.5, 0.4],
    )
    Ω ~ InverseWishart(2, fill(0.9, 1, 1) .* (2 + 1 + 1)) # NONMEM specifies the inverse Wishart in terms of its mode
    σ² ~ Gamma(1.0, 0.388)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @pre begin
    Ka = (SEX == 1 ? θ[1] : θ[4]) + η[1]
    K = θ[2]
    CL = θ[3] * (WT / 70)^θ[5]
    Vc = CL / K
    SC = Vc / (WT / 70)
  end

  @covariates SEX WT

  @vars begin
    conc = Central / SC
  end

  @dynamics Depots1Central1

  @derived begin
    dv ~ @. Normal(conc, sqrt(σ²))
  end
end

pop = read_pumas(
  joinpath(@__DIR__, "../../data/event_data/THEOPP.csv"),
  covariates = [:WT, :SEX],
)

# MCMC using NUTS

alg = Pumas.BayesMCMC()
ft = fit(model, theopp, Pumas.init_params(model), alg)

# Burn-in and thinning

tft = Pumas.truncate(ft, burnin = 1000, ratio = 0.2)

# Convergence diagnostics and statistics

Pumas.gewekediag(tft)
Pumas.heideldiag(tft)
Pumas.rafterydiag(tft)
Pumas.ess(tft)
Pumas.rhat(tft)
Pumas.mean(b)
Pumas.std(b)
Pumas.mcse(b)

# Posterior predictive checks

sims = simobs(tft, subject = 1, samples = 10)
sim_plot(sims)

# Cross-validation

cv_method = Pumas.PSISCrossvalidation(
	Pumas.KFold(K = 5),
	Pumas.BySubject(marginal = Pumas.LLQuad(imaxiters = 100)),
)
cv_result = Pumas.crossvalidate(tft, cv_method);
Pumas.elpd(cv_result)
