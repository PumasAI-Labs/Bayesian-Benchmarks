using Pumas

pk2cpt = @model begin
    @param begin
    tvcl ~ LogNormal(log(10), 0.25)
    tvvc ~ LogNormal(log(35), 0.25)
    tvq ~ LogNormal(log(15), 0.5)
    tvvp ~ LogNormal(log(105), 0.5)
    tvka ~ LogNormal(log(2.5), 1)
    # Cauchy is TDist with ν=1
    σ ~ truncated(TDist(1); lower=0)
    end
end
