using CSV
using DataFramesMeta

ext_file = readlines(joinpath(@__DIR__, "pk2cpt.ext"))
ext_file = strip.(ext_file)
ext_file = map(ext_file) do x 
    replace(x, [" "^i => " " for i in 5:-1:2]...)
end
ext_file = map(ext_file) do x 
    replace(x, [" "^i => " " for i in 5:-1:2]...)
end


coef_pk2cpt = CSV.read(
    IOBuffer(join(ext_file, "\n")), DataFrame;
    delim=" ",
    header=2,
    skipto=3
)

@rsubset! coef_pk2cpt :ITERATION > 0

exp_mean(x) = mean(exp.(x))

@combine coef_pk2cpt begin
    $(Not([:ITERATION, Symbol("SIGMA(1,1)")]) .=> exp_mean)
    $(Symbol("SIGMA(1,1)") => mean)
end 