source("utils/nonmem2arrow.R")
library(furrr)

single = "03-depot_2cmt_linear/NONMEM/depot-2cmt-linear/chains"
multi = "03-depot_2cmt_linear/NONMEM/depot-2cmt-linear-md/chains"

future_walk(
  c(
    single,
    multi
  ),
  convert_nonmem_to_arrow
)
