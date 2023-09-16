source("utils/nonmem2arrow.R")
library(furrr)

single = "03-depot_2cmt_linear/NONMEM/depot-2cmt-linear/chains"
multi = "03-depot_2cmt_linear/NONMEM/depot-2cmt-linear-md/chains"

future_walk(
  c(
    single,
    multi
  ),
  ~ for (r in 1:5) {
    convert_nonmem_to_arrow(.x, run = r)
  } 
)
