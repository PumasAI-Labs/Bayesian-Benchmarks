source("utils/nonmem2arrow.R")
library(furrr)

single = "02-depot_1cmt_linear/NONMEM/depot-1cmt-linear/chains"
multi = "02-depot_1cmt_linear/NONMEM/depot-1cmt-linear-md/chains"

future_walk(
  c(
    single,
    multi
  ),
  ~ for (r in 1:5) {
    convert_nonmem_to_arrow(.x, run = r)
  } 
)
