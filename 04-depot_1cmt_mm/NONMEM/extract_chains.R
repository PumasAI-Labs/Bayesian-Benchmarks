source("utils/nonmem2arrow.R")
library(furrr)

single = "04-depot_1cmt_mm/NONMEM/depot-1cmt-mm/chains"
multi = "04-depot_1cmt_mm/NONMEM/depot-1cmt-mm-md/chains"

future_walk(
  c(
    single,
    multi
  ),
  ~ for (r in 1:5) {
    convert_nonmem_to_arrow(.x, run = r)
  } 
)
