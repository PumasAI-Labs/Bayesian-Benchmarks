source("utils/nonmem2arrow.R")
library(furrr)

single = "01-iv_2cmt_linear/NONMEM/iv-2cmt-linear/chains"
multi = "01-iv_2cmt_linear/NONMEM/iv-2cmt-linear-md/chains"

future_walk(
  c(
    single,
    multi
  ),
  convert_nonmem_to_arrow
)
