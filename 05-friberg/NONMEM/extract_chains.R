source("utils/nonmem2arrow.R")

dir = "05-friberg/NONMEM/friberg/chains"

for (r in 1:5) {
  convert_nonmem_to_arrow(dir, run = r)
} 
