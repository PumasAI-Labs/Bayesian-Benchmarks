# ACoP 13 Bayesian Workshop

## Benchmarks

### `Torsten` Setup

- `Torsten` version 0.90.0, which comes with `CmdStan` 2.29.2.
- `OpenMPI` version 4.1.2.

Create a local file in cmdstan/make folder with:

```bash
STAN_THREADS=true
TORSTEN_MPI=1
```

Then rebuild `CmdStan` in `cmdstanr`:

```R
library(cmdstanr)
set_cmdstan_path("Torsten/cmdstan") # Torsten CmdStan path
rebuild_cmdstan()
```
