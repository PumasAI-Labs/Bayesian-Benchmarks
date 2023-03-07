# Bayesian Benchmarks

With Bayesian methods becoming more popular in the pharmacometrics community, a
thorough comparison of the available software is necessary to determine the 
accuracy and (speed) performance of the three main software packages that 
implement Bayesian methods for pharmacometrics - NONMEM, Pumas, and 
Stan/Torsten.

This repo contains code for all 3 of the above-mentioned software to implement
Bayesian estimation for some commonly used models:
  
+ A one-compartment model with IV infusion and first-order elimination 
+ A one-compartment model with first-order absorption and first-order 
elimination 
+ A two-compartment model with first-order absorption and first-order 
elimination 
+ A two-compartment model with first-order absorption and Michaelis-Menten 
elimination
+ A two-compartment PK model with IV infusion and a Friberg-Karlsson PD model

The code in this repository is meant to be as minimal as possible to compare 
performance. For a more thorough workflow that one might use in practice see
[here for NONMEM](https://github.com/metrumresearchgroup/iu-ctsi-2023-merge),
[here for Pumas](https://docs.pumas.ai/stable/basics/bayesian/), and
[here for Stan/Torsten](stanpmx.github.io).