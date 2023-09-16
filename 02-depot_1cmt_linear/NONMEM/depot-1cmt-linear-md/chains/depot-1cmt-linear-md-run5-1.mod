;Model Desc: 02-depot-1cmt-linear-md
;Project Name: Bayesian Benchmarks

$PROB 
$INPUT ID	AMT	RATE	II	ADDL	CMT	EVID	SS	LLOQ	BLOQ	MDV	TIME	ODV DV
$DATA ../../../data/multiple_dose.csv IGNORE=@

$SUBROUTINES ADVAN2 TRANS2

$PK 
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)

CL=EXP(MU_1+ETA(1))
V=EXP(MU_2+ETA(2))
KA=EXP(MU_3+ETA(3))

K=CL/V
S2=V

$ERROR
IPRED=LOG(F)
Y=IPRED + EPS(1)

; Initial Estimates  
$THETA
1.684 3.926 -0.039

$OMEGA BLOCK(3)
0.45 0.01 0.259 0.01 0.01 0.418

$SIGMA 
0.289

; Prior Specification
$PRIOR NWPRI

$THETAP 1.386 FIX 4.248 FIX 0 FIX 
$THETAPV BLOCK(3) 
1 FIX
0 1
0 0 1

$OMEGAP BLOCK(3) 
0.16  FIX
0 0.16
0 0 0.16
$OMEGAPD 3 FIX

$SIGMAP 0.25 FIX
$SIGMAPD 1 FIX

$EST METHOD=NUTS 
INTERACTION
AUTO=1
PRINT=50 
NBURN=500 
NITER=1000 
OLKJDF=2 
SEED=12345
FILE=depot-1cmt-linear-md-iter-run5-1.tab
NUTS_TEST=1
NUTS_MASS=D

$COV UNCONDITIONAL MATRIX=R PRINT=E

$TABLE ID EVID TIME IPRED NOPRINT ONEHEADER FORMAT=,1PE13.6 FILE=depot-1cmt-linear-md-sdtab-run5-1.tab
$TABLE ID CL V K KA ETA1 ETA2 ETA3 NOPRINT NOAPPEND ONEHEADER FORMAT=,1PE13.6 FILE=depot-1cmt-linear-md-patab-run5-1.tab
