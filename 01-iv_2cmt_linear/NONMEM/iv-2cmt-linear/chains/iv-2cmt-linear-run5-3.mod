;Model Desc: 01-iv-2cmt-linear
;Project Name: Bayesian Benchmarks

$PROB 
$INPUT ID	AMT	RATE	II	ADDL	CMT	EVID	SS	LLOQ	BLOQ	MDV	TIME	DV
$DATA ../single_dose.csv IGNORE=@

$SUBROUTINES ADVAN3 TRANS4

$PK 
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)

CL=EXP(MU_1+ETA(1))
V1=EXP(MU_2+ETA(2))
Q=EXP(MU_3+ETA(3))
V2=EXP(MU_4+ETA(4))
S1=V1


$ERROR
IPRED=A(1)/V1
Y = IPRED*(1+EPS(1))

; Initial Estimates  
$THETA
1.662 4.068 1.379 3.721

$OMEGA BLOCK(4)
0.248 0.01 0.128 0.01 0.01 0.452 0.01 0.01 0.01 0.243

$SIGMA 
0.216

; Prior Specification
$PRIOR NWPRI

$THETAP 1.386 FIX 4.248 FIX 1.386 FIX 3.912 FIX 
$THETAPV BLOCK(4) 
1 FIX
0 1
0 0 1
0 0 0 1

$OMEGAP BLOCK(4) 
0.16  FIX
0 0.16
0 0 0.16
0 0 0 0.16
$OMEGAPD 4 FIX

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
FILE=iv-2cmt-linear-iter-run5-3.tab

$COV UNCONDITIONAL MATRIX=R PRINT=E

$TABLE ID EVID TIME IPRED NOPRINT ONEHEADER FORMAT=,1PE13.6 FILE=iv-2cmt-linear-sdtab-run5-3.tab
$TABLE ID CL V1 Q V2 ETA1 ETA2 ETA3 ETA4 NOPRINT NOAPPEND ONEHEADER FORMAT=,1PE13.6 FILE=iv-2cmt-linear-patab-run5-3.tab