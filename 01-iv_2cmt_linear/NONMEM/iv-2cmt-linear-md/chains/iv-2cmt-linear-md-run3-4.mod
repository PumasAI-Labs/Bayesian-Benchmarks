;Model Desc: 01-iv-2cmt-linear multiple dose
;Project Name: Bayesian Benchmarks

$PROB 
$INPUT ID	AMT	RATE	II	ADDL	CMT	EVID	SS	LLOQ	BLOQ	MDV	TIME	ODV DV
$DATA ../../../data/multiple_dose.csv IGNORE=@

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
IPRED=LOG(F)
Y=IPRED + EPS(1)

; Initial Estimates  
$THETA
2.048 4.213 1.485 3.921

$OMEGA BLOCK(4)
0.38 0.01 0.228 0.01 0.01 0.247 0.01 0.01 0.01 0.312

$SIGMA 
0.129

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
FILE=iv-2cmt-linear-md-iter-run3-4.tab

$COV UNCONDITIONAL MATRIX=R PRINT=E

$TABLE ID EVID TIME IPRED NOPRINT ONEHEADER FORMAT=,1PE13.6 FILE=iv-2cmt-linear-md-sdtab-run3-4.tab
$TABLE ID CL V1 Q V2 ETA1 ETA2 ETA3 ETA4 NOPRINT NOAPPEND ONEHEADER FORMAT=,1PE13.6 FILE=iv-2cmt-linear-md-patab-run3-4.tab
