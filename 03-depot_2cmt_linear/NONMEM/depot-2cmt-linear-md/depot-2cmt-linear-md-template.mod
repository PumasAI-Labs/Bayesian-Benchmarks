;Model Desc: 03-depot-2cmt-linear-md
;Project Name: Bayesian Benchmarks

$PROB 
$INPUT ID	AMT	RATE	II	ADDL	CMT	EVID	SS	LLOQ	BLOQ	MDV	TIME	ODV DV
$DATA ../../../data/multiple_dose.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

$PK 
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
MU_5=THETA(5)

CL=EXP(MU_1+ETA(1))
V2=EXP(MU_2+ETA(2))
Q=EXP(MU_3+ETA(3))
V3=EXP(MU_4+ETA(4))
KA=EXP(MU_5+ETA(5))

K=CL/V2
K23=Q/V2
K32=Q/V3
S2=V2

$ERROR
IPRED=LOG(F)
Y=IPRED + EPS(1)

; Initial Estimates  
$THETA
THETAUPDATE

$OMEGA BLOCK(5)
OMEGAUPDATE

$SIGMA 
SIGMAUPDATE

; Prior Specification
$PRIOR NWPRI

$THETAP 1.386 FIX 4.248 FIX 1.386 FIX 3.689 FIX 0 FIX 
$THETAPV BLOCK(5) 
1 FIX
0 1
0 0 1
0 0 0 1
0 0 0 0 1

$OMEGAP BLOCK(5) 
0.16  FIX
0 0.16
0 0 0.16
0 0 0 0.16
0 0 0 0 0.16
$OMEGAPD 5 FIX

$SIGMAP 0.25 FIX
$SIGMAPD 1 FIX

$EST METHOD=NUTS 
INTERACTION
PRINT=50 
NBURN=500 
NITER=1000 
OLKJDF=2 
SEED=12345
FILE=ITERFILE

$COV UNCONDITIONAL MATRIX=R PRINT=E

$TABLE ID EVID TIME IPRED NOPRINT ONEHEADER FORMAT=,1PE13.6 FILE=SDFILE
$TABLE ID CL V2 Q V3 KA K K23 K32 ETAS(1:LAST) NOPRINT NOAPPEND ONEHEADER FORMAT=,1PE13.6 FILE=PAFILE
