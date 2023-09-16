;Model Desc: 05-friberg
;Project Name: Bayesian Benchmarks

$SIZES LVR=50
$PROB
$INPUT ID AMT RATE II ADDL CMT EVID SS LLOQ BLOQ MDV TIME CP ODV DV
$DATA ../../../data/multiple_dose_nonmem.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=6

$MODEL
NCOMP=8
COMP = (GUT, DEFDOSE)
COMP = (CENT)
COMP = (PERIPH)
COMP = (PROL)
COMP = (TRANS1)
COMP = (TRANS2)
COMP = (TRANS3)
COMP = (CIRC)

$PK
MU_1 = THETA(1)     ; CL
MU_2 = THETA(2)     ; V1
MU_3 = THETA(3)     ; Q
MU_4 = THETA(4)     ; V2
MU_5 = THETA(5)     ; KA
MU_6 = THETA(6)     ; MTT
MU_7 = THETA(7)     ; CIRC0
MU_8 = THETA(8)     ; GAMMA
MU_9 = THETA(9)     ; ALPHA 

CL  = EXP(MU_1 +ETA(1))
V1  = EXP(MU_2 +ETA(2))
Q   = EXP(MU_3 +ETA(3))
V2  = EXP(MU_4 +ETA(4))
KA  = EXP(MU_5 +ETA(5))
MTT = EXP(MU_6 +ETA(6))
CIRC0 = EXP(MU_7 +ETA(7))
GAMMA = EXP(MU_8 +ETA(8))
ALPHA = EXP(MU_9 +ETA(9))

; intermediate calculations
K10 = CL/V1
K12 = Q/V1
K21 = Q/V2
KTR = 4/MTT

KPROL = KTR
KCIRC = KTR

; initial PK conditions
A_0(1) = 0
A_0(2) = 0
A_0(3) = 0

; initial PD conditions
A_0(4) = CIRC0
A_0(5) = CIRC0
A_0(6) = CIRC0
A_0(7) = CIRC0
A_0(8) = CIRC0

S2 = V1
S3 = V2


$DES
CONC = A(2)/V1 
EDRUG = ALPHA * CONC
IF(EDRUG > 1.0) EDRUG = 1.0

DADT(1) = -KA * A(1) 
DADT(2) = KA * A(1) - (K10 + K12) * A(2) + K21 * A(3)
DADT(3) = K12 * A(2) - K21 * A(3) 
DADT(4) = KPROL * A(4) * (1 - EDRUG) * ((CIRC0 / A(8))**GAMMA) - KTR * A(4)
DADT(5) = KTR * (A(4) - A(5))
DADT(6) = KTR * (A(5) - A(6))
DADT(7) = KTR * (A(6) - A(7))
DADT(8) = KTR * A(7) - KCIRC * A(8)


$ERROR 
IPRED = LOG(F)
IND = 0
IF(CMT.EQ.2) IND = 1 
YPK = IPRED + EPS(1) ; Exp Err PK
YPD = IPRED + EPS(2) ; Exp Err PD
Y = YPK*IND + YPD*(1-IND)

; Initial Estimates
$THETA
1.322 3.758 1.724 4.126 -0.215 4.893 1.519 -1.573 -8.335

$OMEGA BLOCK(9) 
0.474 0.01 0.34349 0.01 0.01 0.29704 0.01 0.01 0.01 0.44261 0.01 0.01 0.01 0.01 0.35863 0.01 0.01 0.01 0.01 0.01 0.35433 0.01 0.01 0.01 0.01 0.01 0.01 0.28806 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.42428 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.31372

$SIGMA 
0.33528 0.1529

; Prior Specification
$PRIOR NWPRI

$THETAP       
(1.386 FIX)     ; CL Log(4)  
(4.248 FIX)     ; V1 Log(70)
(1.386 FIX)     ; Q  Log(4)
(3.689 FIX)     ; V2 Log(40)
(0 FIX)         ; KA log(1)
(4.828 FIX)     ; MTT Log(125)
(1.609 FIX)     ; CIRC0 Log(5) 
(-1.772 FIX)    ; GAMMA Log(0.17)
(-8.112 FIX)    ; ALPHA Log(3E-4)

$THETAPV BLOCK(9)   
1   FIX   ; CL 
0.00   1      ; V1
0.00   0.00   1     ; Q
0.00   0.00   0.00   1    ; V2
0.00   0.00   0.00   0.00   1   ; KA 
0.00   0.00   0.00   0.00   0.00   1    ; MTT 
0.00   0.00   0.00   0.00   0.00   0.00   1   ; CIRC0 
0.00   0.00   0.00   0.00   0.00   0.00   0.00   1    ; GAMMA
0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  1    ; ALPHA

$OMEGAP BLOCK(9) 
0.16    FIX   ; CL 
0.00   0.16       ; V1
0.00   0.00   0.16    ; Q
0.00   0.00   0.00   0.16   ; V2
0.00   0.00   0.00   0.00   0.16  ; KA
0.00   0.00   0.00   0.00   0.00   0.16  ; MTT
0.00   0.00   0.00   0.00   0.00   0.00   0.16  ; CIRC0
0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.16  ; GAMMA
0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.16  ; ALPHA
$OMEGAPD (9, FIXED)

$SIGMAP 0.25 FIX 0.25 FIX ; INVERSE WISHART VAR
$SIGMAPD 1 1 FIX

$EST METHOD=NUTS 
INTERACTION
AUTO=1
PRINT=50 
NBURN=500 
NITER=1000 
OLKJDF=2 
SEED=12345
FILE=friberg-iter-run3-4.tab

$COV UNCONDITIONAL MATRIX=R PRINT=E

$TABLE ID EVID TIME IPRED NOPRINT ONEHEADER FORMAT=,1PE13.6 FILE=friberg-sdtab-run3-4.tab
$TABLE ID CL V1 Q V2 KA MTT CIRC0 GAMMA ALPHA ETAS(1:LAST) NOPRINT NOAPPEND ONEHEADER FORMAT=,1PE13.6 FILE=friberg-patab-run3-4.tab
