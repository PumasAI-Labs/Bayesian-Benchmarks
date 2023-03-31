$PROB 1 cmt iv lin
$DATA ../../data/sd_iv_1cmt_lin.csv  IGNORE=@
$INPUT ID TIME CP DV AMT EVID CMT RATE DOSEGRP
$SUB ADVAN6 TOL=8
$MODEL COMP = (CENTR)

$PK
TVCL = THETA(1)
TVV = THETA(2)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
S1 = V/1000

$DES
DADT(1) = -CL*A(1)/V

$ERROR
IPRED=F
Y = IPRED*(1 + ERR(1))
$THETA
(0,4)
(0,70)

$OMEGA
0.09
0.09

$SIGMA
0.04

$EST METHOD=1 INTER MAXEVAL=9999

$COV PRINT=E