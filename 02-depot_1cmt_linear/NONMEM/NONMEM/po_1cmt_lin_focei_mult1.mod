$PROB 1 cmt po lin
$DATA ../../data/sd_po_1cmt_lin.csv  IGNORE=@
$INPUT ID TIME CP DV AMT EVID CMT RATE DOSEGRP
$SUB ADVAN6 TOL=8
$MODEL COMP=(DEPOT) COMP=(CENTRAL, DEFOBS)

$PK
TVCL = THETA(1)
TVV = THETA(2)
TVKA = THETA(3)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
KA = TVKA * EXP(ETA(3))
S2 = V/1000

$DES
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - (CL/V)*A(2)

$ERROR
IPRED=F
Y = IPRED*(1 + ERR(1))
$THETA
(0,4)
(0,70)
(0,1)

$OMEGA
0.09
0.09
0.09

$SIGMA
0.04

$EST METHOD=1 INTER MAXEVAL=9999

$COV PRINT=E