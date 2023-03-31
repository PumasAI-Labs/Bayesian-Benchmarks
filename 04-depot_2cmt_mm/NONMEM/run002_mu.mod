;; 1. Based on: run002
;; 2. Description: 002, model with parallel linear and Michelis-Menten elimination (concentration only)
;; x1. Author: user
$PROBLEM 002_mu, model with parallel linear and Michelis-Menten elimination (concentration only)

$INPUT ID,TIME,AMT,DROP,DROP,DV,EVID,MDV,CMT,DOSE=DROP,STUD=DROP,
       CFR1=DROP,CTO1=DROP,RFR1=DROP,RTO1=DROP,COM1=DROP

$DATA ../../data/tmdd/SimulatedNonmemDataConc.csv IGNORE = @

$SUBROUTINES ADVAN13 TOL=9

$MODEL
 COMP(Depot)
 COMP(Central,DEFDOSE,DEFOBS)
 COMP(Periph)

$PK
; Note: Automatic Mu-referencing by Pirana is currently only implemented for
;       simple linear and exponential parameter equations (par+eta or par*exp(eta)).
;       Please convert more complex equations manually.
  MU_1 = LOG(THETA(1))  ; ** MU-referenced by Pirana
  CL = EXP(MU_1+ETA(1)) ;    Original equation: CL = THETA(1)*EXP(ETA(1)))
  MU_2 = LOG(THETA(2))  ; ** MU-referenced by Pirana
  VC = EXP(MU_2+ETA(2)) ;    Original equation: VC = THETA(2)*EXP(ETA(2)))
  MU_3 = LOG(THETA(3)) ; ** MU-referenced by Pirana
  Q = EXP(MU_3+ETA(3)) ;    Original equation: Q  = THETA(3)*EXP(ETA(3)))
  VP = THETA(4)
  K = CL/VC
  K12 = Q/VC
  K21 = Q/VP

  F1 = THETA(5)
  KA = THETA(6)

  VMAX=THETA(7)
  MU_4 = LOG(THETA(8))   ; ** MU-referenced by Pirana
  KSS = EXP(MU_4+ETA(4)) ;    Original equation: KSS =THETA(8)*EXP(ETA(4)))

$DES
  DADT(1) =-KA*A(1)                                                       ; Free Drug depot
  DADT(2) = KA*A(1)-(K+K12)*A(2)+K21*A(3)-VMAX*A(2)/(KSS+A(2)/VC)         ; Free drug amount
  DADT(3) = K12*A(2)-K21*A(3)                                             ; Free Drug second compartment amount

$ERROR
  CALLFL = 0
  Y =  LOG(A(2)/VC)+EPS(1)

  IPRED = Y
  TY = 0
  IF(IPRED.GT.0) TY = EXP(IPRED)


$THETA
 (0,0.15)    ;1 CL
 (0,3.00)    ;2 V1
 (0,0.45)    ;3 Q
 (0,1.50)    ;4 V2
 (0,0.6)     ;5 F1
 (0,1.00)    ;6 KA
 (0,1.78)    ;7 VMAX
 (0,17)      ;8 KSS

$OMEGA
 0.07     ;1 CL
 0.04     ;2 V1
 0.03     ;3 Q
 0.13     ;5 KSS

$SIGMA
0.0268

$EST METHOD=SAEM
;$EST MAXEVAL=9999 NSIG=3 SIGL=9 METHOD=1 PRINT=10 NOABORT MSFO=002_mu.MSF
;NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST
;$COV PRINT=E UNCONDITIONAL MATRIX=S SIGL=12
