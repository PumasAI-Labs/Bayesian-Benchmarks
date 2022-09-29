;;execute NONMEM/pk2cpt.mod -dir NONMEM/pk2cpt
$PROBLEM pk2cpt Bayesian model

$DATA ../data/pk2cpt.csv IGNORE=@

$INPUT TIME ID II EVID DV CMT ADDL AMT


$SUB ADVAN4 TRANS4 	;; 2-comp with first order absorption

$PK
  TVKA = THETA(5)
  KA = TVKA

  TVCL = THETA(1)
  CL = TVCL
  MU_1 = LOG(TVCL)

  TVV2  = THETA(3)
  V2 = TVV2
  MU_2 = LOG(TVV2)

  TVQ = THETA(2)
  Q = TVQ

  TVV3  = THETA(4)
  V3 = TVV3

  S2 = V2  ;; Scaling compartment (check units of dose and observations)

$THETA
  (0, 7.4367958406427)    	;1 CL
  (0, 28.0799996152587)     ;2 Q
  (0, 78.4460632446725)   	;3 Vc
  (0, 68.1255965629187)     ;4 Vp
  (0, 1.0811298754049)     	;5 Ka

$OMEGA
 0.0 FIXED  ;; one subject

$SIGMA      ;; we are dealing them in the THETA block
  0.589695154260051  ;residual variability

$ERROR  	;; Calculation based on linear (non log-transformed) data
  IPRED = F
  Y= IPRED+ERR(1)

$EST MAXEVAL=99999 SIG=3 PRINT=5 NOABORT INTERACTION  	;; Estimation method

$COV PRINT=E UNCONDITIONAL