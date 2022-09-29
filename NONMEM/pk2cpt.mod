;;execute NONMEM/pk2cpt.mod -dir NONMEM/pk2cpt -parafile=/opt/NONMEM/nm751/run/mpilinux8.pnm -nodes=4
$PROBLEM pk2cpt Bayesian model

$DATA ../data/pk2cpt.csv IGNORE=@

$INPUT TIME ID II EVID DV CMT ADDL AMT

$SUB ADVAN4 TRANS4 	;; 2-comp with first order absorption

$PK
  TVKA = THETA(5)
  KA = TVKA

  TVCL = THETA(1) * EXP(ETA(1))
  CL = TVCL
  MU_1 = LOG(TVCL)

  TVV2  = THETA(3)
  V2 = TVV2

  TVQ = THETA(2)
  Q = TVQ

  TVV3  = THETA(4)
  V3 = TVV3

  S2 = V2  ;; Scaling compartment (check units of dose and observations)

$THETA     ;; initial values
  7.4367958406427  	   ;1 CL
  28.0799996152587     ;2 Q
  78.4460632446725     ;3 Vc
  68.1255965629187     ;4 Vp
  1.0811298754049      ;5 Ka

$OMEGA
 0.0 FIX  ;; one subject

$SIGMA
  0.589695154260051  ;residual variability

$ERROR  	;; Calculation based on linear (non log-transformed) data
  IPRED = F
  Y= IPRED+ERR(1)

;; bayesian part
$PRIOR NWPRI
$THETAP           ;; prior information of THETAS
  (2.3 FIX)          ; CL log(10)
  (2.7 FIX)          ; Q  log(15)
  (3.56 FIX)         ; Vc log(35)
  (4.7 FIX)          ; Vp log(105)
  (0.92 FIX)         ; Ka log(2.5)
$THETAPV BLOCK(5) ;; variances for priors on THETAS (var-cov)
0.25 FIX             ; CL
0.0 0.5              ; Q
0.0 0.0 0.25         ; Vc
0.0 0.0 0.0 0.5      ; Vp
0.0 0.0 0.0 0.0 1.0  ; Ka

$OMEGAP (0.0 FIX)  ;; prior information for OMEGA
$OMEGAPD (1.0 FIX) ;; df for OMEGA prior  
$SIGMAP (0.0 FIX)  ;; prior information for SIGMA
$SIGMAPD (1.0 FIX) ;; df for SIGMA prior

;;$EST MAXEVAL=99999 SIG=3 PRINT=5 NOABORT INTERACTION  	;; Estimation method
$EST METHOD=NUTS AUTO=0 PRINT=20 NBURN=1000 NITER=2000 NUTS_MAXDEPTH=10 NUTS_DELTA=0.8

$COV PRINT=E UNCONDITIONAL