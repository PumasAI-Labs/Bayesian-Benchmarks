;;execute NONMEM/pk2cpt.mod -dir NONMEM/pk2cpt -parafile=/opt/NONMEM/nm751/run/mpilinux8.pnm -nodes=4
$PROBLEM pk2cpt Bayesian model

$DATA ../data/pk2cpt.csv IGNORE=@

$INPUT TIME ID II EVID DV CMT ADDL AMT

$SUB ADVAN4 TRANS4 	;; 2-comp with first order absorption

$PK

  MU_1 = THETA(1)
  TVCL = EXP(MU_1)
  CL = TVCL

  MU_2 = THETA(2)
  TVV2  = EXP(MU_2)
  V2 = TVV2

  MU_3 = THETA(3)
  TVQ = EXP(MU_3)
  Q = TVQ

  MU_4 = THETA(4)
  TVV3  = EXP(MU_4)
  V3 = TVV3
  
  MU_5 = THETA(5)
  TVKA = EXP(MU_5)
  KA = TVKA

  S2 = V2  ;; Scaling compartment (check units of dose and observations)

$THETA     ;; initial values
  2.0           	     ;1 CL log(7.4367958406427)
  3.33                 ;2 Q  log(28.0799996152587)
  4.36                 ;3 Vc log(78.4460632446725)
  4.22                 ;4 Vp log(68.1255965629187)
  0.078                ;5 Ka log(1.0811298754049)

$OMEGA
 0.4 FIX  ;; one subject

$SIGMA
  0.589695154260051  ;residual variability

$ERROR  	;; Calculation based on linear (log-transformed) data
  IPRED = log(F)
  Y = IPRED + ERR(1)

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

;; 4 chains in parallels
$EST METHOD=NUTS AUTO=0 PRINT=100 NBURN=1000 SEED=1 NITER=1000 NUTS_MAXDEPTH=10 NUTS_DELTA=0.8 FILE=/dev/null
;;$EST METHOD=BAYES AUTO=1 PRINT=100 NBURN=1000 NITER=1000

$COV PRINT=E UNCONDITIONAL