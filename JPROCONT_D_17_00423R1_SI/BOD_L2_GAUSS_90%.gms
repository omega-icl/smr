SETS
  VAR            model parameters                     / 1*2 /
  DATALL         all data point                       / 1*16 /
  DISCMAX        maximal SIP discretization           / 1*10 /
  DISC(DISCMAX)  dynamic SIP discretization;

* SELECT THE NUMBER OF DATA POINTS HERE:
SET
  DAT(DATALL)  4 data points                         / 4, 8, 12, 16 /;
* DAT(DATALL)  8 data points                         / 2, 4, 6, 8, 10, 12, 14, 16 /;
* DAT(DATALL)  16 data points                        / 1*16 /;

SCALARS
  CHI2           chi-square 2DOF 0.90                 / 4.61 /
  BODvar         measurement variance                 / 1.0  /
  obj_lseerr     perturbed least-squares value        / 0.0  /
  cvgtol         termination tolerance                / 1e-3 /;

PARAMETERS
  Tm(DATALL)  'time points [day]'
     / 1   0.5
       2   1.0
       3   1.5
       4   2.0
       5   2.5
       6   3.0
       7   3.5
       8   4.0
       9   4.5
      10   5.0
      11   5.5
      12   6.0
      13   6.5
      14   7.0
      15   7.5
      16   8.0 /
  BODm(DATALL)  'BOD data [mg/L]'
     / 1   3.696
       2   6.197
       3  11.289
       4  12.626
       5  13.589
       6  14.702
       7  17.105
       8  17.425
       9  17.810
      10  18.212
      11  18.625
      12  18.329
      13  19.138
      14  19.522
      15  20.768
      16  21.562 /

VARIABLES
  p(VAR)              model parameters
  err(DAT)            measurement error
  lse                 least-squares variable
  gap                 optimality gap variable;

PARAMETERS
  pdisc(VAR,DISCMAX)  SIP set discretization;

p.LO('1') = 0;  p.UP('1') = 50;  p.L('1') = 10;
p.LO('2') = 0;  p.UP('2') = 2;   p.L('2') = 1;
err.LO(DAT)  = -sqrt( CHI2 * BODvar );
err.UP(DAT)  =  sqrt( CHI2 * BODvar );

$MACRO BOD(T)    p('1') * ( 1 - exp( - p('2') * T ) )
$MACRO DBOD_DP1(T) ( 1 - exp( - p('2') * T ) )
$MACRO DBOD_DP2(T) p('1') * T * exp( - p('2') * T )
$MACRO BOD_DISC(DISC,T) pdisc('1',DISC) * ( 1 - exp( - pdisc('2',DISC) * T ) )

$MACRO BOD_VAL(T)      p.L('1') * ( 1 - exp( - p.L('2') * T ) )
$MACRO DBOD_DP1_VAL(T) ( 1 - exp( - p.L('2') * T ) )
$MACRO DBOD_DP2_VAL(T) p.L('1') * T * exp( - p.L('2') * T )

EQUATIONS
  obj_lse              least-squares objective
  err_conf             confidence region for measurement error
  nco_p1               first-order optimality condition for p1
  nco_p2               first-order optimality condition for p2
  optim_gap            optimality gap objective
  ctr_objerr(DISCMAX)  restriction on perturbed least-squares error;

obj_lse..    lse =E= SUM( DAT, sqr( BODm(DAT) - BOD( Tm(DAT) ) ) );
err_conf..   CHI2 * BODvar =G= SUM( DAT, sqr( err(DAT) ) );
nco_p1..     0 =E= SUM( DAT, ( BODm(DAT) + err(DAT) - BOD( Tm(DAT) ) ) * DBOD_DP1( Tm(DAT) ) );
nco_p2..     0 =E= SUM( DAT, ( BODm(DAT) + err(DAT) - BOD( Tm(DAT) ) ) * DBOD_DP2( Tm(DAT) ) );
optim_gap..  gap =E= obj_lseerr - SUM( DAT, sqr( BODm(DAT) + err.L(DAT) - BOD( Tm(DAT) ) ) );
ctr_objerr(DISC).. SUM( DAT, sqr( BODm(DAT) + err(DAT) - BOD_disc( DISC, Tm(DAT) ) ) ) =G= SUM( DAT, sqr( BODm(DAT) + err(DAT) - BOD( Tm(DAT) ) ) );

OPTION NLP = BARON;
*OPTION NLP = ANTIGONE;
OPTION OPTCR = 1e-3;
OPTION OPTCA = 1e-4;
OPTION RESLIM = 7200;


*******************************************************************************
* LEAST-SQUARES REGRESSION

MODEL lsq / obj_lse /;
SOLVE lsq USING NLP MINIMIZING lse;

PARAMETERS mle(VAR); mle(VAR) = p.L(var);
DISPLAY lse.L, mle;

SCALARS nco_p1_val, nco_p2_val;
nco_p1_val = SUM( DAT, ( BODm(DAT) + err.L(DAT) - BOD_VAL( Tm(DAT) ) ) * DBOD_DP1_VAL( Tm(DAT) ) );
nco_p2_val = SUM( DAT, ( BODm(DAT) + err.L(DAT) - BOD_VAL( Tm(DAT) ) ) * DBOD_DP2_VAL( Tm(DAT) ) );
DISPLAY nco_p1_val, nco_p2_val;
SCALAR loglr; loglr = - 0.5 * CHI2 - 0.5 * CARD(DAT) * log( 2 * PI * BODvar ) - 0.5 / BODvar * lse.L;
DISPLAY loglr;


*******************************************************************************
* SET-MEMBERSHIP REGRESSION

$MACRO SOLVE_SIP( lbd, obj, sense ) \
    gap.L = 2*cvgtol; \
    LOOP( DISCMAX$( gap.L > cvgtol ), \
       DISC(DISCMAX) = YES; \
       if( ORD(DISCMAX) = 1, pdisc(VAR,DISCMAX) = mle(VAR); \
       else                  pdisc(VAR,DISCMAX) = p.L(VAR); ); \
       if( sense = 0, SOLVE lbd USING NLP MINIMIZING obj; \
       else           SOLVE lbd USING NLP MAXIMIZING obj; ); \
       DISPLAY obj.L, p.L; \
       obj_lseerr = SUM( DAT, sqr( BODm(DAT) + err.L(DAT) - BOD_VAL( Tm(DAT) ) ) ); \
       DISPLAY obj_lseerr; \
       SOLVE viol USING NLP MAXIMIZING gap; \
       DISPLAY gap.L, p.L; \
    );


* LIKELIHOOD-CONTOUR ENCLOSURE

MODEL lkhcntr / all /;
MODEL viol / optim_gap /;
SOLVE_SIP( lkhcntr, lse, 1 );

SCALAR lce; lce = lse.L;
SCALAR loglce; loglce = - 0.5 * CARD(DAT) * log( 2 * PI * BODvar ) - 0.5 / BODvar * lce;
DISPLAY loglce;


* BOX ENCLOSURE BASED ON LIKELIHOOD-CONTOUR

$MACRO PUTVEC(v,vset) LOOP(vset, PUT ' ' v.L(vset):10:8 )
FILE results /PESMR_GAUSS_L2_90%.out/;
results.AP = 1;
PUT results;
PUT '#Run on ' system.date '  using source file  ' system.ifile /;

VARIABLES
  objvar;

EQUATIONS
  obj_p1  
  obj_p2
  lh_cntr    likelihood contour;

obj_p1..     objvar =E= p('1');
obj_p2..     objvar =E= p('2');
lh_cntr..    lce =G= SUM( DAT, sqr( BODm(DAT) - BOD( Tm(DAT) ) ) );

MODEL bnd_Tmin / obj_p1, lh_cntr /;
SOLVE bnd_Tmin USING NLP MINIMIZING objvar;
p.LO('1') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE bnd_Tmin USING NLP MAXIMIZING objvar;
p.UP('1') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax / obj_p2, lh_cntr /;
SOLVE bnd_Tmax USING NLP MINIMIZING objvar;
p.LO('2') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE bnd_Tmax USING NLP MAXIMIZING objvar;
p.UP('2') = objvar.L;
PUT ' ' objvar.L:12:10;
PUTCLOSE;

DISPLAY p.LO, p.UP;


* POLYHEDRAL ENCLOSURE BASED ON LIKELIHOOD-CONTOUR

$MACRO WID(X,I)   (X.UP(I)-X.LO(I))

SCALARS wid_p1, wid_p2;
wid_p1 = WID(p,'1');
wid_p2 = WID(p,'2'); 
DISPLAY wid_p1, wid_p2;

VARIABLES
  d1(VAR,VAR), d2(VAR,VAR);

EQUATIONS
  obj_p1_p_p2
  obj_p1_m_p2;

obj_p1_p_p2..   objvar =E= p('1')/(2*wid_p1) + p('2')/(2*wid_p2);
obj_p1_m_p2..   objvar =E= p('1')/(2*wid_p1) - p('2')/(2*wid_p2);

MODEL bnd_p1_p_p2 / obj_p1_p_p2, lh_cntr /;
SOLVE bnd_p1_p_p2 USING NLP MINIMIZING objvar;
d1.LO('1','2') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE bnd_p1_p_p2 USING NLP MAXIMIZING objvar;
d1.UP('1','2') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_p1_m_p2 / obj_p1_m_p2, lh_cntr /;
SOLVE bnd_p1_m_p2 USING NLP MINIMIZING objvar;
d2.LO('1','2') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE bnd_p1_m_p2 USING NLP MAXIMIZING objvar;
d2.UP('1','2') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

DISPLAY d1.LO, d1.UP, d2.LO, d2.UP;

