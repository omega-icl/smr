SETS
  VAR            model parameters                    / Tmin, Tmax, b, c /
  DAT            data                                / 1*15 /
  DISCMAX        maximal SIP discretization          / 1*20 /
  DISC(DISCMAX)  dynamic SIP discretization;

SCALARS
  CHI2         chi-square 4DOF 0.95                  / 9.488 /
  MUvar        measurement variance                  / 1e-2  /
  obj_lseerr   perturbed least-squares value         / 0.0   /
  cvgtol       termination tolerance                 / 1e-3  /;

PARAMETERS
  Tm(DAT)  temperature data
     / 1  294
       2  296
       3  298
       4  300
       5  302
       6  304
       7  306
       8  308
       9  310
      10  312
      11  314
      12  316
      13  318
      14  319
      15  320 /
  MUm(DAT)  growth rate data
     / 1  0.25
       2  0.56
       3  0.61
       4  0.79
       5  0.94
       6  1.04
       7  1.16
       8  1.23
       9  1.36
      10  1.32
      11  1.36
      12  1.34
      13  0.96
      14  0.83
      15  0.16 /;

VARIABLES
  p(VAR)   model parameters
  err(DAT) measurement error
  lse      least-squares variable
  gap      optimality gap variable;

PARAMETERS
  pdisc(VAR,DISCMAX)  SIP set discretization;

p.LO('Tmin') = 250;  p.UP('Tmin') = 300;
p.LO('Tmax') = 315;  p.UP('Tmax') = 330;
p.LO('b')    = 0.01; p.UP('b')    = 0.1;
p.LO('c')    = 0.01; p.UP('c')    = 1.0;
err.LO(DAT)  = -sqrt( CHI2 * MUvar );
err.UP(DAT)  =  sqrt( CHI2 * MUvar );

$MACRO MU(T)       sqr( p('b') * ( T - p('Tmin') ) * ( 1 - exp( p('c') * ( T - p('Tmax') ) ) ) )
$MACRO DMUDTMIN(T) ( - 2 * sqr( p('b') ) * ( T - p('Tmin') ) * sqr( 1 - exp( p('c') * ( T - p('Tmax') ) ) ) )
$MACRO DMUDTMAX(T) 2 * sqr( p('b') ) * sqr( T - p('Tmin') ) * ( 1 - exp( p('c') * ( T - p('Tmax') ) ) ) * p('c') * exp( p('c') * ( T - p('Tmax') ) )
$MACRO DMUDB(T)    2 * p('b') * sqr( T - p('Tmin') ) * sqr( 1 - exp( p('c') * ( T - p('Tmax') ) ) )
$MACRO DMUDC(T)    ( - 2 * sqr( p('b') ) * sqr( T - p('Tmin') ) * ( 1 - exp( p('c') * ( T - p('Tmax') ) ) ) * ( T - p('Tmax') ) * exp( p('c') * ( T - p('Tmax') ) ) )
$MACRO MUDISC(DISC,T) sqr( pdisc('b',DISC) * ( T - pdisc('Tmin',DISC) ) * ( 1 - exp( pdisc('c',DISC) * ( T - pdisc('Tmax',DISC) ) ) ) )

$MACRO MUVAL(T)    sqr( p.L('b') * ( T - p.L('Tmin') ) * ( 1 - exp( p.L('c') * ( T - p.L('Tmax') ) ) ) )
$MACRO DMUDTMINVAL(T) ( - 2 * sqr( p.L('b') ) * ( T - p.L('Tmin') ) * sqr( 1 - exp( p.L('c') * ( T - p.L('Tmax') ) ) ) )
$MACRO DMUDTMAXVAL(T) 2 * sqr( p.L('b') ) * sqr( T - p.L('Tmin') ) * ( 1 - exp( p.L('c') * ( T - p.L('Tmax') ) ) ) * p.L('c') * exp( p.L('c') * ( T - p.L('Tmax') ) )
$MACRO DMUDBVAL(T)    2 * p.L('b') * sqr( T - p.L('Tmin') ) * sqr( 1 - exp( p.L('c') * ( T - p.L('Tmax') ) ) )
$MACRO DMUDCVAL(T)    ( - 2 * sqr( p.L('b') ) * sqr( T - p.L('Tmin') ) * ( 1 - exp( p.L('c') * ( T - p.L('Tmax') ) ) ) * ( T - p.L('Tmax') ) * exp( p.L('c') * ( T - p.L('Tmax') ) ) )

EQUATIONS
  obj_lse              least-squares objective
  err_conf             confidence region for measurement error
  nco_Tmin             first-order optimality condition for Tmin
  nco_Tmax             first-order optimality condition for Tmax
  nco_b                first-order optimality condition for b
  nco_c                first-order optimality condition for c
  optim_gap            optimality gap objective
  ctr_objerr(DISCMAX)  restriction on perturbed least-squares error;

obj_lse..    lse =E= SUM( DAT, sqr( MUm(DAT) - MU( Tm(DAT) ) ) );
err_conf..   CHI2 * MUvar =G= SUM( DAT, sqr( err(DAT) ) );
nco_Tmin..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDTMIN( Tm(DAT) ) );
nco_Tmax..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDTMAX( Tm(DAT) ) );
nco_b..      0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDB( Tm(DAT) ) );
nco_c..      0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDC( Tm(DAT) ) );
optim_gap..  gap =E= obj_lseerr - SUM( DAT, sqr( MUm(DAT) + err.l(DAT) - MU( Tm(DAT) ) ) );
ctr_objerr(DISC).. SUM( DAT, sqr( MUm(DAT) + err(DAT) - MUDISC( DISC, Tm(DAT) ) ) ) =G= SUM( DAT, sqr( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) );

OPTION NLP = BARON;
*OPTION NLP = ANTIGONE;
OPTION OPTCR = 1e-4;
OPTION OPTCA = 1e-4;
OPTION RESLIM = 7200;


*******************************************************************************
* LEAST-SQUARES REGRESSION

MODEL lsq /obj_lse/;
SOLVE lsq USING NLP MINIMIZING lse;

PARAMETERS mle(VAR); mle(VAR) = p.L(var);
DISPLAY lse.L, mle;

SCALARS nco_Tmin_val, nco_Tmax_val, nco_b_val, nco_c_val;
nco_Tmin_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDTMINVAL( Tm(DAT) ) );
nco_Tmax_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDTMAXVAL( Tm(DAT) ) );
nco_b_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDBVAL( Tm(DAT) ) );
nco_c_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDCVAL( Tm(DAT) ) );
DISPLAY nco_Tmin_val, nco_Tmax_val, nco_b_val, nco_c_val;

SCALAR mlopt; mlopt = lse.L;
SCALAR loglr; loglr = - 0.5 * CHI2 - 0.5 * CARD(DAT) * log( 2 * PI * MUvar ) - 0.5 / MUvar * mlopt
DISPLAY loglr;


*******************************************************************************
* SET-MEMBERSHIP REGRESSION

$MACRO PUTVEC(v,vset) LOOP(vset, PUT ' ' v.L(vset):10:8 )
FILE results /SMR_RATK_L2_GAUSS_95%.out/;
results.AP = 1;
PUT results;
PUT '#Run on ' system.date '  using source file  ' system.ifile /;

$MACRO SOLVE_SIP( lbd, obj, sense ) \
    gap.L = 2*cvgtol; \
    LOOP( DISCMAX$( gap.L > cvgtol ), \
       DISC(DISCMAX) = YES; \
       if( ORD(DISCMAX) = 1, pdisc(VAR,DISCMAX) = mle(VAR); \
       else                  pdisc(VAR,DISCMAX) = p.L(VAR); ); \
       if( sense = 0, SOLVE lbd USING NLP MINIMIZING obj; \
       else           SOLVE lbd USING NLP MAXIMIZING obj; ); \
       DISPLAY obj.L, p.L; \
       obj_lseerr = SUM( DAT, sqr( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) ); \
       DISPLAY obj_lseerr; \
       SOLVE viol USING NLP MAXIMIZING gap; \
       DISPLAY gap.L, p.L; \
    );


* LIKELIHOOD-CONTOUR ENCLOSURE

MODEL lkhcntr / obj_lse, nco_Tmin, nco_Tmax, nco_b, nco_c, ctr_objerr, err_conf /;
MODEL viol / optim_gap /;
SOLVE_SIP( lkhcntr, lse, 1 );

SCALAR lceopt; lceopt = lse.L;
SCALAR loglce; loglce = - 0.5 * CARD(DAT) * log( 2 * PI * MUvar ) - 0.5 / MUvar * lceopt;
DISPLAY loglce;


* BOX ENCLOSURE

VARIABLES objvar, b(VAR);
EQUATIONS obj_Tmin, obj_Tmax, obj_b, obj_c;
obj_Tmin..   objvar =E= p('Tmin');
obj_Tmax..   objvar =E= p('Tmax');
obj_b..      objvar =E= p('b');
obj_c..      objvar =E= p('c');

MODEL bnd_Tmin /obj_Tmin,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin, objvar, 0 );
b.LO('Tmin') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin, objvar, 1 );
b.UP('Tmin') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax /obj_Tmax,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax, objvar, 0 );
b.LO('Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax, objvar, 1 );
b.UP('Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_b /obj_b,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_b, objvar, 0 );
b.LO('b') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_b, objvar, 1 );
b.UP('b') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_c /obj_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_c, objvar, 0 );
b.LO('c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_c, objvar, 1 );
b.UP('c') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

p.LO(VAR) = b.LO(VAR);
p.UP(VAR) = b.UP(VAR);
DISPLAY p.LO, p.UP;


* POLYHEDRAL ENCLOSURE

$MACRO WID(X,I)   (X.UP(I)-X.LO(I))
SCALARS wid_Tmin, wid_Tmax, wid_b, wid_c;
wid_Tmin = WID(p,'Tmin');
wid_Tmax = WID(p,'Tmax'); 
wid_b = WID(p,'b');
wid_c = WID(p,'c'); 
DISPLAY wid_Tmin, wid_Tmax, wid_b, wid_c;

VARIABLES d1(VAR,VAR), d2(VAR,VAR);
EQUATIONS obj_Tmin_p_Tmax, obj_Tmin_m_Tmax, obj_Tmin_p_b, obj_Tmin_m_b, obj_Tmin_p_c, obj_Tmin_m_c
          obj_Tmax_p_b, obj_Tmax_m_b, obj_Tmax_p_c, obj_Tmax_m_c, obj_b_p_c, obj_b_m_c;
obj_Tmin_p_Tmax..   objvar =E= p('Tmin')/(2*wid_Tmin) + p('Tmax')/(2*wid_Tmax);
obj_Tmin_m_Tmax..   objvar =E= p('Tmin')/(2*wid_Tmin) - p('Tmax')/(2*wid_Tmax);
obj_Tmin_p_b..      objvar =E= p('Tmin')/(2*wid_Tmin) + p('b')/(2*wid_b);
obj_Tmin_m_b..      objvar =E= p('Tmin')/(2*wid_Tmin) - p('b')/(2*wid_b);
obj_Tmin_p_c..      objvar =E= p('Tmin')/(2*wid_Tmin) + p('c')/(2*wid_c);
obj_Tmin_m_c..      objvar =E= p('Tmin')/(2*wid_Tmin) - p('c')/(2*wid_c);
obj_Tmax_p_b..      objvar =E= p('Tmax')/(2*wid_Tmax) + p('b')/(2*wid_b);
obj_Tmax_m_b..      objvar =E= p('Tmax')/(2*wid_Tmax) - p('b')/(2*wid_b);
obj_Tmax_p_c..      objvar =E= p('Tmax')/(2*wid_Tmax) + p('c')/(2*wid_c);
obj_Tmax_m_c..      objvar =E= p('Tmax')/(2*wid_Tmax) - p('c')/(2*wid_c);
obj_b_p_c..         objvar =E= p('b')/(2*wid_b) + p('c')/(2*wid_c);
obj_b_m_c..         objvar =E= p('b')/(2*wid_b) - p('c')/(2*wid_c);

MODEL bnd_Tmin_p_Tmax /obj_Tmin_p_Tmax,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_Tmax, objvar, 0 );
d1.LO('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_Tmax, objvar, 1 );
d1.UP('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_Tmax /obj_Tmin_m_Tmax,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_Tmax, objvar, 0 );
d2.LO('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_Tmax, objvar, 1 );
d2.UP('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmin_p_b /obj_Tmin_p_b,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_b, objvar, 0 );
d1.LO('Tmin','b') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_b, objvar, 1 );
d1.UP('Tmin','b') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_b /obj_Tmin_m_b,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_b, objvar, 0 );
d2.LO('Tmin','b') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_b, objvar, 1 );
d2.UP('Tmin','b') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmin_p_c /obj_Tmin_p_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_c, objvar, 0 );
d1.LO('Tmin','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_c, objvar, 1 );
d1.UP('Tmin','c') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_c /obj_Tmin_m_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_c, objvar, 0 );
d2.LO('Tmin','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_c, objvar, 1 );
d2.UP('Tmin','c') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmax_p_b /obj_Tmax_p_b,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_p_b, objvar, 0 );
d1.LO('Tmax','b') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_p_b, objvar, 1 );
d1.UP('Tmax','b') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax_m_b /obj_Tmax_m_b,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_m_b, objvar, 0 );
d2.LO('Tmax','b') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_m_b, objvar, 1 );
d2.UP('Tmax','b') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmax_p_c /obj_Tmax_p_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_p_c, objvar, 0 );
d1.LO('Tmax','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_p_c, objvar, 1 );
d1.UP('Tmax','c') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax_m_c /obj_Tmax_m_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_m_c, objvar, 0 );
d2.LO('Tmax','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_m_c, objvar, 1 );
d2.UP('Tmax','c') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_b_p_c /obj_b_p_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_b_p_c, objvar, 0 );
d1.LO('b','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_b_p_c, objvar, 1 );
d1.UP('b','c') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_b_m_c /obj_b_m_c,nco_Tmin,nco_Tmax,nco_b,nco_c,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_b_m_c, objvar, 0 );
d2.LO('b','c') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_b_m_c, objvar, 1 );
d2.UP('b','c') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

DISPLAY d1.LO, d1.UP, d2.LO, d2.UP;

