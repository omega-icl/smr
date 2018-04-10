SETS
  VAR            model parameters                    / Tmin, Tmax, Topt, Mopt /
  DAT            data                                / 1*15 /
  DISCMAX        maximal SIP discretization          / 1*10 /
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
p.LO('Topt') = 290;  p.UP('Topt') = 330;
p.LO('Mopt') = 1.0;  p.UP('Mopt') = 2.0;
err.LO(DAT)  = -sqrt( CHI2 * MUvar );
err.UP(DAT)  =  sqrt( CHI2 * MUvar );

$MACRO MU(T)       p('Mopt') * ( 1 - sqr( T - p('Topt') ) / ( sqr( T - p('Topt') ) + T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') ) )
$MACRO DMUDTMIN(T) p('Mopt') * sqr( T - p('Topt') ) * ( T - p('Tmax') ) / sqr( sqr( T - p('Topt') ) + T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') )
$MACRO DMUDTMAX(T) p('Mopt') * sqr( T - p('Topt') ) * ( T - p('Tmin') ) / sqr( sqr( T - p('Topt') ) + T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') )
$MACRO DMUDTOPT(T) 2 * p('Mopt') * ( T - p('Topt') ) * ( T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') ) / sqr( sqr( T - p('Topt') ) + T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') )
$MACRO DMUDMOPT(T) ( 1 - sqr( T - p('Topt') ) / ( sqr( T - p('Topt') ) + T * ( p('Tmin') + p('Tmax') - T ) - p('Tmin')*p('Tmax') ) )
$MACRO MUDISC(DISC,T)  pdisc('Mopt',DISC) * ( 1 - sqr( T - pdisc('Topt',DISC) ) / ( sqr( T - pdisc('Topt',DISC) ) + T * ( pdisc('Tmin',DISC) + pdisc('Tmax',DISC) - T ) - pdisc('Tmin',DISC)*pdisc('Tmax',DISC) ) )

$MACRO MUVAL(T)    p.L('Mopt') * ( 1 - sqr( T - p.L('Topt') ) / ( sqr( T - p.L('Topt') ) + T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') ) )
$MACRO DMUDTMINVAL(T) p.L('Mopt') * sqr( T - p.L('Topt') ) * ( T - p.L('Tmax') ) / sqr( sqr( T - p.L('Topt') ) + T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') )
$MACRO DMUDTMAXVAL(T) p.L('Mopt') * sqr( T - p.L('Topt') ) * ( T - p.L('Tmin') ) / sqr( sqr( T - p.L('Topt') ) + T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') )
$MACRO DMUDTOPTVAL(T) 2 * p.L('Mopt') * ( T - p.L('Topt') ) * ( T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') ) / sqr( sqr( T - p.L('Topt') ) + T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') )
$MACRO DMUDMOPTVAL(T) ( 1 - sqr( T - p.L('Topt') ) / ( sqr( T - p.L('Topt') ) + T * ( p.L('Tmin') + p.L('Tmax') - T ) - p.L('Tmin')*p.L('Tmax') ) )

EQUATIONS
  obj_lse              least-squares objective
  err_conf             confidence region for measurement error
  nco_Tmin             first-order optimality condition for Tmin
  nco_Tmax             first-order optimality condition for Tmax
  nco_Topt             first-order optimality condition for Topt
  nco_Mopt             first-order optimality condition for Mopt
  optim_gap            optimality gap objective
  ctr_objerr(DISCMAX)  restriction on perturbed least-squares error;

obj_lse..    lse =E= SUM( DAT, sqr( MUm(DAT) - MU( Tm(DAT) ) ) );
err_conf..   CHI2 * MUvar =G= SUM( DAT, sqr( err(DAT) ) );
nco_Tmin..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDTMIN( Tm(DAT) ) );
nco_Tmax..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDTMAX( Tm(DAT) ) );
nco_Topt..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDTOPT( Tm(DAT) ) );
nco_Mopt..   0 =E= SUM( DAT, ( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) * DMUDMOPT( Tm(DAT) ) );
optim_gap..  gap =E= obj_lseerr - SUM( DAT, sqr( MUm(DAT) + err.l(DAT) - MU( Tm(DAT) ) ) );
ctr_objerr(DISC).. SUM( DAT, sqr( MUm(DAT) + err(DAT) - MUDISC( DISC, Tm(DAT) ) ) ) =G= SUM( DAT, sqr( MUm(DAT) + err(DAT) - MU( Tm(DAT) ) ) );

OPTION NLP = BARON;
*OPTION NLP = ANTIGONE;
OPTION OPTCR = 1e-3;
OPTION OPTCA = 1e-4;
OPTION RESLIM = 7200;


*******************************************************************************
* LEAST-SQUARES REGRESSION

MODEL lsq /obj_lse/;
SOLVE lsq USING NLP MINIMIZING lse;

PARAMETERS mle(VAR); mle(VAR) = p.L(var);
DISPLAY lse.L, mle;

SCALARS nco_Tmin_val, nco_Tmax_val, nco_Topt_val, nco_Mopt_val;
nco_Tmin_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDTMINVAL( Tm(DAT) ) );
nco_Tmax_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDTMAXVAL( Tm(DAT) ) );
nco_Topt_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDTOPTVAL( Tm(DAT) ) );
nco_Mopt_val = SUM( DAT, ( MUm(DAT) + err.L(DAT) - MUVAL( Tm(DAT) ) ) * DMUDMOPTVAL( Tm(DAT) ) );
DISPLAY nco_Tmin_val, nco_Tmax_val, nco_Topt_val, nco_Mopt_val;

SCALAR mlopt; mlopt = lse.L;
SCALAR loglr; loglr = - 0.5 * CHI2 - 0.5 * CARD(DAT) * log( 2 * PI * MUvar ) - 0.5 / MUvar * mlopt
DISPLAY loglr;


*******************************************************************************
* SET-MEMBERSHIP REGRESSION

$MACRO PUTVEC(v,vset) LOOP(vset, PUT ' ' v.L(vset):10:8 )
FILE results /SMR_CARD_L2_GAUSS_95%.out/;
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

MODEL lkhcntr / obj_lse, nco_Tmin, nco_Tmax, nco_Topt, nco_Mopt, ctr_objerr, err_conf /;
MODEL viol / optim_gap /;
SOLVE_SIP( lkhcntr, lse, 1 );

SCALAR lceopt; lceopt = lse.L;
SCALAR loglce; loglce = - 0.5 * CARD(DAT) * log( 2 * PI * MUvar ) - 0.5 / MUvar * lceopt;
DISPLAY loglce;


* BOX ENCLOSURE

VARIABLES objvar, b(VAR);
EQUATIONS obj_Tmin, obj_Tmax, obj_Topt, obj_Mopt;
obj_Tmin..   objvar =E= p('Tmin');
obj_Tmax..   objvar =E= p('Tmax');
obj_Topt..   objvar =E= p('Topt');
obj_Mopt..   objvar =E= p('Mopt');

MODEL bnd_Tmin /obj_Tmin,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin, objvar, 0 );
b.LO('Tmin') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin, objvar, 1 );
b.UP('Tmin') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax /obj_Tmax,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax, objvar, 0 );
b.LO('Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax, objvar, 1 );
b.UP('Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Topt /obj_Topt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Topt, objvar, 0 );
b.LO('Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Topt, objvar, 1 );
b.UP('Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Mopt /obj_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Mopt, objvar, 0 );
b.LO('Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Mopt, objvar, 1 );
b.UP('Mopt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

p.LO(VAR) = b.LO(VAR);
p.UP(VAR) = b.UP(VAR);
DISPLAY p.LO, p.UP;


* POLYHEDRAL ENCLOSURE

$MACRO WID(X,I)   (X.UP(I)-X.LO(I))
SCALARS wid_Tmin, wid_Tmax, wid_Topt, wid_Mopt;
wid_Tmin = WID(p,'Tmin');
wid_Tmax = WID(p,'Tmax'); 
wid_Topt = WID(p,'Topt');
wid_Mopt = WID(p,'Mopt'); 
DISPLAY wid_Tmin, wid_Tmax, wid_Topt, wid_Mopt;

VARIABLES d1(VAR,VAR), d2(VAR,VAR);
EQUATIONS obj_Tmin_p_Tmax, obj_Tmin_m_Tmax, obj_Tmin_p_Topt, obj_Tmin_m_Topt, obj_Tmin_p_Mopt, obj_Tmin_m_Mopt
          obj_Tmax_p_Topt, obj_Tmax_m_Topt, obj_Tmax_p_Mopt, obj_Tmax_m_Mopt, obj_Topt_p_Mopt, obj_Topt_m_Mopt;
obj_Tmin_p_Tmax..   objvar =E= p('Tmin')/(2*wid_Tmin) + p('Tmax')/(2*wid_Tmax);
obj_Tmin_m_Tmax..   objvar =E= p('Tmin')/(2*wid_Tmin) - p('Tmax')/(2*wid_Tmax);
obj_Tmin_p_Topt..   objvar =E= p('Tmin')/(2*wid_Tmin) + p('Topt')/(2*wid_Topt);
obj_Tmin_m_Topt..   objvar =E= p('Tmin')/(2*wid_Tmin) - p('Topt')/(2*wid_Topt);
obj_Tmin_p_Mopt..   objvar =E= p('Tmin')/(2*wid_Tmin) + p('Mopt')/(2*wid_Mopt);
obj_Tmin_m_Mopt..   objvar =E= p('Tmin')/(2*wid_Tmin) - p('Mopt')/(2*wid_Mopt);
obj_Tmax_p_Topt..   objvar =E= p('Tmax')/(2*wid_Tmax) + p('Topt')/(2*wid_Topt);
obj_Tmax_m_Topt..   objvar =E= p('Tmax')/(2*wid_Tmax) - p('Topt')/(2*wid_Topt);
obj_Tmax_p_Mopt..   objvar =E= p('Tmax')/(2*wid_Tmax) + p('Mopt')/(2*wid_Mopt);
obj_Tmax_m_Mopt..   objvar =E= p('Tmax')/(2*wid_Tmax) - p('Mopt')/(2*wid_Mopt);
obj_Topt_p_Mopt..   objvar =E= p('Topt')/(2*wid_Topt) + p('Mopt')/(2*wid_Mopt);
obj_Topt_m_Mopt..   objvar =E= p('Topt')/(2*wid_Topt) - p('Mopt')/(2*wid_Mopt);

MODEL bnd_Tmin_p_Tmax /obj_Tmin_p_Tmax,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_Tmax, objvar, 0 );
d1.LO('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_Tmax, objvar, 1 );
d1.UP('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_Tmax /obj_Tmin_m_Tmax,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_Tmax, objvar, 0 );
d2.LO('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_Tmax, objvar, 1 );
d2.UP('Tmin','Tmax') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmin_p_Topt /obj_Tmin_p_Topt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_Topt, objvar, 0 );
d1.LO('Tmin','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_Topt, objvar, 1 );
d1.UP('Tmin','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_Topt /obj_Tmin_m_Topt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_Topt, objvar, 0 );
d2.LO('Tmin','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_Topt, objvar, 1 );
d2.UP('Tmin','Topt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmin_p_Mopt /obj_Tmin_p_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_p_Mopt, objvar, 0 );
d1.LO('Tmin','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_p_Mopt, objvar, 1 );
d1.UP('Tmin','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmin_m_Mopt /obj_Tmin_m_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmin_m_Mopt, objvar, 0 );
d2.LO('Tmin','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmin_m_Mopt, objvar, 1 );
d2.UP('Tmin','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmax_p_Topt /obj_Tmax_p_Topt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_p_Topt, objvar, 0 );
d1.LO('Tmax','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_p_Topt, objvar, 1 );
d1.UP('Tmax','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax_m_Topt /obj_Tmax_m_Topt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_m_Topt, objvar, 0 );
d2.LO('Tmax','Topt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_m_Topt, objvar, 1 );
d2.UP('Tmax','Topt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Tmax_p_Mopt /obj_Tmax_p_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_p_Mopt, objvar, 0 );
d1.LO('Tmax','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_p_Mopt, objvar, 1 );
d1.UP('Tmax','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Tmax_m_Mopt /obj_Tmax_m_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Tmax_m_Mopt, objvar, 0 );
d2.LO('Tmax','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Tmax_m_Mopt, objvar, 1 );
d2.UP('Tmax','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

MODEL bnd_Topt_p_Mopt /obj_Topt_p_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Topt_p_Mopt, objvar, 0 );
d1.LO('Topt','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Topt_p_Mopt, objvar, 1 );
d1.UP('Topt','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

MODEL bnd_Topt_m_Mopt /obj_Topt_m_Mopt,nco_Tmin,nco_Tmax,nco_Topt,nco_Mopt,ctr_objerr,err_conf/;
SOLVE_SIP( bnd_Topt_m_Mopt, objvar, 0 );
d2.LO('Topt','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10;

SOLVE_SIP( bnd_Topt_m_Mopt, objvar, 1 );
d2.UP('Topt','Mopt') = objvar.L;
PUT ' ' objvar.L:12:10 /;
PUTCLOSE;

DISPLAY d1.LO, d1.UP, d2.LO, d2.UP;

