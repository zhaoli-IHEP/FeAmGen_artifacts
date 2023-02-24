CFunctions DeltaFun(symmetric), DeltaAdj(symmetric), SUNTrace(cyclesymmetric);
CFunctions SUNT, SUNF, SUNTConj, sunTraceConj, sunTrace;
CFunctions SUNTChain, SUNTraceChain;

symbols cla0,...,cla100;
symbols claC0,...,claC100;
symbols claM0,...,claM100;
symbols clb0,...,clb100;
symbols clbC0,...,clbC100;
symbols clbM0,...,clbM100;
symbols clv0,...,clv10000;
symbols clw0,...,clw10000;

Symbols colorX, colorY;

Symbols n0,...,n100;
Symbols m0,...,m100;
Symbols I, im, ca, cf;

*----------------------------------------
#procedure calc1_CF()
**** we use the definition f(a,b,c) = -2*i*Tr(a,b,c)+2*i*Tr(c,b,a) as same as defined in MadGraph etc.

repeat id SUNTConj(?n1,n2?,n3?) = SUNT(reverse_(?n1),n3,n2);
repeat id sunTraceConj(?n0) = sunTrace(reverse_(?n0));
.sort

*
* Trim the cyclesymmetric of SUNTrace
* NB: But we should have replaced SUNTrace by sunTrace in MIRACLE master code.
*multiply replace_(SUNTrace,sunTrace);
*.sort

repeat;
  id SUNF(m1?,m2?,m3?) = -2*I*sunTrace(m1,m2,m3)+2*I*sunTrace(m3,m2,m1);
  id SUNT(m0?,n1?,n2?)*SUNT(m0?,n3?,n4?) = 1/2*( DeltaFun(n1,n4)*DeltaFun(n2,n3)-1/ca*DeltaFun(n1,n2)*DeltaFun(n3,n4) );
endrepeat;
.sort

repeat;
  id DeltaFun(n1?,n0?)*SUNT(m0?,n0?,n2?) = SUNT(m0,n1,n2);

  id DeltaFun(n1?,n0?)*DeltaFun(n0?,n2?) = DeltaFun(n1,n2);
  id DeltaFun(n0?,n0?) = ca;
  id DeltaAdj(n1?,n0?)*DeltaAdj(n0?,n2?) = DeltaAdj(n1,n2);
  id DeltaAdj(n0?,n0?) = 2*ca*cf;

  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?,?m2) = sunTrace(?m1,m3,?m2);
  id DeltaAdj(m0?,m3?)*SUNT(?m1,m0?,?m2,n1?,n2?) = SUNT(?m1,m3,?m2,n1,n2);

  id SUNT(?m3,m0?,?m4,n1?,n2?)*sunTrace(?m1,m0?,?m2) = 1/2*( SUNT(?m3,?m2,?m1,?m4,n1,n2)-1/ca*SUNT(?m3,?m4,n1,n2)*sunTrace(?m2,?m1) );

  id SUNT(?m1,m0?,?m2,n1?,n2?)*SUNT(?m3,m0?,?m4,n3?,n4?) 
    = 1/2*( SUNT(?m1,?m4,n1,n4)*SUNT(?m3,?m2,n3,n2) -1/ca*SUNT(?m1,?m2,n1,n2)*SUNT(?m3,?m4,n3,n4) );
  id SUNT(n1?,n2?) = DeltaFun(n1,n2);

  id SUNT(?m1,n1?,n2?)*SUNT(?m2,n2?,n3?) = SUNT(?m1,?m2,n1,n3);
  id SUNT(?m1,n1?,n1?) = sunTrace(?m1);
  id SUNT(m0?,m0?,n1?,n2?) = cf*DeltaFun(n1,n2);
  id SUNT(?m1,m0?,m0?,?m2,n1?,n2?) = cf*SUNT(?m1,?m2,n1,n2);

  id SUNT(m0?,?m2,m0?,n1?,n2?) = 1/2*( DeltaFun(n1,n2)*sunTrace(?m2) - 1/ca*SUNT(?m2,n1,n2) );

  id SUNT(?m1,m0?,?m2,m0?,?m3,n1?,n2?) = 1/2*( SUNT(?m1,?m3,n1,n2)*sunTrace(?m2) - 1/ca*SUNT(?m1,?m2,?m3,n1,n2) );
  id sunTrace(m0?) = 0;

  id sunTrace(m1?,m1?) = ca*cf;
  id sunTrace(m1?,m2?) = 1/2*DeltaAdj(m1,m2);
  id sunTrace(?m1,m0?,?m2)*sunTrace(?m3,m0?,?m4) = 1/2*( sunTrace(?m1,?m4,?m3,?m2)-1/ca*sunTrace(?m2,?m1)*sunTrace(?m4,?m3) );
endrepeat;
.sort

*
* Recover the cyclesymmetric of SUNTrace
*
multiply replace_(sunTrace,SUNTrace);
.sort


#endprocedure




*----------------------------------------
#procedure calc2_CF()

*
* Trim the cyclesymmetric of SUNTrace
*
multiply replace_(SUNTrace,sunTrace);
.sort

repeat id SUNF(m1?,m2?,m3?) = -2*I*sunTrace(m1,m2,m3)+2*I*sunTrace(m3,m2,m1);
.sort

repeat;
  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?) = sunTrace(?m1,m3);
  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?,?m2) = sunTrace(?m1,m3,?m2);
  id DeltaAdj(m0?,m3?)*SUNT(?m1,m0?,?m2) = SUNT(?m1,m3,?m2);
  id DeltaAdj(m3?,m0?)*SUNT(?m1,m0?,?m2) = SUNT(?m1,m3,?m2);

  id SUNT(?m1,n1?,n2?)*DeltaFun(n2?,n3?) = SUNT(?m1,n1,n3);
  id SUNT(?m1,n1?,n2?)*DeltaFun(n3?,n2?) = SUNT(?m1,n1,n3);

  id DeltaFun(n1?,n2?)*DeltaFun(n2?,n3?) = DeltaFun(n1,n3);
  id DeltaFun(n1?,n1?) = ca;
  id DeltaAdj(n0?,n0?) = 2*ca*cf;

  id SUNT(?m1,n1?,n2?)*SUNT(?m2,n2?,n3?) = SUNT(?m1,?m2,n1,n3);
  id SUNT(?m0,m1?,n1?,n1?) = sunTrace(?m0,m1);

  id sunTrace(m1?,m1?) = ca*cf;
  id sunTrace(m1?,m2?) = 1/2*DeltaAdj(m1,m2);

  id sunTrace(?m1,m0?,m0?,?m2) = cf*sunTrace(?m1,?m2);

  id sunTrace(?m1,m0?,?m2,m0?,?m3) = 1/2*( sunTrace(?m1,?m3)*sunTrace(?m2) - 1/ca*sunTrace(?m1,?m2,?m3) );

  id sunTrace(m0?) = 0;

  id DeltaAdj(n1?,n0?)*DeltaAdj(n0?,n2?) = DeltaAdj(n1,n2);
  id sunTrace(?m1,m0?,?m2)*sunTrace(?m3,m0?,?m4) = 1/2*( sunTrace(?m1,?m4,?m3,?m2)-1/ca*sunTrace(?m2,?m1)*sunTrace(?m4,?m3) );
endrepeat;
.sort

id sunTrace() = ca;
.sort

*
* Recover the cyclesymmetric of SUNTrace
*
multiply replace_(sunTrace,SUNTrace);
.sort

#endprocedure

