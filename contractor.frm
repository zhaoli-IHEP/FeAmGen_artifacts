
symbol diim, epsilon;
dimension diim;

vector q1,q2,q3,q4;
vector k1,...,k10;
vector K1,...,K10;
vector mom,mom0,mom1,...,mom10;
vector unity;

Set Kn: K1,...,K10;
Set kn: k1,...,k10;

set LOOP: q1,q2,q3,q4;
set NonLOOP: k1,...,k10,K1,...,K10;
set MOM: k1,...,k10,K1,...,K10,q1,q2,q3,q4;

symbol pi, im, sqrt2, shat;
auto symbol ver, gc;

***Dirac indices
auto symbol spa, spb, spv spw;

***Lorentz indices
auto index mua, mub, muv, muw, rho, dum, epMU;

Set RHO: rho,rho0,...,rho200;
Set EPMU: epMU1,...,epMU100;
Set DUM: dum1,...,dum100;
Set NonEPMU: mua,mua0,...,mua100,mub,mub0,...,mub100,muv, muv0,...,muv10000,muw, muw0,...,muw10000,dum1,...,dum100;

Set LOR: mua,mua0,...,mua100,mub,mub0,...,mub100,
         muv,muv0,...,muv10000,muw,muw0,...,muw10000,
         rho,rho0,...,rho100,dum,dum0,...,dum100,
         epMU1,...,epMU100;


***For matching both indices and momenta
auto index var;
***For matching numbers
auto symbol int, setint;

symbol mass, mass1,...,mass4, width;

CFunction Den, PowDen, FV, FermionLoopPow, GhostLoopPow;
CFunction VecEpsilon, VecEpsilon1,...,VecEpsilon4, VecEp, VecEpC;
CFunction Spinor, Spinor1,...,Spinor4, FermionChain;

CFunction SpUB, SpVB, SpU, SpV;
Set SPSET: SpUB, SpVB, SpU, SpV;
Set LSPSET: SpUB, SpVB;
Set RSPSET: SpU, SpV;

CFunction UB, VB, U, V;
Set INSPSET: UB, VB, U, V;
Set ILSPSET: UB, VB;
Set IRSPSET: U, V;

CFunction SP(symmetric), LMT(symmetric), Levi(antisymmetric);
CFunction SPC(symmetric);
CFunction GAij, PLij, PRij, ONEij(symmetric), Trace, Trace5, Trace5sym;
CFunction GA, GA5, PL, PR;

CFunction JJ(antisymmetric), FF(antisymmetric);

symbol dZ3x1, dZ3x2; 
symbol dtZ3x1, dtZ3x2; 
symbol dZ2tx1, dZ2tx2; 
symbol dZmtx1, dZmtx2; 
symbol dZ2x1, dZ2x2; 
symbol dZ1x1, dZ1x2; 
symbol dZ4x1, dZ4x2; 
symbol dtZ1x1, dtZ1x2; 
symbol dZ1Fx1, dZ1Fx2; 
symbol dZ1Ftx1, dZ1Ftx2; 
symbol dZgx1, dZgx2;
symbol d2s1, d2s2, d2s3, d2s4;
symbol f2s1, f2s2, f2s3, f2s4;



*----------------------------------------
#procedure Simplification()

id im^2 = -1;
id sqrt2^2 = 2;
id sqrt2^(-2) = 1/2;

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort
*
* trivial Lorentz indices contraction
*
repeat;
  id LMT(rho1?,rho1?) = diim;
  id LMT(rho1?,rho2?)^2 = diim;
  id LMT(rho1?,rho2?)*LMT(rho2?,rho3?) = LMT(rho1,rho3);

  id LMT(rho1?,rho2?)*FV(mom1?,rho1?) = FV(mom1,rho2);
  id LMT(rho1?,rho2?)*GAij(spa1?,spa2?,rho1?) = GAij(spa1,spa2,rho2);
  id LMT(rho1?,rho2?)*Levi(rho1?,rho3?,rho4?,rho5?) = Levi(rho2,rho3,rho4,rho5);
  id LMT(rho1?,rho2?)*VecEpsilon?{VecEp,VecEpC}(int1?,rho2?,mom1?,mass1?) = VecEpsilon(int1,rho1,mom1,mass1);
endrepeat;

id GAij(spa1?,spa2?,rho?)*FV(mom?,rho?) = GAij(spa1,spa2,mom);
id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);

repeat;
  id FV(mom?,rho?)*Levi(var1?,var2?,var3?,rho?) = Levi(var1,var2,var3,mom);
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho2?),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
endrepeat;

repeat id LMT(rho?NonEPMU,rho0?)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = LMT(EPMU[int],rho0)*VecEpsilon(int,EPMU[int],?vars);
id GAij(spa1?,spa2?,rho?NonEPMU)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = GAij(spa1,spa2,EPMU[int])*VecEpsilon(int,EPMU[int],?vars);
.sort

*
* vanishing 
*
id FV(mom?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom?,mass?) = 0;
.sort

*
* use EPMU indices
*
id SP( FV(mom1?,0), VecEpsilon?{VecEp,VecEpC}(int?,0,mom2?,mass?) )
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id FV(mom1?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom2?,mass?)
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id VecEpsilon1?{VecEp,VecEpC}(int1?, rho?, ?vars1) * VecEpsilon2?{VecEp,VecEpC}(int2?, rho?, ?vars2)
  = SP( VecEpsilon1(int1, 0, ?vars1), VecEpsilon2(int2, 0, ?vars2) );
.sort

*
*Linearly expand momentum polynomial in FV and GA
*
id FV(var?,rho?) = FV(var,rho);
id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
id SP(rho1?,rho2?) = SP(rho1,rho2);
id Levi(rho1?,rho2?,rho3?,rho4?) = Levi(rho1,rho2,rho3,rho4);
.sort

*
* vanishing momentum scalar product also vanishes
*
id SP(mom?NULL,mom?NULL) = 0;

*
* Explain FermionLoopPow and GhostLoopPow
*
id FermionLoopPow(-1,int?) = (-1)^int;
id GhostLoopPow(-1,int?) = (-1)^int;
.sort

id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho?ALLLOR), ?vars2 ) = diim*FermionChain(?vars1,?vars2);
.sort


#endprocedure








*----------------------------------------
*** without expansion, we do not expand momentum in FV or GA.
#procedure SimplificationNoExpand()

id im^2 = -1;
id sqrt2^2 = 2;
id sqrt2^(-2) = 1/2;

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort
*
* trivial Lorentz indices contraction
*
repeat;
  id LMT(rho1?,rho1?) = diim;
  id LMT(rho1?,rho2?)^2 = diim;
  id LMT(rho1?,rho2?)*LMT(rho2?,rho3?) = LMT(rho1,rho3);

  id LMT(rho1?,rho2?)*FV(mom1?,rho1?) = FV(mom1,rho2);
  id LMT(rho1?,rho2?)*GAij(spa1?,spa2?,rho1?) = GAij(spa1,spa2,rho2);
  id LMT(rho1?,rho2?)*Levi(rho1?,rho3?,rho4?,rho5?) = Levi(rho2,rho3,rho4,rho5);
  id LMT(rho1?,rho2?)*VecEpsilon?{VecEp,VecEpC}(int1?,rho2?,mom1?,mass1?) = VecEpsilon(int1,rho1,mom1,mass1);
endrepeat;

id GAij(spa1?,spa2?,rho?)*FV(mom?,rho?) = GAij(spa1,spa2,mom);
id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);

repeat;
  id FV(mom?,rho?)*Levi(var1?,var2?,var3?,rho?) = Levi(var1,var2,var3,mom);
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho2?),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
endrepeat;

repeat id LMT(rho?NonEPMU,rho0?)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = LMT(EPMU[int],rho0)*VecEpsilon(int,EPMU[int],?vars);
id GAij(spa1?,spa2?,rho?NonEPMU)*VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,?vars)
  = GAij(spa1,spa2,EPMU[int])*VecEpsilon(int,EPMU[int],?vars);
.sort

*
* vanishing 
*
id FV(mom?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom?,mass?) = 0;
.sort

*
* use EPMU indices
*
id SP( FV(mom1?,0), VecEpsilon?{VecEp,VecEpC}(int?,0,mom2?,mass?) )
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id FV(mom1?,rho?) * VecEpsilon?{VecEp,VecEpC}(int?,rho?,mom2?,mass?)
  = FV(mom1,EPMU[int]) * VecEpsilon(int,EPMU[int],mom2,mass);
id VecEpsilon1?{VecEp,VecEpC}(int1?, rho?, ?vars1) * VecEpsilon2?{VecEp,VecEpC}(int2?, rho?, ?vars2)
  = SP( VecEpsilon1(int1, 0, ?vars1), VecEpsilon2(int2, 0, ?vars2) );
.sort

*
*Linearly expand momentum polynomial in FV and GA
*
***id FV(var?,rho?) = FV(var,rho);
***id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
***id SP(rho1?,rho2?) = SP(rho1,rho2);
***id Levi(rho1?,rho2?,rho3?,rho4?) = Levi(rho1,rho2,rho3,rho4);
***.sort

*
* vanishing momentum scalar product also vanishes
*
id SP(mom?NULL,mom?NULL) = 0;

*
* Explain FermionLoopPow and GhostLoopPow
*
id FermionLoopPow(-1,int?) = (-1)^int;
id GhostLoopPow(-1,int?) = (-1)^int;
.sort


#endprocedure




























*----------------------------------------
#procedure ArrangeTrace()

repeat;
  id Trace(?vars1,PL,PL,?vars2) = Trace(?vars1,PL,?vars2);
  id Trace(?vars1,PR,PR,?vars2) = Trace(?vars1,PR,?vars2);
  id Trace(?vars1,PL,PR,?vars2) = 0;
  id Trace(?vars1,PR,PL,?vars2) = 0;

  id Trace(?vars1,GA(rho?ALLLOR),PL,?vars2) = Trace(?vars1,PR,GA(rho),?vars2);
  id Trace(?vars1,GA(rho?ALLLOR),PR,?vars2) = Trace(?vars1,PL,GA(rho),?vars2);

  id Trace(?vars1,GA(mom?),PL,?vars2) = Trace(?vars1,PR,GA(mom),?vars2);
  id Trace(?vars1,GA(mom?),PR,?vars2) = Trace(?vars1,PL,GA(mom),?vars2);
endrepeat;
.sort


repeat;
  id FV(mom?,rho?ALLLOR)*Trace(?vars1,GA(rho?ALLLOR),?vars2) = Trace(?vars1,GA(mom),?vars2);
  id LMT(rho1?ALLLOR,rho2?ALLLOR)*Trace(?vars1,GA(rho2?ALLLOR),?vars2) = Trace(?vars1,GA(rho1),?vars2);

  id Trace(?vars1,GA(mom?ALLMOM),GA(mom?ALLMOM),?vars2) = SP(mom,mom)*Trace(?vars1,?vars2);
  id SP(mom?NULL,mom?NULL) = 0;

  id Trace(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = Trace(?vars1,?vars2)*diim;
endrepeat;
.sort
 
id Trace(PL,?vars) = (1+sign_(nargs_(?vars)))/2*Trace(PL,?vars);
id Trace(PR,?vars) = (1+sign_(nargs_(?vars)))/2*Trace(PR,?vars);
id Trace(GA(mom?),?vars) = (1+sign_(1+nargs_(?vars)))/2*Trace(GA(mom),?vars);
id Trace(GA(rho?ALLLOR),?vars) = (1+sign_(1+nargs_(?vars)))/2*Trace(GA(rho),?vars);
.sort

id Trace(PL,?vars) = 1/2*Trace(?vars)-1/2*Trace5(?vars);
id Trace(PR,?vars) = 1/2*Trace(?vars)+1/2*Trace5(?vars);
.sort

***
***
repeat;
  id once Trace5(?vars) = 1/24*e_(rho100,rho101,rho102,rho103)*Trace(GA(rho100),GA(rho101),GA(rho102),GA(rho103),?vars);
  sum rho100;
  sum rho101;
  sum rho102;
  sum rho103;
endrepeat;
.sort
***
***

repeat id Trace(?vars1, GA(var?), ?vars2) = Trace(?vars1, var, ?vars2);
.sort

repeat id Trace(?vars1, unity, ?vars2 ) = Trace( ?vars1, ?vars2 );
.sort

repeat;
  id once Trace(?vars) = g_(1,?vars);
  tracen, 1;
endrepeat;
.sort

contract;
.sort

id VecEpsilon?{VecEp,VecEpC}(int?,mom0?,mom?,mass?) = FV(mom0,EPMU[int])*VecEpsilon(int,EPMU[int],mom,mass);
id mom?NULL.mom?NULL = 0;
id mom1?.mom2? = SP(mom1,mom2);
id mom?(rho?ALLLOR) = FV(mom,rho);
id e_(rho1?,rho2?,rho3?,rho4?) = -im*Levi(rho1,rho2,rho3,rho4);
.sort

id d_(rho1?,rho2?) = LMT(rho1,rho2);
repeat id LMT(rho1?ALLLOR,rho2?ALLLOR)*FermionChain(?vars1,GA(rho2?ALLLOR),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
.sort

id VecEpsilon?{VecEp,VecEpC}(int?,rho?NonEPMU,mom?,mass?)*FermionChain(?vars1,GA(rho?NonEPMU),?vars2) 
  = VecEpsilon(int,EPMU[int],mom,mass)*FermionChain(?vars1,GA(EPMU[int]),?vars2);
.sort

#endprocedure









  
#procedure SimpleOrdering()

repeat;
  id FermionChain(?vars1,PL,PL,?vars2) = FermionChain(?vars1,PL,?vars2);
  id FermionChain(?vars1,PR,PR,?vars2) = FermionChain(?vars1,PR,?vars2);
  id FermionChain(?vars1,PL,PR,?vars2) = 0;
  id FermionChain(?vars1,PR,PL,?vars2) = 0;

  id FermionChain(?vars1,GA(rho?),PL,?vars2) = FermionChain(?vars1,PR,GA(rho),?vars2);
  id FermionChain(?vars1,GA(rho?),PR,?vars2) = FermionChain(?vars1,PL,GA(rho),?vars2);
endrepeat;
.sort

*
* handle the simple structure in FermionChain
*
repeat;
  id FermionChain( ?vars1, GA(mom?), GA(mom?), ?vars2 ) = FermionChain( ?vars1, ?vars2 )*SP(mom,mom);
  id SP(mom?NULL,mom?NULL) = 0;
  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho?ALLLOR), ?vars2 ) = FermionChain(?vars1,?vars2)*diim;
endrepeat;
.sort

***pull momentum of right-side spinor to the relevant spinor
repeat;
  id FermionChain( ?vars1, GA(mom?), GA(rho?ALLLOR), ?vars2, Spinor?IRSPSET(int?,mom?,mass?) )
    = FermionChain( ?vars1, ?vars2, Spinor(int,mom,mass) )*2*FV(mom,rho)
     -FermionChain( ?vars1, GA(rho), GA(mom), ?vars2, Spinor(int,mom,mass) );

  id FV(mom?NULL,rho?)^2 = 0;
  id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);

  repeat id FermionChain( ?vars1, GA(mom?), GA(mom?), ?vars2 ) = FermionChain( ?vars1, ?vars2 )*SP(mom,mom);
  id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);
  id SP(mom?NULL,mom?NULL) = 0;

  id FermionChain( ?vars1, GA(mom1?), GA(mom2?ALLMOM), ?vars2, Spinor?{U,V}(int?,mom1?,mass?) )
    = FermionChain( ?vars1, ?vars2, Spinor(int,mom1,mass) )*2*SP(mom1,mom2)
     -FermionChain( ?vars1, GA(mom2), GA(mom1), ?vars2, Spinor(int,mom1,mass) );

  id FermionChain( ?vars, GA(mom?), U(int?,mom?,mass?) ) = mass*FermionChain( ?vars, U(int,mom,mass) );
  id FermionChain( ?vars, GA(mom?), V(int?,mom?,mass?) ) = -mass*FermionChain( ?vars, V(int,mom,mass) );
endrepeat;
.sort

***pull momentum of left-side spinor to the relevant spinor
repeat;
  id FermionChain( Spinor?ILSPSET(int?,mom?,mass?), ?vars1, GA(rho?ALLLOR), GA(mom?), ?vars2 )
    = FermionChain( Spinor(int,mom,mass), ?vars1, ?vars2 )*2*FV(mom,rho)
     -FermionChain( Spinor(int,mom,mass), ?vars1, GA(mom), GA(rho), ?vars2 );

  id FV(mom?NULL,rho?)^2 = 0;
  id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);

  id FermionChain( Spinor?ILSPSET(int?,mom?,mass?), ?vars1, PR, GA(mom?), ?vars2 )
    = FermionChain( Spinor(int,mom,mass), ?vars1, GA(mom), PL, ?vars2 );

  id FermionChain( Spinor?ILSPSET(int?,mom?,mass?), ?vars1, PL, GA(mom?), ?vars2 )
    = FermionChain( Spinor(int,mom,mass), ?vars1, GA(mom), PR, ?vars2 );

  repeat id FermionChain( ?vars1, GA(mom?), GA(mom?), ?vars2 ) = FermionChain( ?vars1, ?vars2 )*SP(mom,mom);
  id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);
  id SP(mom?NULL,mom?NULL) = 0;

  id FermionChain( Spinor?ILSPSET(int?,mom1?,mass?), ?vars1, GA(mom2?ALLMOM), GA(mom1?), ?vars2 )
    = FermionChain( Spinor(int,mom1,mass), ?vars1, ?vars2 )*2*SP(mom1,mom2)
     -FermionChain( Spinor(int,mom1,mass), ?vars1, GA(mom1), GA(mom2), ?vars2 );

  id FermionChain( UB(int?,mom?,mass?), GA(mom?), ?vars ) 
    = mass*FermionChain( UB(int,mom,mass), ?vars );
  id FermionChain( VB(int?,mom?,mass?), GA(mom?), ?vars ) 
    = -mass*FermionChain( VB(int,mom,mass), ?vars );
endrepeat;
.sort


***move momentums to right-side of lorentz indices
repeat;
  id FermionChain( ?vars1, GA(mom?ALLMOM), GA(mom?ALLMOM), ?vars2 ) = FermionChain(?vars1,?vars2)*SP(mom,mom);
  id SP(mom?NULL,mom?NULL) = 0;

  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho?ALLLOR), ?vars2 ) = FermionChain(?vars1,?vars2)*diim;

  id FermionChain( ?vars1, GA(mom?), GA(rho?ALLLOR), ?vars2 )
    = 2*FV(mom,rho)*FermionChain(?vars1,?vars2)
     - FermionChain( ?vars1, GA(rho), GA(mom), ?vars2 );
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) 
    = FermionChain(?vars1,GA(mom),?vars2);
endrepeat;
.sort


***
*** Contract dummy indices in single fermion chain.
***
repeat;
  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho?ALLLOR), ?vars2 )
    = diim*FermionChain(?vars1,?vars2);

  id FermionChain( ?vars1, GA(rho2?ALLLOR), GA(rho1?), GA(rho2?ALLLOR), ?vars2 )
    = (2-diim)*FermionChain( ?vars1, GA(rho1), ?vars2 );

  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho1?ALLLOR), GA(rho2?ALLLOR), GA(rho?ALLLOR), ?vars2 )
    = 4*LMT(rho1,rho2)*FermionChain( ?vars1, ?vars2 )
     + (diim-4)*FermionChain( ?vars1, GA(rho1), GA(rho2), ?vars2 );

  id LMT(rho1?,rho2?)*FermionChain( ?vars1, GA(rho2?), ?vars2 ) = FermionChain( ?vars1, GA(rho1), ?vars2 );

  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho1?), GA(rho2?), GA(rho3?), GA(rho?ALLLOR), ?vars2 )
    = -2*FermionChain( ?vars1, GA(rho3), GA(rho2), GA(rho1), ?vars2 )
     + (4-diim)*FermionChain( ?vars1, GA(rho1), GA(rho2), GA(rho3), ?vars2 );

  id FermionChain( ?vars1, GA(rho?ALLLOR), GA(rho1?), GA(rho2?), GA(rho3?), GA(rho4?), GA(rho?ALLLOR), ?vars2 )
    = 2*FermionChain( ?vars1, GA(rho2), GA(rho3), GA(rho4), GA(rho1), ?vars2 )
     + 2*FermionChain( ?vars1, GA(rho1), GA(rho4), GA(rho3), GA(rho2), ?vars2 )
     + (diim-4)*FermionChain( ?vars1, GA(rho1), GA(rho2), GA(rho3), GA(rho4), ?vars2 );
endrepeat;
.sort

*
* Substitute dummy indices (except epMU) using system Nm_?, and look for largest number of dummy index Nm_? by iterating over each term
*
while( match(FermionChain(?vars1,GA(rho?NonEPMU$LORENTZ),?vars2)) );
  sum $LORENTZ;
endwhile;
.sort

*
*move epMU to right-side of the other lorentz indices 
*
repeat;
  id FermionChain(?vars1,GA(rho1?EPMU),GA(rho2?dummyindices_),?vars2)
    = 2*LMT(rho1,rho2)*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(rho2),GA(rho1),?vars2);
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho2?),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
endrepeat;
.sort

*
* rearrange the sequence of epMU in FermionChain
* We assume epMU is no more than 10
*
repeat;
  #do MUIDX1 = 1, 9, 1 
    #do MUIDX2 = `MUIDX1'+1, 10, 1 
      id FermionChain(?vars1,GA(epMU`MUIDX2'),GA(epMU`MUIDX1'),?vars2)
        = 2*LMT(epMU`MUIDX1',epMU`MUIDX2')*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(epMU`MUIDX1'),GA(epMU`MUIDX2'),?vars2);
    #enddo
  #enddo
endrepeat;
.sort

*
* renumber the dummy indices Nm_? in order to have better format.
* this may lead to some cost, but further check has to be made.
*
renumber 0;
.sort

*
* rearrange the sequence of dummy Lorentz indices in FermionChain
* We assume dummyindices_ no more than 20
*
repeat;
  #do MUIDX1 = 1, 19, 1
    #do MUIDX2 = `MUIDX1'+1, 20, 1
      id FermionChain(?vars1,GA(N`MUIDX2'_?),GA(N`MUIDX1'_?),?vars2)
        = 2*LMT(N`MUIDX1'_?,N`MUIDX2'_?)*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(N`MUIDX1'_?),GA(N`MUIDX2'_?),?vars2);
    #enddo
  #enddo
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho1?),?vars2) = FermionChain(?vars1,GA(rho2),?vars2);
endrepeat;
.sort
*
* Replace system dummy indices Nm_? by our dummy indices dum in case to read back to GiNaC.
* We assume this should give the canonical form of FermionChain, 
*   since it seems dummy indices Nm_? can make canonical form of an expression automatically.
*
#do MUIDX = 1, 20, 1
  Multiply replace_(N`MUIDX'_?,dum`MUIDX');
#enddo
.sort


*
* move ki to the left of Ki in FermionChain
*
repeat;
  id FermionChain(?vars1,GA(K1?MASSIVE),GA(k1?NULL),?vars2)
    = 2*SP(k1,K1)*FermionChain(?vars1,?vars2) - FermionChain(?vars1,GA(k1),GA(K1),?vars2);
endrepeat;
*
* Rearrange the sequence of ki in FermionChain
* We assume no more than 10
*
repeat;
  #do MOMIDX1 = 1, 9, 1
    #do MOMIDX2 = `MOMIDX1'+1, 10, 1
      id FermionChain(?vars1, GA(k`MOMIDX2'), GA(k`MOMIDX1'), ?vars2)
        = 2*SP(k`MOMIDX1',k`MOMIDX2')*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(k`MOMIDX1'),GA(k`MOMIDX2'),?vars2);
    #enddo
  #enddo
endrepeat;
*
* Rearrange the sequence of Ki in FermionChain
* We assume no more than 10
*
repeat;
  #do MOMIDX1 = 1, 9, 1
    #do MOMIDX2 = `MOMIDX1'+1, 10, 1
      id FermionChain(?vars1, GA(K`MOMIDX2'), GA(K`MOMIDX1'), ?vars2)
        = 2*SP(K`MOMIDX1',K`MOMIDX2')*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(K`MOMIDX1'),GA(K`MOMIDX2'),?vars2);
    #enddo
  #enddo
endrepeat;

*
* Rearrange the sequence of qi in FermionChain
* We assume no more than 4
*
repeat;
  #do MOMIDX1 = 1, 3, 1
    #do MOMIDX2 = `MOMIDX1'+1, 4, 1
      id FermionChain(?vars1, GA(q`MOMIDX2'), GA(q`MOMIDX1'), ?vars2)
        = 2*SP(q`MOMIDX1',q`MOMIDX2')*FermionChain(?vars1,?vars2)-FermionChain(?vars1,GA(q`MOMIDX1'),GA(q`MOMIDX2'),?vars2);
    #enddo
  #enddo
endrepeat;


*
* Contract adjacent momentum slash or lorent indices in FermionChain
*
repeat;
  id FermionChain(?vars1,GA(mom?ALLMOM),GA(mom?ALLMOM),?vars2) = FermionChain(?vars1,?vars2)*SP(mom,mom);
  id FermionChain(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,?vars2)*diim;
endrepeat;
.sort


#endprocedure
  
  
  
  






















































*----------------------------------------
#procedure contractDiracIndices()

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id PLij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PLij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*PRij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PRij(spa2?,spa3?) = 0;
  id PRij(spa1?,spa2?)*PLij(spa2?,spa3?) = 0;


  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort

***id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
***.sort



*
* now chainin the Dirac objects in FermionChain according to Dirac indices
*
***Here var? could be mom? or rho?
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,mom?) = FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );

***flip
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,mom?) = -FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );


id Spinor1?LSPSET[setint1](int1?,spa?,?var1)*Spinor2?RSPSET[setint2](int2?,spa?,?var2) = FermionChain( ILSPSET[setint1](int1,?var1), IRSPSET[setint2](int2,?var2) );
.sort

repeat;
***Here var? could be mom? or rho?
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,mom?) = FermionChain( ?vars, GA(mom), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa1?,spa2?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa1?,spa2?) = FermionChain( ?vars, PR, spa2 );

***flip
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,mom?) = -FermionChain( ?vars, GA(mom), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa2?,spa1?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa2?,spa1?) = FermionChain( ?vars, PR, spa2 );
endrepeat;

id FermionChain(?vars1,spa?)*Spinor?RSPSET[setint](int?,spa?,?vars2) = FermionChain( ?vars1, IRSPSET[setint](int,?vars2) );
.sort


*
* Look for Trace
*
repeat;
  id once, GAij(spa1?,spa2?,var?) = Trace(GA(var),spa1,spa2);
  repeat;
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,mom?) = Trace(?vars,GA(mom),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,rho?ALLLOR) = Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa2?,spa3?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa2?,spa3?) = Trace(?vars,PR,spa1,spa3);

    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,mom?) = -Trace(?vars,GA(mom),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,rho?ALLLOR) = -Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa3?,spa2?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa3?,spa2?) = Trace(?vars,PR,spa1,spa3);
  endrepeat;
  id Trace(?vars,spa?symbol_,spa?symbol_) = Trace(?vars);
  id PLij(spa?symbol_,spa?symbol_) = 2;
  id PRij(spa?symbol_,spa?symbol_) = 2; 
  id ONEij(spa?symbol_,spa?symbol_) = 4;
endrepeat;
.sort

#call ArrangeTrace();

***move PL and PR to left-side of FermionChain, right-side of spinor
repeat;
  id FermionChain( ?vars1, PL, PR, ?vars2 ) = 0;
  id FermionChain( ?vars1, PR, PL, ?vars2 ) = 0;
  id FermionChain( ?vars1, PL, PL, ?vars2 ) = FermionChain( ?vars1, PL, ?vars2 );
  id FermionChain( ?vars1, PR, PR, ?vars2 ) = FermionChain( ?vars1, PR, ?vars2 );

  id FermionChain( ?vars1, GA(mom?), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(mom), ?vars2 );
  id FermionChain( ?vars1, GA(mom?), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(mom), ?vars2 );

  id FermionChain( ?vars1, GA(rho?ALLLOR), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(rho), ?vars2 );
  id FermionChain( ?vars1, GA(rho?ALLLOR), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(rho), ?vars2 );
endrepeat;
.sort


repeat;
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);
  id LMT(rho1?ALLLOR,rho2?ALLLOR)*FermionChain(?vars1,GA(rho1?ALLLOR),?vars2) = FermionChain(?vars1,GA(rho2),?vars2);
endrepeat;
.sort

repeat;
  id FermionChain(?vars1,GA(mom?),GA(mom?),?vars2) = SP(mom,mom)*FermionChain(?vars1,?vars2);
  id SP(mom?NULL,mom?NULL) = 0;
  id FermionChain(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = diim*FermionChain(?vars1,?vars2);

*** Dirac equation for U and V (UB and VB) 
  id FermionChain( ?vars, GA(mom?), U(int?,mom?,0) ) = 0;
  id FermionChain( ?vars, GA(mom?), V(int?,mom?,0) ) = 0;

  id FermionChain( UB(int?,mom?,0), GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;


endrepeat;
.sort

#endprocedure













*----------------------------------------
*** without expansion, the commutative relation between PL/PR and GA is uncertain, since GA[mom] may have mass*unity here.
#procedure contractDiracIndicesNoExpand()

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA5)/2, GA(6)
* PL-->tagF, (1-GA5)/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id PLij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PLij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*PRij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PRij(spa2?,spa3?) = 0;
  id PRij(spa1?,spa2?)*PLij(spa2?,spa3?) = 0;


  id Spinor?SPSET(int?,spa1?,mom?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,mass);
endrepeat;
.sort

***id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
***.sort



*
* now chainin the Dirac objects in FermionChain according to Dirac indices
*
***Here var? could be mom? or rho?
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,mom?) = FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );

***flip
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,mom?) = -FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );


id Spinor1?LSPSET[setint1](int1?,spa?,?var1)*Spinor2?RSPSET[setint2](int2?,spa?,?var2) = FermionChain( ILSPSET[setint1](int1,?var1), IRSPSET[setint2](int2,?var2) );
.sort

repeat;
***Here var? could be mom? or rho?
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,?expr) = FermionChain( ?vars, GA(?expr), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa1?,spa2?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa1?,spa2?) = FermionChain( ?vars, PR, spa2 );

***flip
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,?expr) = -FermionChain( ?vars, GA(?expr), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa2?,spa1?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa2?,spa1?) = FermionChain( ?vars, PR, spa2 );
endrepeat;

id FermionChain(?vars1,spa?)*Spinor?RSPSET[setint](int?,spa?,?vars2) = FermionChain( ?vars1, IRSPSET[setint](int,?vars2) );
.sort


*
* Look for Trace
*
repeat;
  id once, GAij(spa1?,spa2?,?expr) = Trace(GA(?expr),spa1,spa2);
  repeat;
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,?expr) = Trace(?vars,GA(?expr),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,rho?ALLLOR) = Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa2?,spa3?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa2?,spa3?) = Trace(?vars,PR,spa1,spa3);

    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,?expr) = -Trace(?vars,GA(?expr),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,rho?ALLLOR) = -Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa3?,spa2?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa3?,spa2?) = Trace(?vars,PR,spa1,spa3);
  endrepeat;
  id Trace(?vars,spa?symbol_,spa?symbol_) = Trace(?vars);
  id PLij(spa?symbol_,spa?symbol_) = 2;
  id PRij(spa?symbol_,spa?symbol_) = 2; 
  id ONEij(spa?symbol_,spa?symbol_) = 4;
endrepeat;
.sort

*** #call ArrangeTrace();

***move PL and PR to left-side of FermionChain, right-side of spinor
repeat;
  id FermionChain( ?vars1, PL, PR, ?vars2 ) = 0;
  id FermionChain( ?vars1, PR, PL, ?vars2 ) = 0;
  id FermionChain( ?vars1, PL, PL, ?vars2 ) = FermionChain( ?vars1, PL, ?vars2 );
  id FermionChain( ?vars1, PR, PR, ?vars2 ) = FermionChain( ?vars1, PR, ?vars2 );

**id FermionChain( ?vars1, GA(mom?), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(mom), ?vars2 );
**id FermionChain( ?vars1, GA(mom?), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(mom), ?vars2 );

**id FermionChain( ?vars1, GA(rho?ALLLOR), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(rho), ?vars2 );
**id FermionChain( ?vars1, GA(rho?ALLLOR), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(rho), ?vars2 );
endrepeat;
.sort


repeat;
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);
  id LMT(rho1?ALLLOR,rho2?ALLLOR)*FermionChain(?vars1,GA(rho1?ALLLOR),?vars2) = FermionChain(?vars1,GA(rho2),?vars2);
endrepeat;
.sort

repeat;
  id FermionChain(?vars1,GA(mom?),GA(mom?),?vars2) = SP(mom,mom)*FermionChain(?vars1,?vars2);
  id SP(mom?NULL,mom?NULL) = 0;
  id FermionChain(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = diim*FermionChain(?vars1,?vars2);

*** Dirac equation for U and V (UB and VB) 
  id FermionChain( ?vars, GA(mom?), U(int?,mom?,0) ) = 0;
  id FermionChain( ?vars, GA(mom?), V(int?,mom?,0) ) = 0;

  id FermionChain( UB(int?,mom?,0), GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PL, GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,0), PR, GA(mom?), ?vars ) = 0;
endrepeat;
.sort

#endprocedure




