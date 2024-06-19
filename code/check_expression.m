



N:=61;
d:=12;
Nd:=N*d;
k:=2;
SB:=Floor(k*Index(Gamma0(Nd))/12)+2;


load "../data/big_expressions/732.txt";


HeckeKernels:=[
[2,1]
];



M1:=ModularForms(N,k);
N1:=NewSubspace(CuspidalSubspace(M1));
M2:=ModularForms(Nd,k);

if #HeckeKernels eq 0 then
  f:=N1.1;
else
  f:=M1!N1!((&meet[Kernel(HeckeOperator(N1,u[1]) + u[2]) : u in HeckeKernels ]).1);
end if;



P<q>:=PowerSeriesRing(Rationals(),SB);

ETA:=function(N,list)
  prod:=P!1;
  D:=Divisors(N);
  for u in [1..#D] do
    prod*:=&*[1-q^(D[u]*i) : i in [1..SB]]^list[u];
  end for;
  return prod*q^(&+[list[i]*D[i] : i in [1..#D]]/24);
end function;



eta_forms:=[M2!ETA(Nd,qq) : qq in quots];

vecs:=[Vector(ElementToSequence(u)) : u in eta_forms];
v:=Vector(ElementToSequence(M2!f));

vecs:=[ChangeRing(u,Rationals()) : u in vecs];
v:=ChangeRing(v,Rationals());


&+[vecs[i]*coeffs[i] : i in [1..#vecs]] eq v;







//
