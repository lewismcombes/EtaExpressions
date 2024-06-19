

// ensure you have already loaded your list of eta quotients of level N*d

// N is the level of the modular form you care about, Nd is the level the quotients live in
// these need not be the same. d is the factor added to the level.
// k is the weight

// HeckeKernels are used to find the exact form f you care about in M_k(N)
// in general your form is determined by a series of kernels of Hecke operators
// the entry in HeckeKernels [p,a] means f is in Kernel( T_p + a )
// if the list is empty, the first basis vector of the newspace is used, so only leave it empty
// when your f generates the whole newspace.

// SB is the sturm bound for the space


N:=94;
d:=10;
Nd:=N*d;
k:=2;
SB:=Floor(k*Index(Gamma0(Nd))/12)+2;

// unfortunately magma doesn't support dynamic loading so we can't make this string and then
// feed it to the load command. you have to update this for yourself every time.
load "..data/eta_quotients/940.txt";


HeckeKernels:=[
[3,0]
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




RR:=RSpace(CoefficientRing(vecs[1]),Degree(vecs[1]));
v in sub<RR|vecs>;






//
