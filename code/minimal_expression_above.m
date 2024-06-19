

// N is the level of the modular form you care about, Nd is the level the quotients live in
// these need not be the same. d is the factor added to the level.
// k is the weight

// HeckeKernels are used to find the exact form f you care about in M_k(N)
// in general your form is determined by a series of kernels of Hecke operators
// the entry in HeckeKernels [p,a] means f is in Kernel( T_p + a )
// if the list is empty, the first basis vector of the newspace is used, so only leave it empty
// when your f generates the whole newspace.

// SB is the sturm bound for the space


// computs the total height of all the rationals in the vector vec
h:=function(vec)
  AA:=ElementToSequence(vec);
  return &+[Log(Abs(Numerator(u))) + Log(Abs(Denominator(u))) : u in AA | u ne 0];
end function;

// returns the number of non-zero entries of vec
nz:=function(vec)
  AA:=ElementToSequence(vec);
  return [i : i in [1..#AA] | AA[i] ne 0];
end function;

// returns only those basis elements that correspond to non-zero entries of vec
NonZeroifyBasis:=function(basis,vec)
  a:=Solution(Matrix(basis),vec);
  return [basis[i] : i in nz(a)];
end function;


N:=37;
d:=8;
Nd:=N*d;
k:=2;
SB:=Floor(k*Index(Gamma0(Nd))/12)+2;


load "../data/eta_quotients/296.txt";


HeckeKernels:=[
[2,2]
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


// we find a bad solution to start with
basis:=[];
is_in_space:=false;
i:=1;
while not is_in_space do
  Append(~basis,vecs[i]);
  i+:=1;
  is_in_space:=v in sub<RR|basis>;
end while;

// these variables keep track of the current record holder for minimality
record:=#basis;
best_vecs:=basis;
sol:=Solution(Matrix(basis),v);


// just in case you happen to immediately find the best possible solution
[quots[Index(vecs,u)] : u in best_vecs];
ElementToSequence(sol);



// this number is silly and large. making it very large will enure it spends
// a long time trying to make the expression smaller.
// from experience, this seems to have diminishing returns after a while
for i in [1..10^9] do

  basis:=[];

  space:=sub<RR|basis>;

  is_in_space:=false;

  while not is_in_space do
    new_vec:=Random([u : u in vecs | not u in basis]);
    Append(~basis,new_vec);
    space:=sub<RR|basis>;
    is_in_space:=v in space;
  end while;

  basis:=NonZeroifyBasis(basis,v);
  a:=Solution(Matrix(basis),v);

  // if we have found a better expression, we print it to the terminal
  if (#basis lt record) or (#basis eq record and h(a) lt h(sol)) then
    record:=#basis;
    best_vecs:=basis;
    sol:=a;
    record;
    h(sol);
    "";
    a:=Solution(Matrix(best_vecs),v);
    [quots[Index(vecs,u)] : u in best_vecs];
    ElementToSequence(a);
    "";
  end if;

end for;



//
