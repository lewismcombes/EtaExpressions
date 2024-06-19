

// N is the level of the modular form you care about, Nd is the level the quotients live in
// these need not be the same. d is the factor added to the level.
// k is the weight

// HeckeKernels are used to find the exact form f you care about in M_k(N)
// in general your form is determined by a series of kernels of Hecke operators
// the entry in HeckeKernels [p,a] means f is in Kernel( T_p + a )
// if the list is empty, the first basis vector of the newspace is used, so only leave it empty
// when your f generates the whole newspace.

// SB is the sturm bound for the space

N:=43;
d:=8;
Nd:=N*d;
k:=2;
SB:=Floor(k*Index(Gamma0(Nd))/12)+2;


load "../data/eta_quotients/344.txt";


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





// lists all the binomial choices of n Choose k.
// this function very quickly makes magma very slow, probably
// because it makes a list with many elements in it. there is
// surely a better way to do this. if you think you know it,
// please email me (or make a pull request!)
Binom:=function(n,k)
  if k eq 0 then
    return [[]];
  elif k eq 1 then
    return [[i] : i in [1..n]];
  elif n lt k then
    return [];
  else
    old1:=$$(n-1,k);
    old2:=$$(n-1,k-1);

    return [u cat [n] : u in old2] cat [u : u in old1 | not n in u];
  end if;
end function;

HH:=function(rat)
  a:=Numerator(rat);
  b:=Denominator(rat);
  return Log(Abs(a*b));
end function;




// this is here to make sure things don't go awry
// if you ask for too many combinations, magma gives up
// even if it didn't give up, it would use a tremendous amount of memory
BigNumberLimit:=10^8;

found:=false;

expressions:=[];
expression_choices:=[];

TopDim:=Index([Binomial(#vecs,t) lt BigNumberLimit : t in [1..#vecs]],false)-1;
if TopDim eq -1 then
  TopDim:=#vecs;
end if;

print "checking up to dimension", TopDim;
print "total vectors:",#vecs;

for d in [2..TopDim] do

  RR:=RSpace(CoefficientRing(vecs[1]),Degree(vecs[1]));

  choices:=Binom(#vecs,d);

  for u in choices do
    ss:=sub<RR|[vecs[vv] : vv in u]>;
    if v in ss then
      sol:=Solution(ChangeRing(Matrix([vecs[vv] : vv in u]),Rationals()),ChangeRing(v,Rationals()));
      Append(~expressions,sol);
      Append(~expression_choices,u);
      found:=true;
    end if;
  end for;

  if found then
    break d;
  end if;
  d;

end for;

// sorts all the expressions by height
heights:=[&+[HH(v) : v in ElementToSequence(u)]  : u in expressions];
heights_again:=heights;
ParallelSort(~heights,~expressions);
ParallelSort(~heights_again,~expression_choices);

// prints the expression of smallest height to the terminal
if #expressions ne 0 then
  e:=ElementToSequence(expressions[1]);
  for i in [1..#e] do
    print quots[expression_choices[1][i]];
    print expressions[1][i];
    "";
  end for;
end if;

// this returns vitality to magma
delete choices;





/*

// in case you want to see all the expressions found

for j in [1..#expressions] do
  e:=ElementToSequence(expressions[j]);
  for i in [1..#e] do
    print quots[expression_choices[j][i]];
    print expressions[j][i];
    "";
  end for;
  "";
end for;

*/

//
