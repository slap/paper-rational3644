randi := proc(n)
  # random integer
  rand(-n..n)();
end proc:

dg := 3:
alias(t=RootOf(Z^dg-2)):
randomize(1234): #seed to reproduce the results


p := seq(add(a[i,j]*t^(j-1),j=1..dg), i=1..dg):  
F := add(p[i]^2, i=1..dg):
F:= evala(expand(F)):  
eqns := [seq(coeff(F,t,i),i=1..dg-1)]:
sub := [seq(a[i,1]=0,i=2..dg)]:
c := subs(sub, eqns):

sol := solve(eqns,[seq(a[i,1],i=2..dg)])[1]:
p := seq(subs(sol,p[i]),i=1..dg):

a:=Matrix(dg,dg):
a[2..dg,2..dg] := Matrix([[x,z+y],[z-y,x]]):
adj := linalg:-adj(%):
adj := linalg:-transpose(%):  
#force the quadratic denominator to be the polynomial defining a circle
adj := convert(adj,list):
q := linalg:-det(a[2..dg,2..dg]):

k:= 0:
#parametrize the linear forms of the first row by unknowns b1 to b9
for i from 1 to dg do
  a[1,i]:=cat(b,(k+1))*x+cat(b,(k+2))*y+cat(b,(k+3))*z:
  k:=k+3:
od:
  
while true do:
  pts := NULL:
  for i from 1 to 6 do  
    tn := randi(5):
    den := 1+tn^2:
    pts := pts, [(1-tn^2)/den,2*tn/den,1]: #rational points on a circle
  od:

  eqns := NULL:
  for P in pts do:
    cP := subs([x=P[1],y=P[2],z=P[3]],c):
    adj1 := subs([x=P[1],y=P[2],z=P[3]],adj[1]):
    eqns := eqns, <cP>^%T.<adj1>:
  od:
  #6 equations, 9 unknowns
  eqns := [eqns]:
  select(tmp->tmp<> 0, %):
  if nops(%)>5 then break; fi:
od:

b7 := randi(3): b8:= randi(3): b9:=randi(3): 

#increase the accuracy to convert to rational, if necessary 
fsol := fsolve(eqns):
sol := convert(fsol,rational, 4):
pNew := seq(factor(subs(sol2,p[i])*q), i=1..dg):

# This is faster than using Groebner bases (Maple 15)
evalf(allvalues(solve({pNew})));

#check if it has max (5) real zeros
#gb := Groebner:-Basis([pNew,z-1],plex(x,y,z)):
#fsolve(gb[2]); 
#printf("%a\n",[p]);

# The resulting polynomial
pp := pNew[1]^2+pNew[2]^2+pNew[3]^2;
pp := simplify(pp);


# Verification with rationalSOS and SEDUMI

# Uncomment and set the path to rationalSOS.mpl file
#currentdir("C:/Users/Usuario/rationalSOS");

#######################################################################
# Load "Rational SOS" procedures
#######################################################################
read("rationalSOS.mpl");
with(rationalSOS);

# Display tables of any size
interface(rtablesize=infinity);

out := exactSOS(pp, realPolynomials = [pNew], facial = "yes", traceEquations = "no", objFunction = "eig", forceRational = "yes"):
eig(out[3]);

# The input for SEDUMI has rank 5 and 2 indeterminates.
# Sedumi finds a positive definite solution, indicating that a rational decomposition exists.
