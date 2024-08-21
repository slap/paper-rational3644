# Uncomment and set the path to rationalSOS.mpl file
#currentdir("C:/Users/User/rationalSOS");

#######################################################################
# Load "Rational SOS" procedures
#######################################################################
read("rationalSOS.mpl");
with(rationalSOS);

# Display tables of any size
interface(rtablesize=infinity);

#######################################################################
# Construction of the example
#######################################################################

# We define a polynomial z as the sum of three squares in an algebraic
# extension of degree 3 with generic coefficients.

mp := t^3-2; 
p1 := c1*t^2+b1*t+a1; 
p2 := c2*t^2+b2*t+a2; 
p3 := c3*t^2+b3*t+a3; 

fGeneric := p1^2+p2^2+p3^2; 
fGeneric := expand(fGeneric);

# We choose the parameters so that all terms are multiple by x, and we can reduce the degree.
b2 := -b1;  c2 := b2;  c1 := b2; 
b1 := x0; b3 := x1; a3 := x2; c3 := x3;

# We solve the coefficients a1 and a2 so that the polynomial is in Q,
f2 := NormalForm(fGeneric, [mp], plex(a1, a2, x0, x1, x2, x3, t)); 
f3 := collect(f2, t); 
lf := CoefficientList(f3, t); 
ss := solve({lf[2], lf[3]}, {a1, a2});

# We plug in the solutions found for a1 and a2 and compute the resulting polynomial
ssDen := denom(rhs(ss[1])); 
p1s := simplify(subs(ss, p1)*ssDen); 
p2s := simplify(subs(ss, p2)*ssDen); 
p3s := simplify(subs(ss, p3)*ssDen);

p1ss := subs({t = RootOf(x^3-2)}, p1s); 
p2ss := subs({t = RootOf(x^3-2)}, p2s); 
p3ss := subs({t = RootOf(x^3-2)}, p3s); 

f := simplify(p1ss^2+p2ss^2+p3ss^2);
# f := 40*x0^4+8*x0^2*x1^2+32*x0^2*x1*x2+64*x0^2*x1*x3+16*x0^2*x2^2+16*x0^2*x2*x3+32*x0^2*x3^2+2*x1^4+8*x1^2*x2^2+8*x1^2*x2*x3+16*x1*x2*x3^2+8*x2^2*x3^2+8*x3^4

# We verify that the polynomial is absolutely irreducible
evala(AIrreduc(f));
## -> true

# Matrix Q associated to the problem (parametrization of the space L)
Q, QVars, v := polyToMatrix(f):

# Matrix associated to the original decomposition (for verifications)
MNEW := decompositionToMatrix([p1ss, p2ss, p3ss], v):


# We start from Q and go step by step.
nops(indets(Q));  # 20
randomRank(Q);    # 10

# Real solutions
sSym := solve({p1ss = 0, p2ss = 0, p3ss = 0});
# Alternative: 
#sSym := solve({f=0, diff(f, x0)=0, diff(f,x1)=0, diff(f,x2)=0, diff(f,x3)=0});

# We replace the RootOf(Z^3+2) by -RootOf(Z^3-2) in the second branch,
# so that we have all solutions in terms of that root
alias(r1 = RootOf(_Z^3-2));
s2 := {x0 = 0, x1 = -2*x2*((-r1)^2-(-r1)+1)/((-r1)^2-(-r1)+2), x2 = x2, x3 = (-r1)*x2};

#####################################################

## Plain reduction

## sSym[3] plain equations - reduction to 14 variables and rank 9
v1 := eval(Vector(v), sSym[3]):
v11 := eval(v1, {x2=1}):
Q1 := reduceByLinearEquation(Q, v11):
nops(indets(Q1)); # 14
randomRank(Q1); # 9

## sSym[2] plain equations - reduction to 9 variables and rank 8
v2 := eval(Vector(v), s2):
v21 := eval(v2, {x2=1}):
Q2 := reduceByLinearEquation(Q1, v21):
nops(indets(Q2)); # 9
randomRank(Q2); #8

## sSym[3] plain equations - reduction to 5 variables and rank 7
v3 := eval(Vector(v), sSym[1]):
v31 := eval(v3, {x0=1}):
Q3 := reduceByLinearEquation(Q2, v31):
nops(indets(Q3)); # 5
randomRank(Q3); # 7

## sSym[3]b plain equations - reduction to 2 variables and rank 6
# The first branch contain two real points. Here we take the second real point.
sS3b := eval(sSym[1], {x0=1}):
v3b := eval(Vector(v), sS3b):
v3b := eval(v3b, {x0=-1}):
simplify(LinearAlgebra[Transpose](v3b).MNEW.v3b);  # Verification
Q4 := reduceByLinearEquation(Q3, v3b):
nops(indets(Q4)); # 2
randomRank(Q4); # 6

# Using the real points, we reduce to only two parameters, but we still don't get the uniqueness of the Gram matrix.

# We now look for ghost solutions

# We try to find some principal submatrix Q4S3 of Q4 for which there is a vector v with real coordiantes such that vt . Q4S3 . v = 0.
# This will imply that also vt . Q4 . v = 0 (filling v with zeros in the other entries).
# and hence Q4 . v = 0, a new linear restriction. 
vv:=[3, 4, 7]; # The indexes are chosen so that the system we obtain has real solutions
# I have tried other triplets and this was the first one I found containing real points.

Q4S3 := SubMatrix(Q4, vv,vv):
#indets(Q4S3); # {a_0[4, 6], a_0[7, 9]}
vx := Vector([z0, z1, z2]):
vxt := LinearAlgebra[Transpose](vx):
ff := expand(vxt . Q4S3 . vx):

# We must find the values of (z0, z1, z2) for which ff is 0 for any choice of the 
# two parameters {a_0[4, 6], a_0[7, 9]} of the matrix Q4.
# For this, we compute the coefficient of each of the two parameters and the independent term.

fco := coeffs(ff, indets(Q4)):
# We get a sytem of the polynomial equations in (z0, z1, z2) and we solve it using command solve.
newS := solve(Equate([fco], Vector(nops([fco]))));

# The system has unique solution (up to conjugacy) and one of the conjugated points is real.
# We can use it as a new restriction.

v10 := [0, 0, z0, z1, 0, 0, z2, 0, 0, 0]:
v10s := eval(v10, newS):
v10s1 := eval(v10s, {z1=1}):
v10sv := Vector(v10s1):

# We verify that this vector contains real points
evalf(allvalues(v10sv));

# We verify that indeed vt . Q4 . v = 0
v10svt := LinearAlgebra[Transpose](v10sv):
ff := expand(v10svt . Q4 . v10sv):
simplify(ff);
# 0 -> We found a new vector that must be in the kernel

# We add this restriction
Q5 := reduceByLinearEquation(Q4, v10sv):
nops(indets(Q5)); # 0
randomRank(Q5); # 3

# We proved unique solution!

Q5;    # The matrix found
MNEW;  # The original matrix corresponding to the real SOS decomposition

# The matrix obtained corresponds to the real SOS decomposition.
