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
b1 := x; b3 := y; a3 := z; c3 := w;

# We solve the coefficients a1 and a2 so that the polynomial is in Q,
f2 := NormalForm(fGeneric, [mp], plex(a1, a2, x, y, z, w, t)); 
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
# f := 8*w^4+32*w^2*x^2+16*w^2*y*z+8*w^2*z^2+64*w*x^2*y+16*w*x^2*z+8*w*y^2*z+40*x^4+8*x^2*y^2+32*x^2*y*z+16*x^2*z^2+2*y^4+8*y^2*z^2

# We verify that the polynomial is absolutely irreducible
evala(AIrreduc(f));
## -> true

# Matrix Q associated to the problem (parametrization of the space L)
Q, QVars, v := polyToMatrix(f):

# Matrix associated to the original decomposition (for verifications)
MNEW := decompositionToMatrix([p1ss, p2ss, p3ss], v):

#-----------------------------------

# We start from Q and go step by step.
nops(indets(Q));  # 20
randomRank(Q);    # 10

# Real solutions
sSym := solve({f=0, diff(f, x)=0, diff(f,y)=0, diff(f,z)=0, diff(f,w)=0});
# Alternative: 
# sSym := solve({p1ss = 0, p2ss = 0, p3ss = 0});

# We verify that the third branch contains real points:
evalf(allvalues(RootOf(27*_Z^12+108*_Z^10+36*_Z^8-2272*_Z^6-48*_Z^4+192*_Z^2-64)));
## -> 1.817516153, .4774054007+.1915725139*I, .8359620401+2.083246622*I, .6353178959*I, -.8359620401+2.083246622*I, -.4774054007+.1915725139*I, -1.817516153, -.4774054007-.1915725139*I, -.8359620401-2.083246622*I, -.6353178959*I, .8359620401-2.083246622*I, .4774054007-.1915725139*I

## sSym[1] plain equations - reduction to 14 variables and rank 9
v1 := eval(Vector(v), sSym[1]):
v11 := eval(v1, {z=1}):
simplify(LinearAlgebra[Transpose](v11).MNEW.v11);  # Verification
Q1 := reduceByLinearEquation(Q, v11):
nops(indets(Q1)); # 14
randomRank(Q1); # 9

## sSym[2] plain equations - reduction to 9 variables and rank 8
v2 := eval(Vector(v), sSym[2]):
v21 := eval(v2, {w=1}):
evalf(allvalues(v21));
simplify(LinearAlgebra[Transpose](v21).MNEW.v21);  # Verification
Q2 := reduceByLinearEquation(Q1, v21):
nops(indets(Q2)); # 9
randomRank(Q2); #8

## sSym[3] trace equations - reduction to 1 variable and rank 7
v3 := eval(Vector(v), sSym[3]):
v31 := eval(v3, {x=1}):
v31t := vectorTrace(v31);
Q3 := reduceByLinearEquation(Q2, v31t):
nops(indets(Q3)); # 1
randomRank(Q3); # 7

## Matrix Q3 cannot be positive semidefinite:

Q3[2,2];  # 0
Q3[2,6];  # 32

# The determinant of submatrix using rows and columns 2 and 6 is negative.

# This finishes the proof that there does not exist a rational SOS for f.

################################################################################

## If we use the trace of the second solutions, we get a further reduction:

v2 := eval(Vector(v), sSym[2]):
v21 := eval(v2, {w=1}):
v21t := vectorTrace(v21);
Q2 := reduceByLinearEquation(Q1, v21t):
nops(indets(Q2));
randomRank(Q2);

v3 := eval(Vector(v), sSym[3]):
v31 := eval(v3, {x=1}):
v31t := vectorTrace(v31);
Q3 := reduceByLinearEquation(Q2, v31t):
nops(indets(Q3));
randomRank(Q3);

# In this case the matrix found is unique and it is not positive semidefinite:
LinearAlgebra[Eigenvalues](Q3);
# Vector(10, {(1) = 0, (2) = 0, (3) = 0, (4) = 0, (5) = 0, (6) = -32, (7) = 16, (8) = 32, (9) = 33-sqrt(481), (10) = 33+sqrt(481)})


