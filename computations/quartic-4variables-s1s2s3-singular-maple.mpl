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

mp := t^3+s1*t^2s2*t+s3; 
p1 := c1*t^2+b1*t+a1; 
p2 := c2*t^2+b2*t+a2; 
p3 := c3*t^2+b3*t+a3; 

fGeneric := p1^2+p2^2+p3^2; 
fGeneric := expand(fGeneric);

# We choose the parameters so that all terms are multiple by x, and we can reduce the degree.
b2 := -b1;  c2 := b2;  c1 := b2; 
b1 := x0; b3 := x1; a3 := x2; c3 := x3;

# We solve the coefficients a1 and a2 so that the polynomial is in Q,
f2 := NormalForm(fGeneric, [mp], plex(a1, a2, x0, x1, x2, x3, t, s1, s2, s3)); 
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

fs := simplify(p1s^2+p2s^2+p3s^2);
fs2 := NormalForm(fs, [mp], plex(t, a1, a2, x0, x1, x2, x3, s1, s2, s3)); 

LIB"primdec.lib";
ring R = (0,s1, s2, s3), (x0, x1, x2, x3), dp;
poly f = 8*s1^4*x0^4+8*s1^4*x0^2*x3^2+2*s1^4*x3^4+8*s1^2*s2^2*x0^4+8*s1^2*s2^2*x0^2*x3^2+2*s1^2*s2^2*x3^4-16*s1^3*x0^2*x1*x3-8*s1^3*x1*x3^3-16*s1^2*s2*x0^4-16*s1^2*s2*x0^2*x3^2-4*s1^2*s2*x3^4-16*s1*s2^2*x0^2*x1*x3-8*s1*s2^2*x1*x3^3-16*s1*s2*s3*x0^4-16*s1*s2*s3*x0^2*x3^2-4*s1*s2*s3*x3^4+16*s1^2*x0^4+8*s1^2*x0^2*x1^2+16*s1^2*x0^2*x2*x3+8*s1^2*x0^2*x3^2+12*s1^2*x1^2*x3^2+8*s1^2*x2*x3^3+16*s1*s2*x0^2*x1*x2+16*s1*s2*x0^2*x1*x3+8*s1*s2*x1*x2*x3^2+8*s1*s2*x1*x3^3+32*s1*s3*x0^4+16*s1*s3*x0^2*x3^2+8*s2^2*x0^4+8*s2^2*x0^2*x3^2+8*s2^2*x1^2*x3^2+2*s2^2*x3^4+16*s2*s3*x0^2*x1*x3+8*s2*s3*x1*x3^3+8*s3^2*x0^4+8*s3^2*x0^2*x3^2+2*s3^2*x3^4-16*s1*x0^2*x1*x3-8*s1*x1^3*x3-16*s1*x1*x2*x3^2-16*s2*x0^4-8*s2*x0^2*x1^2-16*s2*x0^2*x2*x3-8*s2*x0^2*x3^2-16*s2*x1^2*x2*x3-4*s2*x1^2*x3^2-8*s2*x2*x3^3-16*s3*x0^2*x1*x2-32*s3*x0^2*x1*x3-8*s3*x1*x2*x3^2+8*x0^4+8*x0^2*x1^2+16*x0^2*x2^2+16*x0^2*x2*x3+2*x1^4+8*x1^2*x2^2+8*x1^2*x2*x3+8*x2^2*x3^2;
ideal I = f, diff(f, x0), diff(f, x1), diff(f, x2), diff(f, x3);
groebner(I);
radical(I);
list l = minAssGTZ(I);
l;
dim(groebner(l[1]));
dim(groebner(l[2]));
dim(groebner(l[3]));
ideal J = intersect(l[1], l[3]);
groebner(J);
