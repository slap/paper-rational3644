from giacpy import giac, rem, coeff, idn

"""
Let Q<= K be a finite field extension of the rationals Q. Also suppose that K=Q(x) for some algebraic number x defined by the zero of an irred. poly f 
over Q. K is then a finite Q v.s. of dimension d=deg(f), with canonical basis 1,x,x^2,..,x^(d-1)
Also, for each a in K (not 0), there is a v.s. automorphism (linear transform) K->K given by b|-> a.b which is associated to a dxd matrix with entries in Q.
The method below returns this matrix (given, a,f and x). The determinant of this matrix is the K/Q norm of a. The trace of this matrix is the K/Q trace of a.
"""
def NM(a,f,x):
  n=f.degree(x)
  #col=i represents vector rep. of a*x^i
  #mat = newMat(n,n) #problem with shallow copy! so we use idn for a moment
  mat = idn(n)
  for i in xrange(0,n):
    g = rem(a*(x**i), f, x)
    for j in xrange(0,n):
      mat[j,i] = coeff(g,x,j)
  return mat

"""
x=giac("x")
y=giac("[x0,x1,x2]")
f=giac("x^4+x^2*t^2+t")
Nm=det(y[0]*NM(1,f)+y[1]*NM(x,f)+y[2]*NM(x**2,f))
#as expected (see proof by Flanders) the above is 
print factor(Nm)
"""
