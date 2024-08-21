#todo collect denominators when reducing
from giacpy import giac,linsolve, gcd, trace, gbasis
from sos_decomp import *
from three_sos import *
from norm import *

giac("printpow(1)")

#substitution of arbitrary polynomials for : a3,b1,b2,b3,c1,c2,c3. In the L. construction of 3 sos 
def assign_abc(fcn):
  global X,a,b,c
  fcn=fcn.subst(b,[X[0],-X[0],X[1]])  
  fcn=fcn.subst(c,[-X[0],-X[0],X[3]])
  fcn=fcn.subst(a[2],X[2])
  return fcn.ratnormal()

a=giac("[a1,a2,a3]")
b=giac("[b1,b2,b3]")
c=giac("[c1,c2,c3]")

X=[]
n=4 #no. of variables 
for i in xrange(n):
  s="x"+str(i)
  X.append(giac(s))
d=2
lst = lexico_fixed_sum(len(X),d)

t,s = giac("t,s") #this is also the primitive root
f = t**3 -s

F,p1,p2,p3 = create_3_sos(a,b,c,t,f)

print "Debuging the polynomials: "
print p1.factors() #cubic
print p2.factors() #cubic
print p3.factors() #conic and line
#7 indeterminates = P^6, we want to cut with a (2Dim) plane to get many real points
"""
[-2*a3*b2*c3+2*a3*b3*c2-b2^3-b2*b3^2-b2*b1^2+2*b2*b1*c1*t+2*b2*c1^2*t^2+2*c3^2*c2+2*c2^3-2*c2*b1^2*t-2*c2*b1*c1*t^2+2*c2*c1^2,1]
[2*a3*b1*c3-2*a3*b3*c1+b1^3+b1*b3^2+b1*b2^2-2*b1*b2*c2*t-2*b1*c2^2*t^2-2*c3^2*c1-2*c1^3+2*c1*b2^2*t+2*c1*b2*c2*t^2-2*c1*c2^2,1]
[2,1,b2*c1-c2*b1,1,a3+b3*t+t^2*c3,1] #to make a quadratic factor, b1!=b2, c1!=c2 random should work

rem b2*c1-c2*b1:
p1: -b2^3+2*c2^3-b2*b3^2+2*c2*c3^2-b1^2*b2+2*c1^2*c2-2*a3*b2*c3+2*a3*b3*c2
p2: b1^3-2*c1^3+b1*b2^2+b1*b3^2-2*c1*c2^2-2*c1*c3^2+2*a3*b1*c3-2*a3*b3*c1
"""

p1=assign_abc(p1)
p2=assign_abc(p2)
p3=assign_abc(p3)
com = gcd(p1,p2,p3)
p1 = (p1/com).ratnormal()
p2 = (p2/com).ratnormal()
p3 = (p3/com).ratnormal()

print "The polynomials: "
print p1
print p2
print p3
"""
-2*s*x0^2-s*x3^2-4*t^2*x0^2+4*t*x0^2+2*x0^2+x1^2-2*x1*x2+2*x2*x3
2*s*x0^2+s*x3^2-4*t^2*x0^2-4*t*x0^2+2*x0^2+x1^2+2*x1*x2+2*x2*x3
4*t^2*x0*x3+4*t*x0*x1+4*x0*x2
"""

F=assign_abc(F)
F = F/(com**2)
F= F.ratnormal()
print "The result sos: "
print F
print "Number of monomials: " , len(F) #13
"""
8*s^2*x0^4+8*s^2*x0^2*x3^2+2*s^2*x3^4+16*s*x0^2*x1*x2+32*s*x0^2*x1*x3+8*s*x1*x2*x3^2+8*x0^4+8*x0^2*x1^2+16*x0^2*x2^2+16*x0^2*x2*x3+2*x1^4+8*x1^2*x2^2+8*x1^2*x2*x3+8*x2^2*x3^2
"""

"""
in maple:
solve([p1,p2,p3,f],[x0,x1,x2,x3,t])
Solutions of p1,p2,p3:
[x0 = 0, x1 = t*x3, x2 = -(1/2)*x3*(s+t^2)/(t+1), x3 = x3], 
[x0 = 0, x1 = 0, x2 = x2, x3 = 0], 
[x0 = RootOf(2*s*Z^2+s+2*RootOf(Z^2-1+(2*t-2/t)*Z)*t^2)*x3, x1 = RootOf(Z^2-1+(2*t-2/t)*Z)*x3, x2 = -t^2*x3-t*RootOf(Z^2-1+(2*t-2/t)*Z)*x3, x3 = x3, t = t]
"""
#checking the last solution
t1,t2 = giac("t1,t2")
f1=t*t1**2-t+(2*t**2-2)*t1
f2=2*s*t2**2+s+2*t1*t**2
#[t2,t1,-t**2-t*t1,t] must be a zero by the above computation
test_zerop(X,[t2,t1,-t**2-t*t1,1],p1,p2,p3,[f2,f1,f],[t2,t1,t])

"""
Get the evaluation of P1 =[0,1,0,0], P2=[0,1,t,-t**2], P3=[t2,t1,-t**2-t*t1,t]
by the monomial vector
"""
mvec = mon_vectors(X,lst)
mvec = giac(mvec)
P1 = mvec.subst(X,[0,1,0,0]).ratnormal()
P2 = mvec.subst(X,[0,1,t,-t**2]).ratnormal()
P3 = mvec.subst(X,[t2,t1,-t**2-t*t1,t]).ratnormal()
modout_vector(P1,[f2,f1,f],[t2,t1,t])
modout_vector(P2,[f2,f1,f],[t2,t1,t])
modout_vector(P3,[f2,f1,f],[t2,t1,t])

print P1
print P2
print P3

"""
[0,0,0,0,1,0,0,0,0,0]
[0,0,0,0,1,t,-t^2,t^2,-s,s*t]
[t2^2,t1*t2,-t^2*t2-t*t1*t2,t*t2,t1^2,t^2*t1-t-2*t1,t*t1,s*t+t^2+2*t*t1,-s-t^2*t1,t^2]

[t2^2,t1*t2,-t^2*t2-t*t1*t2,t*t2,-,t^2*t1-t-2*t1,t*t1,s*t+t^2+2*t*t1,-,t^2]
"""
