#todo collect denominators when reducing
from giacpy import giac,linsolve, gcd, trace, gbasis
from sos_decomp import *
from three_sos import *
from norm import *

giac("printpow(1)")

#substitution of arbitrary polynomials for : a3,b1,b2,b3,c1,c2,c3. In the L. construction of 3 sos 
def assign_abc(fcn):
  global X,a,b,c
  fcn=fcn.subst(b,[X[0],X[0],-X[1]])  
  fcn=fcn.subst(c,[X[0],X[0]+3*X[2],X[1]])
  fcn=fcn.subst(a[2],X[1])
  return fcn.ratnormal()

a=giac("[a1,a2,a3]")
b=giac("[b1,b2,b3]")
c=giac("[c1,c2,c3]")

X=[]
n=3 #no. of variables 
for i in xrange(n):
  s="x"+str(i)
  X.append(giac(s))
d=3
lst = lexico_fixed_sum(len(X),d)

t = giac("t") #this is also the primitive root
f = t**3 -2
f0=f
t0 = f.proot(t,32)[2] #float approx of real root of f

F,p1,p2,p3 = create_3_sos(a,b,c,t,f)

print "Debuging the polynomials: "
print p1.factors() #cubic
print p2.factors() #cubic
print p3.factors() #conic and line
#7 indeterminates = P^6, we want to cut with a (2Dim) plane to get many real points
"""
[-2*a3*b2*c3+2*a3*b3*c2-b2^3-b2*b3^2-b2*b1^2+2*b2*b1*c1*t+2*b2*c1^2*t^2+2*c3^2*c2+2*c2^3-2*c2*b1^2*t-2*c2*b1*c1*t^2+2*c2*c1^2,1]
[2*a3*b1*c3-2*a3*b3*c1+b1^3+b1*b3^2+b1*b2^2-2*b1*b2*c2*t-2*b1*c2^2*t^2-2*c3^2*c1-2*c1^3+2*c1*b2^2*t+2*c1*b2*c2*t^2-2*c1*c2^2,1]
[2,1,b2*c1-c2*b1,1,a3+b3*t+t^2*c3,1]
"""

p1=assign_abc(p1)
p2=assign_abc(p2)
p3=assign_abc(p3)

print "The polynomials: "
print p1
print p2
print p3
"""
-6*t^2*x0^2*x2-6*t*x0^2*x2+2*x0^3+24*x0^2*x2-3*x0*x1^2+54*x0*x2^2+54*x2^3
-6*t^2*x0^2*x2-18*t^2*x0*x2^2-6*t*x0^2*x2-2*x0^3-12*x0^2*x2+3*x0*x1^2-18*x0*x2^2
-6*t^2*x0*x1*x2+6*t*x0*x1*x2-6*x0*x1*x2
"""

F=assign_abc(F)
print "The result sos: "
print F
print "number of monomials: ", len(F) #12 yes!
#8*x0^6+144*x0^5*x2-24*x0^4*x1^2+1296*x0^4*x2^2-216*x0^3*x1^2*x2+3672*x0^3*x2^3+18*x0^2*x1^4-540*x0^2*x1^2*x2^2+5832*x0^2*x2^4-324*x0*x1^2*x2^3+5832*x0*x2^5+2916*x2^6

P=[p1.subst(t,t0),p2.subst(t,t0),p3.subst(t,t0)]
G0 = get_psd(P,X,lst)
print "One Gram matrix is: ", G0

"""
from giacpy import gbasis
vars = [X[1],X[2],t]
ideal = giac([p1,p2,p3,f0])
gb = gbasis(ideal.subst(X[0],1),vars,"plex")
file = open("./gb_n3_d6_x0.txt", "w")
file.write(str(gb))
file.close()
#[-3*x1^2+2,-x2,t^3-2]
exit(0)
"""

print "Factors of p1,p2,p3.."
print p1.factors()
print p2.factors()
print p3.factors()
"""
[-6*t^2*x0^2*x2-6*t*x0^2*x2+2*x0^3+24*x0^2*x2+54*x0*x2^2-3*x0*x1^2+54*x2^3,1]
[x0,1,6*t^2*x0*x2+18*t^2*x2^2+6*t*x0*x2+2*x0^2+12*x0*x2+18*x2^2-3*x1^2,1,-1,1]
[6,1,t^2-t+1,1,x2,1,x1,1,x0,1,-1,1]
"""

print "Factors of p1,p2,p3 in the hyperplane x=0"
print p1.subst(X[0],0).factors()
print p2.subst(X[0],0).factors()
print p3.subst(X[0],0).factors()
"""
[54,1,x2,3]
[0,1]
[0,1]
# (0,1,0) is a soltion!
"""

print "Factors of p1,p2,p3 in the hyperplane 3y+z=0"
print p1.subst(X[1],0).factors()
print p2.subst(X[1],0).factors()
print p3.subst(X[1],0).factors()
"""
[2,1,3*t^2*x0^2*x2+3*t*x0^2*x2-x0^3-12*x0^2*x2-27*x0*x2^2-27*x2^3,1,-1,1]
[2,1,x0,1,3*t^2*x0*x2+9*t^2*x2^2+3*t*x0*x2+x0^2+6*x0*x2+9*x2^2,1,-1,1]
[0,1]
"""

nlst,b,sol = gram_solve(2*d,F,X,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b)
#27 unknowns!
G,b = reduce_zero_psubmatrix(G,b)
print "Dimension: ", len(b)
#20 unknowns
G,b = reduce_zero_diag(G,b)
print "Dimension: ", len(b)
#14 unknowns
G,b = reduce_zero_diag(G,b)
print "Dimension: ", len(b)
#9 unknowns

mvec = mon_vectors(X,lst)
mvec = giac(mvec)

#[0,1,0] is a solution
#(+/-sqrt(3/2),1,0) are solutions 
s = giac("s")
fs = 2*s**2-3
#test_zerop(X,[s,1,0],p1,p2,p3,[fs],[s]) #OK
vec = mvec.subst(X,[s,1,0]); 
eqns = modout_eqns((G*vec).ratnormal(),[fs],[s])
sol = linsolve(eqns, b)
G = G.subst(b,sol)
b=G.indets().remove(s)
print "Dimension: ", len(b)
#4 unknowns 

G=G.ratnormal()
G,b = reduce_zero_diag(G,b)
G=G.ratnormal()
b=b.remove(s)
print "Dimension: ", len(b)
#2 unknowns 

vec = mvec.subst(X,[0,1,0]); 
eqns = modout_eqns((G*vec).ratnormal(),[fs],[s])
sol = linsolve(eqns, b)
print sol #this is not good, 2 free variables!

