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
com = gcd(p1,p2,p3)
p1 = (p1/com).ratnormal()
p2 = (p2/com).ratnormal()
p3 = (p3/com).ratnormal()

print "The polynomials: "
print p1
print p2
print p3
"""
-4*t^2*x0^2+4*t*x0^2-2*x0^2+x1^2-2*x1*x2+2*x2*x3-2*x3^2
-4*t^2*x0^2-4*t*x0^2+6*x0^2+x1^2+2*x1*x2+2*x2*x3+2*x3^2
4*t^2*x0*x3+4*t*x0*x1+4*x0*x2
"""

F=assign_abc(F)
F = F/(com**2)
F= F.ratnormal()
print "The result sos: "
print F
print "Number of monomials: " , len(F) #13
#40*x0^4+8*x0^2*x1^2+32*x0^2*x1*x2+64*x0^2*x1*x3+16*x0^2*x2^2+16*x0^2*x2*x3+32*x0^2*x3^2+2*x1^4+8*x1^2*x2^2+8*x1^2*x2*x3+16*x1*x2*x3^2+8*x2^2*x3^2+8*x3^4

P=[p1.subst(t,t0),p2.subst(t,t0),p3.subst(t,t0)]
G0 = get_psd(P,X,lst)
print "One Gram matrix is: ", G0
"""
from giacpy import gbasis
vars = [X[1],X[2],X[3],t]
ideal = giac([p1,p2,p3,f0])
gb = gbasis(ideal.subst(X[0],1),vars,"plex")
file = open("./gb_n4_d4_x0.txt", "w")
file.write(str(gb))
file.close()
"""

print "Factors of p1,p2,p3.."
print p1.factors()
print p2.factors()
print p3.factors()
"""
[-2*x2*x1+2*x2*x3-4*x0^2*t^2+4*x0^2*t-2*x0^2+x1^2-2*x3^2,1]
[2*x2*x1+2*x2*x3-4*x0^2*t^2-4*x0^2*t+6*x0^2+x1^2+2*x3^2,1]
[4,1,x3*t^2+t*x1+x2,1,x0,1]
"""

print "Factors of p1,p2,p3 in the hyperplane x=0"
print p1.subst(X[0],0).factors()
print p2.subst(X[0],0).factors()
print p3.subst(X[0],0).factors()
"""
[-2*x2*x1+2*x2*x3+x1^2-2*x3^2,1]
[2*x2*x1+2*x2*x3+x1^2+2*x3^2,1]
[0,1]
hyperplane section zero dim.
"""
#p1!=p2, p3 == 0, this hyperplane section is zero-dim

#looking at the other hyperplane section 
print p1.subst(X[2],-(X[3]*(t**2)+t*X[1])).rem(f0,t).factors()
print p2.subst(X[2],-(X[3]*(t**2)+t*X[1])).rem(f0,t).factors()
print p3.subst(X[2],-(X[3]*(t**2)+t*X[1])).rem(f0,t).factors()
"""
[-4*x0^2*t^2+4*x0^2*t-2*x0^2+2*x1^2*t+x1^2+2*x1*x3*t^2-2*x1*x3*t-2*x3^2*t^2-2*x3^2,1]
[-4*x0^2*t^2-4*x0^2*t+6*x0^2-2*x1^2*t+x1^2-2*x1*x3*t^2-2*x1*x3*t-2*x3^2*t^2+2*x3^2,1]
[0,1]
hyperplane section zero dim
"""

nlst,b,sol = gram_solve(2*d,F,X,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b)
#20 unknowns

mvec = mon_vectors(X,lst)
mvec = giac(mvec)

#[0,0,1,0] is a solution
#test_zerop(X,[0,0,1,0],p1,p2,p3)
vec = mvec.subst(X,[0,0,1,0])
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
print "(Debug) unknowns are: ", b
print "(Debug) Gram matrix determined as: ", G
#14 indets for G

#[0,2*t,-t**2,2] is a solution
#test_zerop(X,[0,2*t,-(t**2),2],p1,p2,p3,[f0],[t])
vec = mvec.subst(X,[0,2*t,-t**2,2])
eqns = modout_eqns((G*vec).ratnormal(),[f0],[t])
sol = linsolve(eqns, b)
#here linsolve is possibly missing solution

G = G.subst(b,sol).ratnormal()
b=G.indets().remove(t)
print "Dimension:", len(b)
print "(Debug) unknowns are: ", b
print "(Debug) Gram matrix determined as: ", G
#7 indets for G

G,b = reduce_zero_diag(G,b)
b = b.remove(t)
print "Dimension:" , len(b)
print "(Debug) unknowns are: ", b
print "(Debug) Gram matrix determined as: ", G
#3 indets for G

G,b = reduce_zero_psubmatrix(G,b)
b = b.remove(t)
print "Dimension:" , len(b)
print "(Debug) unknowns are: ", b
print "(Debug) Gram matrix determined as: ", G
#2 indets for G

"""
#uncomment this for indirect proof
#We take tower of fields Q<Q(t)<Q(t,t1), Q(t,t1)/Q(t) has minpoly ft1, Q(t)/Q has minpoly f0
for i in xrange(len(vec)):
  vec[i] = trace(NM(vec[i],ft1,t1))
  vec[i] = trace(NM(vec[i],f0,t)) #trace is transitive
#print vec #[192,0,0,0,-128,-256,0,0,-128,-128]
#print "Debug: ", G0*vec #is 0, but this is because we are proving by contradiction
eqns = modout_eqns((G*vec).ratnormal(),[f0],[t])
sol = linsolve(eqns, b)
G = G.subst(b,sol).ratnormal()
b = G.indets()
print "Dimension: ", len(b)
print G
print "This is not positive semidefinite!"
#there are two indeterminates and we are pretty sure these are non-rational algebraic numbers -- todo
#Above shows there is no Q Gram matrix, uncomment this and comment below to see 
"""

#below proves directly uniqueness, comment this if you want to use indirect proof
s = giac("s")
fs = -4+3*s**4+(-4*t+4)*s**2
#[4,t1*t*(-3*t1^2+4*t-6), t1*(3*t1^2*t^2+2*t^2-8),4*t1] is a solution
pz = [4,s*t*(-3*(s**2)+4*t-6), s*(3*(s**2)*(t**2)+2*(t**2)-8),4*s]
#test_zerop(X,pz,p1,p2,p3,[fs,f0],[s,t]) #OK
vec = mvec.subst(X,pz); #[0,2*t,-t**2,2])
eqns = modout_eqns((G*vec).ratnormal(),[fs,f0],[s,t])
print "(Debug) last system of equation is now: ", eqns
#sol = linsolve(giac(eqns), b) #giac is pretty stupid in solving linear equations

"""
print "We now need to solve the following ", eqns
print "where t,s are algebraic numbers that satisfy":
print f0
print fs
print "respectively\n"
"""
print "Unknowns: " , b
eqns.append(fs); eqns.append(f0)
vars = giac("[b_0_8,b_0_9,s,t]")
gb = gbasis(eqns,vars,"plex")
#print linsolve(eqns, b)
print gb
#[b_0_8+16*t^2-8,b_0_9+16*t-16,3*s^4-4*s^2*t+4*s^2-4,t^3-2]
