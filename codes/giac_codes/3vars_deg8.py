#todo collect denominators when reducing
from giacpy import giac,linsolve, gcd, trace, gbasis
from sos_decomp import *
from three_sos import *
from norm import *

giac("printpow(1)")

#substitution of arbitrary polynomials for : a3,b1,b2,b3,c1,c2,c3. In the L. construction of 3 sos 
def assign_abc(fcn):
  global X,a,b,c
  fcn=fcn.subst(b[0:2],[X[0]*X[1],-X[0]*X[1]])  
  fcn=fcn.subst(c[0:2],[-X[0]*X[1],-X[0]*X[1]])  
  #fcn=fcn.subst(a[2],2*X[1]**2)
  return fcn.ratnormal()

a=giac("[a1,a2,a3]")
b=giac("[b1,b2,b3]")
c=giac("[c1,c2,c3]")

X=[]
n=3 #no. of variables 
for i in xrange(n):
  s="x"+str(i)
  X.append(giac(s))
d=4
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
print com
p1 = (p1/com).ratnormal()
p2 = (p2/com).ratnormal()
p3 = (p3/com).ratnormal()

print "The polynomials: "
print p1
print p2
print p3
"""
-4*t^2*x0^4+4*t*x0^4-2*x0^4+x0^2*x1^2-2*x0*x1^3-2*x0*x1*x2^2+2*x1^3*x2+x2^4
-4*t^2*x0^4-4*t*x0^4+6*x0^4+x0^2*x1^2+2*x0*x1^3-2*x0*x1*x2^2+2*x1^3*x2+x2^4
4*t^2*x0^2*x1*x2+4*t*x0^3*x1-4*t*x0^2*x2^2+4*x0^2*x1^2
"""

print "Factors of p1,p2,p3.."
print p1.factors()
print p2.factors()
print p3.factors()
exit(0)

F=assign_abc(F)
F = F/(com**2)
F= F.ratnormal()
print "The result sos: "
print F
print "Number of monomials: " , len(F) #18
#40*x0^8+8*x0^6*x1^2+32*x0^5*x1^3+64*x0^5*x1^2*x2-16*x0^5*x1*x2^2+18*x0^4*x1^4+16*x0^4*x1^3*x2-64*x0^4*x1*x2^3+8*x0^4*x2^4-8*x0^3*x1^3*x2^2+8*x0^2*x1^6+8*x0^2*x1^5*x2+12*x0^2*x1^2*x2^4-16*x0*x1^4*x2^3-8*x0*x1*x2^6+8*x1^6*x2^2+8*x1^3*x2^5+2*x2^8

P=[p1.subst(t,t0),p2.subst(t,t0),p3.subst(t,t0)]
G0 = get_psd(P,X,lst)
print "One Gram matrix is: ", G0

"""
from giacpy import gbasis
vars = [X[0],X[2],t]
ideal = giac([p1,p2,p3,f0])
gb = gbasis(ideal.subst(X[1],1),vars,"plex")
file = open("./gb_n3_d8_x1.txt", "w")
file.write(str(gb))
file.close()
#[-x0,-x2^4-2*x2,t^3-2]
"""

"""
[-4*t^2*x0^4+4*t*x0^4-2*x0^4+x0^2*x1^2-2*x0*x1^3-2*x0*x1*x2^2+2*x1^3*x2+x2^4,1]
[-4*t^2*x0^4-4*t*x0^4+6*x0^4+x0^2*x1^2+2*x0*x1^3-2*x0*x1*x2^2+2*x1^3*x2+x2^4,1]
[4,1,x0,2,t^2*x1*x2+t*x0*x1-t*x2^2+x1^2,1]
"""

print "Factors of p1,p2,p3 in the hyperplane x=0"
print p1.subst(X[0],0).factors()
print p2.subst(X[0],0).factors()
print p3.subst(X[0],0).factors()
"""
[x2,1,12*x1^3-14*x1^2*x2+x2^3,1]
[x2,1,12*x1^3+14*x1^2*x2+x2^3,1]
[0,1]
# (0,1,0) is a solution!
# (0,-1,t) is a solution!
"""
exit(0)

nlst,b,sol = gram_solve(2*d,F,X,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b)
#75 unknowns!
G,b = reduce_zero_psubmatrix(G,b)
print "Dimension: ", len(b)
#52 unknowns!

mvec = mon_vectors(X,lst)
mvec = giac(mvec)

vec = mvec.subst(X,[1,0,0])
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
exit(0)
#[0,-1,t] is a solution
#test_zerop(X,[0,1,0],p1,p2,p3) #OK
#I'm having problem with this, so I gave it to maple
#vec = mvec.subst(X,[0,-1,t])
#eqns = modout_eqns((G*vec).ratnormal(),[f0],[t])
#print eqns #for maple 
#sol = linsolve(eqns, b) 
sol = "[b_0_11, b_0_12, b_0_14, b_12_14, b_1_11, b_1_14, b_2_11, b_3_11, b_3_8, b_3_9, b_4_11, b_4_12, b_4_7, b_4_8, b_5_11, b_5_12, b_5_6, b_5_7, b_5_9, b_6_12, b_6_8, b_7_13, b_7_7, b_7_9, b_8_12, b_9_11, b_9_12, b_12_14*t^2, 4*t, 4+b_12_14*t, (1/2)*b_5_11*t-(1/2)*b_5_12*t^2, (1/2)*b_6_12*t^2-b_7_13*t, t*b_8_12+b_9_11*t, -4+(1/2)*b_9_11*t-(1/2)*b_9_12*t^2, -(1/2)*b_0_11*t+(1/2)*b_0_12*t^2+b_0_14*t, 6+(1/2)*b_4_11*t-(1/2)*b_4_12*t^2-(1/2)*b_5_11*t^2-b_7_9, -(1/2)*b_1_11*t+b_1_14*t-(1/2)*b_2_11*t^2-(1/2)*b_3_8*t^2-(1/2)*b_4_7*t^2-(1/2)*b_5_6*t^2-2*t^2, -(1/2)*b_3_11*t+6*t-b_4_12-b_5_11-t*b_7_9-(1/2)*b_6_8*t^2-(1/4)*b_7_7*t^2, (1/2)*b_1_11-b_1_14+(1/2)*t*b_3_8+(1/2)*t*b_4_7+(1/2)*t*b_5_6+2*t-(1/2)*b_3_9*t^2-(1/2)*b_4_8*t^2-(1/2)*b_5_7*t^2-b_5_9*t]"
b="[b_0_11, b_0_12, b_0_13, b_0_14, b_11_13, b_11_14, b_12_14, b_1_11, b_1_13, b_1_14, b_2_11, b_2_13, b_3_11, b_3_13, b_3_14, b_3_8, b_3_9, b_4_11, b_4_12, b_4_14, b_4_7, b_4_8, b_5_11, b_5_12, b_5_6, b_5_7, b_5_9, b_6_12, b_6_13, b_6_8, b_7_13, b_7_14, b_7_7, b_7_9, b_8_11, b_8_12, b_8_14, b_9_11, b_9_12]"

G = G.subst(b,sol).ratnormal()
b=G.indets().remove(t)
print "Dimension: ", len(b)
#51 unknowns 

#sol = linsolve(eqns, b)
#print sol
eqns.append(f0)
#vars = giac(G.indets())
gb = gbasis(eqns,G.indets(),"plex")
#print linsolve(eqns, b)
print gb
exit(0)

G = G.subst(b,sol).ratnormal()
b=G.indets().remove(t)
print "Dimension: ", len(b)
#50 unknowns 
#G,b = reduce_zero_diag(G.ratnormal(),b)
#print b
#print "Dimension: ", len(b)
print G
exit(0)
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

