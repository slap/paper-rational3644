#todo collect denominators when reducing
from giacpy import giac,linsolve
from sos_decomp import *
from three_sos import *

giac("printpow(1)")


def test_zerop(X,P,p1,p2,p3,modF=None,varF=None):
  print "Testing if the summands are zero at ", P
  
  if modF:
    q1 = p1.subst(X,P)
    q2 = p2.subst(X,P)
    q3 = p3.subst(X,P)
    
    for i in xrange(len(modF)):
      q1 = q1.rem(modF[i],varF[i])
      q2 = q2.rem(modF[i],varF[i])
      q3 = q3.rem(modF[i],varF[i])
    
    print q1
    print q2
    print q3
    
  else:
    print p1.subst(X,P).ratnormal()
    print p2.subst(X,P).ratnormal()
    print p3.subst(X,P).ratnormal()

def modout_eqns(eqns, modF, varF):
  out=[]
  for e in eqns:
    temp = e.numer()
    for i in xrange(len(modF)):
      temp = temp.rem(modF[i],varF[i])
    if not temp==0:
      out.append(temp)
  return out

#substitution of arbitrary polynomials for : a3,b1,b2,b3,c1,c2,c3. In the L. construction of 3 sos 
def assign_abc(fcn):
  global X,a,b,c
  #fcn=fcn.subst(b,[X[0],-X[0],X[1]])
  #fcn=fcn.subst(c,[X[0],X[0],3*X[0]+X[3]])
  #fcn=fcn.subst(a[2],X[0]+2*X[1]+5*X[2])
  fcn=fcn.subst(b,[X[0],X[0],3*X[2]])  
  fcn=fcn.subst(c,[X[1],X[2],X[1]])
  fcn=fcn.subst(a[2],-3*X[2])
  return fcn.ratnormal() #verified this to be the same as in the paper

"""
coeffs is a list of triples, coeffs[i] contain the coeffs of t^i for p1, p2 and p3 resp.
"""
def testfcn(coeffs,tosolve,t,minpoly):
  d = len(coeffs)
  p1=0; p2=0; p3=0
  for i in xrange(d):
    p1 = p1+coeffs[i][0]*(t**i)
    p2 = p2+coeffs[i][1]*(t**i)
    p3 = p3+coeffs[i][2]*(t**i)
    
  F = (p1**2 + p2**2 + p3**2).ratnormal()
  F=F.rem(minpoly,t)

  A=[]
  #we want to get rid of coeff of t^i, i>0
  #there are d-1 equations, so we can solve d-1 unknowns, these are in tosolve
  for i in xrange(1,d):
    A.append(F.coeff(t,i))
  sol=linsolve(A,tosolve)
  
  p1 = p1.subst(tosolve,sol)
  p2 = p2.subst(tosolve,sol)
  den= sol[0].denom()
  p1=p1*den
  p2=p2*den
  p3=p3*den
  F=F.subst(tosolve,sol)
  F=F.ratnormal().numer()
  
  return F, p1,p2,p3

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
#F,p1,p2,p3 = testfcn([a,b],[a[0]],t,f)
  
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

P=[p1.subst(t,t0),p2.subst(t,t0),p3.subst(t,t0)]
G0 = get_psd(P,X,lst)
#print "One Gram matrix is: ", G0 #this seems to have 3 eigenvalues so assumedly it should be sum of 3 squares

print "The polynomials: "
print p1
print p2
print p3
"""
-2*t^2*x0^2*x2-8*t^2*x0*x2^2-2*t*x0^3-8*t*x0^2*x2-26*x0^3+38*x0^2*x2-3*x0*x1^2+6*x0*x1*x2+296*x0*x2^2+42*x1*x2^2+700*x2^3,
-2*t^2*x0^3-22*t^2*x0^2*x2-56*t^2*x0*x2^2-6*t*x0^3-24*t*x0^2*x2+10*x0^3+2*x0^2*x2+x0*x1^2-28*x0*x2^2-6*x1*x2^2-100*x2^3,
-2*t^2*x0^3-8*t^2*x0^2*x2-2*t*x0^2*x1-8*t*x0*x1*x2-6*x0^2*x2-24*x0*x2^2

-2*t^2*x0^2*x2-8*t^2*x0*x2^2-2*t*x0^3-8*t*x0^2*x2-26*x0^3-10*x0^2*x2-3*x0*x1^2+22*x0*x1*x2+296*x0*x2^2+154*x1*x2^2+700*x2^3,
-2*t^2*x0^3-22*t^2*x0^2*x2-56*t^2*x0*x2^2-6*t*x0^3-24*t*x0^2*x2+10*x0^3+18*x0^2*x2+x0*x1^2-28*x0*x2^2-22*x1*x2^2-100*x2^3,
-2*t^2*x0^3-8*t^2*x0^2*x2-2*t*x0^2*x1-8*t*x0*x1*x2-22*x0^2*x2-88*x0*x2^2
"""
print "Factors of p1,p2,p3.."
print p1.factors()
print p2.factors()
print p3.factors()

#p3: [2,1,x0,1,x2-x0,1,x1*t+x1+x0*t^2+x0+21*x2,1]
print "Factors of p1,p2,p3 in the hyperplane x=0"
print p1.subst(X[0],0).factors()
print p2.subst(X[0],0).factors()
print p3.subst(X[0],0).factors()
"""
test_zerop(X,"[0,1,0]",p1,p2,p3)
test_zerop(X,"[0,0,1]",p1,p2,p3)
test_zerop(X,"[0,2,1]",p1,p2,p3)
test_zerop(X,"[0,-2,1]",p1,p2,p3)
"""

F=assign_abc(F)
F= F.numer().ratnormal()

nlst,b,sol = gram_solve(2*d,F,X,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b)
#27 unknowns

mvec = mon_vectors(X,lst)
mvec = giac(mvec)

#[0,1,0] is a solution
vec = mvec.subst(X,[0,1,0])

print (G*vec).ratnormal()
print b
sol = linsolve((G*vec).ratnormal(), b)
print G
print sol
exit(0)
G = G.subst(b,sol)

b=G.indets()
print "Dimension:" , len(b)
#20 indets for G
exit(0)
t1 = giac("t1")
ft1 = t1**2+21*t1+5
#test_zerop(X,[0,t1,1],p1,p2,p3,[ft1],[t1])

#[0,t1,1] is a solution
vec = mvec.subst(X,[0,t1,1])
eqns = modout_eqns((G*vec).ratnormal(),[ft1],[t1])
sol = linsolve(eqns, b)
G = G.subst(b,sol)
b=G.indets().remove(t1)
print "Dimension:" , len(b)
#13 indets for G
G,b = reduce_zero_psubmatrix(G,b)
b=b.remove(t1)
print "Dimension:" , len(b)
#9 indets for G

t2 = giac("t2")
ft2 = giac("t2^2+42*t2-32")
#test_zerop(X,[1,t2,1],p1,p2,p3,[ft2],[t2])
#[1,t2,1] is a solution.. no use :(
vec = mvec.subst(X,[1,t2,1])
eqns = modout_eqns((G*vec).ratnormal(),[ft2,ft1],[t2,t1])
sol = linsolve(eqns, b)
G = G.subst(b,sol)
