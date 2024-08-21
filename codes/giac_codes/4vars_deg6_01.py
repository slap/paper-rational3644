from giacpy import giac,linsolve
from sos_decomp import *
from three_sos import *

giac("printpow(1)")

#substitution of arbitrary polynomials for : a3,b1,b2,b3,c1,c2,c3. In the L. construction of 3 sos 
def assign_abc(fcn):
  global X,a,b,c
  p,q=giac("p,q")
  fcn=fcn.subst(b,[X[0]-p*X[1],3*X[0]+p*X[3],X[1]])
  fcn=fcn.subst(c,[X[2]+p*X[3],X[0]+7*X[2]+p*X[1],X[3]])
  fcn=fcn.subst(a[2],21*X[2])
  return fcn.ratnormal() #verified this to be the same as in the paper

a=giac("[a1,a2,a3]")
b=giac("[b1,b2,b3]")
c=giac("[c1,c2,c3]")

X=[]
n=4 #no. of variables 
for i in xrange(n):
  s="x"+str(i)
  X.append(giac(s))
d=3
lst = lexico_fixed_sum(len(X),d)

t = giac("t") #this is also the primitive root
f = t**3 -2
t0 = f.proot(t,32)[2] #float approx of real root of f

F,p1,p2,p3 = create_3_sos(a,b,c,t,f)

F=assign_abc(F)
F= F.numer().ratnormal() #verified this to be the same as in the paper
print F
"""
121*x0^6+132*x0^5*x2+42*x0^4*x1^2-32*x0^4*x1*x2+16*x0^4*x1*x3+1627*x0^4*x2^2+84*x0^4*x2*x3-32*x0^4*x3^2+144*x0^3*x1^2*x2-48*x0^3*x1*x2^2-160*x0^3*x1*x2*x3-2984*x0^3*x2^3+288*x0^3*x2^2*x3-48*x0^3*x2*x3^2+5*x0^2*x1^4-8*x0^2*x1^3*x2-108*x0^2*x1^2*x2^2+20*x0^2*x1^2*x2*x3-8*x0^2*x1^2*x3^2+564*x0^2*x1*x2^3+352*x0^2*x1*x2^2*x3+8*x0^2*x1*x2*x3^2+4751*x0^2*x2^4-224*x0^2*x2^3*x3+584*x0^2*x2^2*x3^2-16*x0^2*x2*x3^3+4*x0^2*x3^4+4*x0*x1^4*x2+16*x0*x1^3*x2^2+80*x0*x1^2*x2^3+16*x0*x1^2*x2^2*x3+16*x0*x1^2*x2*x3^2-896*x0*x1*x2^4+192*x0*x1*x2^3*x3-48*x0*x1*x2^2*x3^2-7124*x0*x2^5+208*x0*x2^4*x3-880*x0*x2^3*x3^2+32*x0*x2^2*x3^3-24*x0*x2*x3^4+x1^4*x2^2+12*x1^3*x2^3+162*x1^2*x2^4+4*x1^2*x2^3*x3+12*x1^2*x2^2*x3^2+812*x1*x2^5+40*x1*x2^4*x3+80*x1*x2^3*x3^2+4077*x2^6+244*x2^5*x3+816*x2^4*x3^2+24*x2^3*x3^3+40*x2^2*x3^4
"""

#below are used just to solve 0's so it is ok if we mangle them a bit (numer)
p1=assign_abc(p1)
p2=assign_abc(p2)
p3=assign_abc(p3)

#Groebner Write: below is the code to get the Groebner basis, uncomment to generate (this is saved in gb_5_1_x0.txt)
"""
from giacpy import gbasis
vars = [X[1],X[2],X[3],t]
ideal = giac([p1,p2,p3,f])
gb = gbasis(ideal.subst(X[0],1),vars,"plex")
file = open("./temp.txt", "w")
file.write(str(gb))
file.close()
"""
#End of Groebner Write

P=[p1.subst(t,t0),p2.subst(t,t0),p3.subst(t,t0)]
G0 = get_psd(P,X,lst)
#print "One Gram matrix is: ", G0

print "The polynomials: "
print p1
print p2
print p3
"""
-p^3*x2^3-9*p^2*x0*x2^2+2*p*t^2*x2^3+2*p*t*x0*x2^2-28*p*x0^2*x2-p*x1^2*x2-42*p*x2^2*x3-2*t^2*x0^2*x2-8*t^2*x0*x2^2-2*t*x0^3-8*t*x0^2*x2-28*x0^3+42*x0^2*x2-3*x0*x1^2+42*x0*x1*x2+296*x0*x2^2-126*x0*x2*x3+2*x0*x3^2+294*x1*x2^2+700*x2^3+14*x2*x3^2
2*p^2*t*x2^3+p^2*x0*x2^2+2*p*t^2*x0*x2^2+14*p*t^2*x2^3-2*p*t*x0^2*x2-2*p*t*x0*x2^2+6*p*x0^2*x2-2*t^2*x0^3-22*t^2*x0^2*x2-56*t^2*x0*x2^2-6*t*x0^3-24*t*x0^2*x2+10*x0^3-2*x0^2*x2+x0*x1^2-28*x0*x2^2+42*x0*x2*x3-42*x1*x2^2-100*x2^3-2*x2*x3^2
2*p*t^2*x2^2*x3+2*p*t*x1*x2^2+42*p*x2^3-2*t^2*x0^2*x3-8*t^2*x0*x2*x3-2*t*x0^2*x1-8*t*x0*x1*x2-42*x0^2*x2-168*x0*x2^2
#p=1 works!
"""

print "Factors of p1,p2,p3.."
print p1.factors()
print p2.factors()
print p3.factors()
"""
[2*t^2*x2^3*p-8*t^2*x2^2*x0-2*t^2*x2*x0^2+2*t*x2^2*x0*p-8*t*x2*x0^2-2*t*x0^3-x2^3*p^3+700*x2^3-9*x2^2*x0*p^2+296*x2^2*x0-42*x2^2*p*x3+294*x2^2*x1-28*x2*x0^2*p+42*x2*x0^2+42*x2*x0*x1-126*x2*x0*x3-x2*p*x1^2+14*x2*x3^2-28*x0^3-3*x0*x1^2+2*x0*x3^2,1]
[2*p^2*t*x2^3+p^2*x2^2*x0+14*p*t^2*x2^3+2*p*t^2*x2^2*x0-2*p*t*x2^2*x0-2*p*t*x2*x0^2+6*p*x2*x0^2-56*t^2*x2^2*x0-22*t^2*x2*x0^2-2*t^2*x0^3-24*t*x2*x0^2-6*t*x0^3-100*x2^3-28*x2^2*x0-42*x2^2*x1-2*x2*x0^2+42*x2*x0*x3-2*x2*x3^2+10*x0^3+x0*x1^2,1]
[2,1,21*x2+x3*t^2+x1*t,1,-p*x2^2+4*x2*x0+x0^2,1,-1,1]
"""

exit(0)

print "Factors of p1,p2,p3 in the hyperplane x=0"
print p1.subst(X[0],0).factors()
print p2.subst(X[0],0).factors()
print p3.subst(X[0],0).factors()

print "Factors of p1,p2,p3 in the hyperplane x+z=0"
print p1.subst(X[0],5*X[2]).factors()
print p2.subst(X[0],5*X[2]).factors()
print p3.subst(X[0],5*X[2]).factors()
"""
[2,1,x2,1,5*x1^2-2*x1*x2+615*x2^2+10*x2*x3-2*x3^2,1,-1,1]
[x2,1,5*x1^2-2*x1*x2+615*x2^2+10*x2*x3-2*x3^2,1]
[0,1]
"""
#STOP: from the factors we should have enough real points to continue (to easily conclude non-Q-sos), otherwise retry

nlst,b,sol = gram_solve(2*d,F,X,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b)
#126 indets for G, just as in the paper!

mvec = mon_vectors(X,lst)
mvec = giac(mvec)

###procedure to find as many element as possible in the kernel of G
"""
Note: A Groebner basis solution should be first made before coming to this part of the program!
We did this part in a separate file, because this could take time (although for giac this was rather fast!)
The results is in gb_5_1.txt
"""

vec = mvec.subst(X,[0,1,0,1])
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
#110 indets for G, just as in the paper!

vec = mvec.subst(X,[0,1,0,0]) #still in the branch, 
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
#95 indets for G, just as in the paper!

vec = mvec.subst(X,[0,0,0,1]) #still in the branch, 
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
#81 indets for G, just as in the paper!

vec = mvec.subst(X,[0,-1,0,1]) #still in the branch, 
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
#71 indets for G, just as in the paper!

#dimension reduction from zero-diagonals
G=G.ratnormal() #need to ratnormal to get out the 0's, some equations are just artifacts but actually 0
G,b = reduce_zero_diag(G,b)
print "Dimension:" , len(b)
#here the dimension is only reduced to 59
print "Trying again.."
G,b = reduce_zero_diag(G,b)
print "Dimension:" , len(b)
#here we get the same as in paper: 48 indets for G
#print "Debug, parametrized vs first Gram: ", G.subst(t,t0)-G0

#dimension reduction using principal submatrix
G,b = reduce_zero_psubmatrix(G,b)
print "Dimension:" , len(b)
#here we get the same as in paper: 29 indets for G

#f1 = -21*t-168*X[3]+4*(t**2)+8*(X[3]**2)+165

#if X[3]=1 f1(t,3) = (4*t-1)*(t-5) => root t=1/4 , t=5
#by result of gb_5_1.txt, X[3]=1 (real root for f1) =>  X=[1,(1/2)t,-1/4,1]=[1,1/8,-1/4,1]
vec = mvec.subst(X,giac("[1,1/8,-1/5,1]")) #still in the branch, 
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets()
print "Dimension:" , len(b)
#21 indets for G, just as in the paper!
exit(0) #jcapco now

#introduced alg. number t1
#by result of gb_5_1.txt, X[3]=2 (real root for f1), we are now dealing with algebraic numbers 
f = f1.subst(X[3],2).ratnormal()
vec = mvec.subst(X,giac("[1,(1/2)*t,-1/4,2]")) #still in the branch
for i in xrange(len(vec)):
  vec[i] = vec[i].rem(f,t)
sol = linsolve((G*vec).ratnormal(), b)
sol = sol.subst(t,"t1")
G = G.subst(b,sol)
b=G.indets().remove("t1") #remove t1
ft1 = f.subst(t,"t1") #poly for t1
print "Dimension:" , len(b)
#12 indets for G, NOT as in the paper (2 less), but maybe this will improve
#t1s = f.proot(t,32)
#print "Debug, parametrized vs first Gram: ", (G.subst("[t,t1]",[t0,t1s[0]])-G0).ratnormal()
#check_constants(G.subst("[t,t1]",[t0,t1s[0]])-G0)

#introduced alg. number t2
#X[3]=3 
f = f1.subst(X[3],3).ratnormal()
vec = mvec.subst(X,giac("[1,(1/2)*t,-1/4,3]")) #still in the branch
for i in xrange(len(vec)):
  vec[i] = vec[i].rem(f,t)
sol = linsolve((G*vec).ratnormal(), b)
sol = sol.subst(t,"t2")
G = G.subst(b,sol)
b=G.indets().remove("t1").remove("t2")
ft2 = f.subst(t,"t2") #poly for t2
print "Dimension:" , len(b)
#6 indets ofr G, NOT as in the paper (2 less)
#t2s = f.proot(t,32)
#print "Debug, parametrized vs first Gram: "
#check_constants(G.subst("[t,t1,t2]",[t0,t1s[0],t2s[0]])-G0)

## Fourth branch when X[0]=0
#x1 = -(1/21)*(50*x2^2+x3^2)/x2, x2 = x2, x3 = x3 => [0,-51/21,1,1] (X[2]=X[3]=1)
vec = mvec.subst(X,giac("[0,-51/21,1,1]"))
sol = linsolve((G*vec).ratnormal(), b)
G = G.subst(b,sol)
b=G.indets().remove("t1").remove("t2")
print "Dimension:" , len(b)

G,b = reduce_zero_diag(G,b)
b = b.remove("t1").remove("t2")
print "Dimension:" , len(b)
#3 indets, just as in the paper!!
#print "Debug, parametrized vs first Gram: ", (G.subst("[t,t1]",[t0,t1s[0]])-G0).ratnormal()
#check_constants(G.subst("[t,t1]",[t0,t1s[0]])-G0)

t1 = giac("t1")
t2 = giac("t2")
ft0 = t**3 -2
#recall that ft1 and ft2 were already defined

#Groebner Get: comment this section if you uncommented the generation of Groebner basis 
"""
file = open("./gb_5_1_x0.txt", "r")
data = file.readline()
file.close()
gb = giac(data)
"""
#end of Groebner Get

#upon inspection we just need only a few of the equation from the first factor of gb[3]
fx3 = gb[3].factors()[0]
#print "Degree in x3: ", fx3.degree(X[3]) #degree 4
fx2 = gb[2].factors()[2] #gb[2] is linear wrt to x2
x2 = -fx2.coeff("x2",0)/fx2.coeff("x2",1)
fx1 = gb[1].factors()[2].subst("x2",x2).ratnormal() #gb[1] is linear wrt to x1
x1 = -fx1.coeff("x1",0)/fx1.coeff("x1",1)

vec = mvec.subst("[x0,x1,x2]",[1,x1,x2])
eqns = (G*vec).ratnormal() #G has non-zero denominators in t1 and t2 , so we need to get rid of them before

eqns2=[]
for e in eqns:
  temp = e.numer().rem(fx3,X[3]).rem(ft0,t).rem(ft1,t1).rem(ft2,t2)
  if not temp==0:
    eqns2.append(temp)

#we cannot seem to handle this final set of equations in Giac so we will save it and pass it to another CAS
file = open("./eqns_5_1.txt", "w")
file.write(str(eqns2))
file.close()
#maple result: {t1 = t1, t2 = t2, x3 = x3, b_2_8 = -5376-336*t^2, b_7_9 = 592-16*t^2, b_8_8 = 17640+128*t}