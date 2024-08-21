from sys import path
path.insert(0, "../") 

from giacpy import giac, linsolve
from poly_tools import lexico_fixed_sum, mon_vectors

giac("printpow(1)")

P=[[0,1,1],[1,2,1],[-1,1,1],[-3,-1,1],[5,-1,1],[0,0,1],[5,0,1],[-7,3,1]]
lst = lexico_fixed_sum(3,3)
x,y,z = giac("x,y,z")
mvec = mon_vectors([x,y,z],lst)
mvec = giac(mvec)

coeffs=giac("[-233/5250*C_0-33/350*C_1,313/1050*C_0+23/70*C_1,1331/5250*C_0+181/350*C_1,807/875*C_0+246/175*C_1,-883/5250*C_0-83/350*C_1,-83/525*C_0-8/35*C_1,-C_0-C_1,C_0,C_1,0]")

f1 = mvec*coeffs.subst("[C_0,C_1]",[1,-10])
f1=f1.subst(z,1)
f2 = mvec*coeffs.subst("[C_0,C_1]",[10,-10])
f2=f2.subst(z,1)
print f1
print f2
#here plane curves from f1 and f2 resp. intersect at 9 distinct real points, but neither of them are m-curves