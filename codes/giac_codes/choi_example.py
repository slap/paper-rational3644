from giacpy import giac
from sos_decomp import *
from three_sos import *

#we are working with \sum_{n,2d}
n=3
d=3

giac("printpow(1)")
x=giac("[x,y,z]")

S=giac("x^6+y^6+z^6-(x^4*y^2+x^4*z^2+y^4*x^2+y^4*z^2+z^4*x^2+z^4*y^2)+3*x^2*y^2*z^2")
M=giac("x^6+y^6+z^6")
f = S+M/8
f=2*f #its easier if we deal with twice of f because of choi's sos decomposition

#choi's sos decomposition
p=[]
for i in xrange(n):
  tmp = 3*x[i]**2/2
  for j in xrange(n):
    if j!=i:
      tmp = tmp - x[j]**2
  tmp = x[i]*tmp
  p.append(tmp)

#check if sos of p is the same as f
F = 0
for i in xrange(n):
  F = (F+p[i]**2).ratnormal()

print "We check if CLR decomposition is correct, for 2(S+M/8):"
print F
print f.ratnormal()

"""
p=[]
p.append(giac("x*(x-z)*(x+z)"))
p.append(giac("y*(y-z)*(y+z)"))
p.append(giac("z*(3*x^2+3*y^2-4*z^2)"))
F=(p[0]**2+p[1]**2+p[2]**2).ratnormal()
print F
"""

#Time to get Gram matrix
lst = lexico_fixed_sum(n,d)
G0 = get_psd(p,x,lst)
print "One Gram matrix is : \n", G0 
#rank of G0 is 3, so many vectors are in the kernel!

nlst,b,sol = gram_solve(2*d,F,x,lst)
G=convert_gmatrix(nlst,b,sol)
b=G.indets()
print "Dimension: " , len(b) #27

G,b = reduce_zero_diag(G,b)
print "Dimension:" , len(b)

print "Debug: \n", G