#from sys import path
#path.insert(0, "../") 

from poly_tools import lexico_fixed_sum, mon_vectors
from giacpy import giac, degree, coeff, linsolve, matrix, subst

#special subtraction of vectors, if b-a has a negative entry it returns empty vector, otherwise it returns the difference
def subt(b,a):
  out=[]
  for i in xrange(0,len(b)):
    t = b[i]-a[i]
    if t<0 : return []
    out.append(t)
  return out

"""
let a be a tuple of size n (number of variables) and total even sum d (degree of form f), find two tuples of size
n and total sum d/2 that adds up to a. lst are all tuples with total sum d/2. Outputs are pair of indices of the tuples in the list
preserving the order of the list in the pair.
"""  
def find_half_tuples(a,lst):
  d = (len(lst)/2)    
  out = [] 
  for i in xrange(0,len(lst)):
    temp=subt(a,lst[i])
    if len(temp)>0:      
      #if i<d: left = 1 #less efficient, but ok
      j = lst.index(temp)
      if i>j: continue #this gaurantees unique values in out
      out.append([i,j])
      if i>d and temp==lst[i]: break
  return out

"""
find b's (see power's algorithm)
such that their sum is a and create a linear equation with unknowns as b
"""  
def make_eqn(a,lst,f,vars,b):
  tuples = find_half_tuples(a,lst)
  out = 0
  for t in tuples:
    if t[0]!=t[1]: out = out + 2*b[(t[0],t[1])]
    else: out = out + b[(t[0],t[1])]
  return out - coeff(f,vars,a)

def gram_solve(d,f, vars, lst=None):
  if not lst: 
    lst = lexico_fixed_sum(len(vars),d/2)
  n = len(lst)
  b={} #we use a dictionary because we just want to store the upper half triangle
  for i in xrange(0,n):
    for j in xrange(i,n):
      st="b_"+str(i)+"_"+str(j)
      b[(i,j)]=giac(st)
  lst2 = lexico_fixed_sum(len(vars),d)
  eqns =[]

  for l in lst2:
    #if coeff(f,vars,l)!=0: #this is a bug, 0 equations should be included to get rid of some indets.
    eqns.append(make_eqn(l,lst,f,vars,b))
    
  #print "Basis vectors are: "
  #print mon_vectors(vars,lst)
  print "Length of basis: " , n
  sol= linsolve(eqns,b.values())  
  #remove the 0 contents
  """
  for e in eqns:
    if e != 0: out.append(e)
  """
  return n,b, sol #note: this is not yet the gram matrix.. maybe we need a short way of defining upper-half triangular matrices    
  
def convert_gmatrix(n,b,sol):
  out = matrix(n,n)
  cnt = 0
  for t in b.keys():
    b[t] = sol[cnt]
    cnt = cnt+1
  for i in xrange(n):
    for j in xrange(i,n):
      out[i,j] = b[(i,j)]
      if i!=j: out[j,i] = b[(i,j)]
  return out

#find the psd (Gram matrix) if the sos rep. is already given
def get_psd(P,vars,lst):
  out  = matrix(len(lst),len(lst))
  m = matrix(1,len(lst))
  for p in P:
    for i in xrange(len(lst)):
      m[0,i] = p.coeff(vars,lst[i])
    out = out + (m.transpose())*m
  return out

#test
#f = giac("2*x^6-4*x^5*z+2*x^4*y^2+2*x^4*z^2-10*x^3*y*z^2-6*x^3*z^3+10*x^2*y^3*z+27*x^2*y^2*z^2+6*x^2*y*z^3+6*x^2*z^4-2*x*y^5-18*x*y^4*z-2*x*y^3*z^2-12*x*y^2*z^3+2*x*z^5+3*y^6+3*y^4*z^2-3*y*z^5+5*z^6")
"""
f=giac("10*x^4+6*x^3*y-22*x^3*z+39*x^2*y^2-24*x^2*y*z+33*x^2*z^2-20*x*y^2*z+8*x*y*z^2-20*x*z^3+25*y^4+10*y^3*z+y^2*z^2+4*z^4")
x = giac("x")
y=giac("y")
z=giac("z")

n,b,sol,eqns = gram_matrix(f,[x,y,z])
gM = convert_gmatrix(n,b,sol)
print gM
print eqns
"""