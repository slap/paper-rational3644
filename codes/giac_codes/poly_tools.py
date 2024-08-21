from random import seed, uniform
from giacpy import giac, re, im, ratnormal
#lexicographic ordering with fixed sum

local_var=giac("lcl")
seed(123) #default seed

"""
create a list of n-dim vector with fixed sum d and order them
this is useful when finding all the monomials (ordered) of a multivariate homogeneous polynomial with n variable and of degree d
"""
def lexico_fixed_sum(n,d):
  x=n*[0]
  x[0]=d
  out=[x[:]]
  k=0
  temp = 0
  while x[n-1]!=d:
    for i in xrange(0,n-1):
      if x[n-1]==d: break
      for j in xrange(n-2,i-1,-1):
        if x[j]!=0:
          if x[n-1]!=0:
            temp=x[n-1]
            x[j+1]=x[j+1]+temp
            x[n-1]=x[n-1]-temp
          x[j]=x[j]-1
          x[j+1]=x[j+1]+1;
          out.append(x[:])
          break;      
  return out
#test
#lst=lexico_fixed_sum(3,4)
#print lst
#print find_half_tuple([2,1,5],lst)

#given a polynomial f (as giac gen type), get a list of non-zero monomials of f
#just type list(f.op())

#test
#from giacpy import giac
#f = giac("2*x*y+3*x**2-4*x**2+y-11*x**2*y")
#print list(f.op())

  
#convert f to a univariate polynomial 
def convert_univariate(f,vars):
  global local_var
  g=f
  for v in vars:
    g=g.subst(v,local_var)
  return g
  
#giac is unable to find the total degree.. we will help
def total_degree(f,vars):
  global local_var
  g = convert_univariate(f,vars)
  return g.degree(local_var)

#random form of total degree d with vars variables 
def random_form(vars,lst, coeff_range = 10):
  #lst = lexico_fixed_sum(len(vars),d)
  out = 0
  for powers in lst:
    r = int(uniform(0,1)*2*coeff_range-coeff_range)
    if abs(r)>coeff_range: continue
    m = r
    for i in xrange(len(vars)):
      m=m*(vars[i]**powers[i])
    out += m
  return out

#test
#vars=giac("[x,y,z]")
#print random_form(vars,6)

#Let vars = [x1,..,xn] and powers = [k1,...,kn]
#then this outputs the monomial : vars^powers , i.e. x1^k1*x2^k2*...*xn^kn
def monomial_from_powers(vars,powers, symb=True):
  out = 1
  for i in xrange(len(vars)):
    if symb : out *= (vars[i]**powers[i])
    else: out *= pow(vars[i],powers[i])
  return out

"""
given a list of list l (say each element of the form [a1,a2,a3]) of size vars (say vars=[x,y,z])
create a vector (output) of monomials with powers given by the list elements (so for the example x^a1*y^a2*z^a3)
"""
def mon_vectors(vars, lst):
  out = []
  for l in lst:
    temp = 1
    for i in xrange(0, len(l)):      
      temp=temp*(vars[i]**l[i])
    out.append(temp)
  return out
  
"""
lst: gives a monomial basis for R[x]_d. So lst is a list of list of integers (consisting of powers of vars
v: is a vector in \R^k where k=Choose(n+d-1,d) (where n = |vars|)
vars: the variable names e.g. [x,y,z,..] as a giac list format
dehomog : If True substitutes 1 for z
Output: a form (or polynomial) that v represents (wrt to basis lst)
"""
def form_from_vector(v,lst,vars):
  out = 0
  for i in xrange(len(v)):
    out += v[i]*monomial_from_powers(vars,lst[i])
  return out

"""
gets a list whose elements are linear complex factors of sum of two squares over reals, say {z_1,z_2,..,z_k}
returns two square reals, say a, such that a^2+b^2 = prod |z_i|
uses the formula: (a^2+b^2)(c^2+d^2)=(ac+bd)^2 + (ad-bc)^2
"""
def two_squares(lst):
  a=re(lst[0])
  b=im(lst[0])
  for i in xrange(1,len(lst)):
    a,b = ratnormal(a*re(lst[i])+b*im(lst[i])), ratnormal(a*im(lst[i])-b*re(lst[i]))
  return a,b

"""
Given two bivariate numeric polynomials , f and g, with 2 variables vars. Find out if they have common factors. 
"""
def is_almost_coprime(f,g,vars, eps=1E-4):
  r=uniform(-1,1)
  f1 = f.subst(vars[0],r)
  f2 = g.subst(vars[0],r)
  return f1.resultant(f2,vars[1]).abs()>eps
