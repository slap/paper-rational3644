from itertools import combinations
from itertools import product as it_prod
from giacpy import giac, linsolve, list2mat, matrix

#remove all the zeros in the vector
def trim_zeros(vec):
  out=[]
  for v in vec:
    v=v.ratnormal()
    if v!=0: out.append(v)
  return out

#checking if the constant entries in M are zero
def check_constants(M):
  for i in xrange(M.rowDim()):
    for j in xrange(M.colDim()):
      if len(M[i,j].indets())==0:
        print M[i,j]

def test_zerop(X,P,p1,p2,p3,modF=None,varF=None):
  print "Testing if the summands are zero at ", P
  
  if modF:
    q1 = p1.subst(X,P)
    q2 = p2.subst(X,P)
    q3 = p3.subst(X,P)
    
    for i in xrange(len(modF)):
      q1 = q1.rem(modF[i],varF[i]).ratnormal().numer()
      q2 = q2.rem(modF[i],varF[i]).ratnormal().numer()
      q3 = q3.rem(modF[i],varF[i]).ratnormal().numer()
    
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
      temp = temp.rem(modF[i],varF[i]).ratnormal()
    if not temp==0:
      out.append(temp)
  return out

def modout_vector(v, modF, varF):
  for i in xrange(len(v)):
    for k in xrange(len(modF)):
      v[i]=v[i].rem(modF[k],varF[k]).ratnormal()

#this can be expensive!
def modout_matrix(m, modF, varF):
  out = matrix(m.rowDim(),m.colDim())
  for i in xrange(m.rowDim()):
    for j in xrange(m.colDim()):
      n = m[i,j].numer()
      d = m[i,j].denom()
      for k in xrange(len(modF)):
          n = n.rem(modF[k],varF[k]).ratnormal()
          d = d.rem(modF[k],varF[k]).ratnormal()
      out[i,j] = (n/d).ratnormal()
  return out
  
def reduce_zero_diag(G,b, only_eqn=False):
  zdiag = []
  mat=G
  indets = b
  n = int(mat.rowDim()) #you cannot use directly giac_gen to create a vector n*[0]
  for i in xrange(n):
    if mat[i,i]==0: zdiag.append(i)    
  eqns = []
  for i in zdiag:      
    temp = trim_zeros(mat.row(i))    
    if len(temp)>0:
      eqns = eqns+temp
  if len(eqns)>0:
    if only_eqn:
      print "Only eqns.."
      return eqns
    sol = linsolve(eqns, indets)
    mat = mat.subst(indets,sol)
    indets=mat.indets()
  else: print "Nothing new after zero_diag"
  return mat,indets

#rows and cols is a list of indices
def subMat(m,rows,cols):
  l = []
  for pr in it_prod(rows,cols):
    l.append(m[pr[0],pr[1]])
  return list2mat(l,len(cols))
  
"""
computing kernel of G from the "principal 2x2 submatrices"
here principal submatrices are defined by removing same rows and columns, matrix is diagonal and upper left part 
touches the matrix diagonal.
#todo: for arbitrary submatrices 
"""
def reduce_zero_psubmatrix(G,b,dim=2):
  mat=G
  indets=b
  n=int(mat.rowDim())
  L = range(n)
  eqns = giac([])
  for C in combinations(L,dim): #combination here is ordered!
  #here 2 should be changed if you need more than 2x2 principal submatrix    
    l = []
    for pr in it_prod(C,C):#product here is ordered lexicographically
      l.append(mat[pr[0],pr[1]])
    pmat = list2mat(l,dim) #here 2 is the dim of principal submat
    if pmat.det().ratnormal() == 0:
      k = pmat.ker()[0]
      v = n*[0]
      j=0
      for i in C:
        v[i]=k[j]
        j=j+1
      temp=trim_zeros(G*v)
      if len(temp)>0:
        eqns = eqns.concat(temp)
  if len(eqns)>0:
    print "We got some kernel.."
    sol = linsolve(eqns, indets)
    mat = mat.subst(indets,sol)
    indets=mat.indets()
  return mat,indets

#create three sos assuming field extension defined by degree 3 minpoly with a real root t
#a=coeffs of 0, b=coeffs of t, c=coeffs of t^2
#todo: for pi = deg d we need to write it with coeff up to t**d, automatize to any number
def create_3_sos(a,b,c,t,minpoly):
  p1 = a[0]+b[0]*t+c[0]*(t**2)
  p2 = a[1]+b[1]*t+c[1]*(t**2)
  p3 = a[2]+b[2]*t+c[2]*(t**2)
    
  F = (p1**2 + p2**2 + p3**2).ratnormal()
  F=F.rem(minpoly,t)

  B = F.coeff(t,1)
  C = F.coeff(t,2)
  
  #we can get rid of two variables, a1,a2 and a3 are linear in B,C so we can get rid of one
  sol=linsolve([B,C],a[0:2])
  p1 = p1.subst(a[0:2],sol)
  p2 = p2.subst(a[0:2],sol)
  print "Debug: printing denominator of solution"
  for i in xrange(len(sol)):
    print sol[i].denom()
  den= sol[0].denom()
  p1=p1*den
  p2=p2*den
  p3=p3*den
  F=F.subst(a[0:2],sol)
  F=F.ratnormal().numer()
  #remove common factors
  return F, p1,p2,p3
