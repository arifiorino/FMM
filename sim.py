import cmath, math, random

p = 100

#Theorem 2.1
def naive_potential(qs,zs,z):
  return sum(qi*cmath.log(z-zi) for qi,zi in zip(qs,zs))
def create_potential_me(qs,zs):
  return (sum(qs), [sum(-qi*pow(zi,k)/k for (qi,zi) in zip(qs,zs)) for k in range(1,p+1)])
def run_potential_me(Q,aks,z):
  return Q*cmath.log(z)+sum(aks[k-1]/pow(z,k) for k in range(1,p+1))


def naive_gravitational(qs,zs,z):
  return sum(qi/(z-zi) for qi,zi in zip(qs,zs))
create_gravitational_me = create_potential_me
def run_gravitational_me(Q,aks,z):
  return Q/z+sum(-k*aks[k-1]/pow(z,k+1) for k in range(1,p+1))


zs=[(1,1),(-1,1),(1,-1),(-1,-1)]
qs=[1,1,1,1]
z=(10,10)

zs=[complex(*x) for x in zs]
z=complex(*z)

Q,aks = create_gravitational_me(qs,zs)
print(run_gravitational_me(Q,aks,z))
print(naive_gravitational(qs,zs,z))


#Lemma 2.3
def shift_outer_me(a0, aks, z0):
  return [sum(aks[k-1]*pow(z0,l-k)*math.comb(l-1,k-1) for k in range(1,l+1)) - a0*pow(z0,l)/l for l in range(1,p+1)]
aks=shift_outer_me(Q,aks,5+5j)
print(run_gravitational_me(Q,aks,z))
zs2=[(6,6),(4,6),(6,4),(4,4)]
zs2=[complex(*x) for x in zs2]
# we have aks == create_gravitational_me(qs,zs2)
print(naive_gravitational(qs,zs2,z))


#Lemma 2.4
def shift_inner_me(a0,aks,z0):
  return [sum(aks[k-1]/pow(z0,k)*pow(-1,k) for k in range(1,p+1)) + a0*cmath.log(-z0)] +\
         [(1/pow(z0,l) * sum(aks[k-1]/pow(z0,k)*math.comb(l+k-1,k-1)*pow(-1,k) for k in range(1,p+1))) - a0/(l*pow(z0,l))
             for l in range(1,p+1)]
def run_potential_inner_me(bs,z):
  return sum(bs[l]*pow(z,l) for l in range(p+1))
def run_gravitational_inner_me(bs,z):
  return sum(l*bs[l]*pow(z,l-1) for l in range(p+1))
_,aks = create_gravitational_me(qs,zs)
aks=shift_inner_me(Q,aks,5+5j)
z=1+0j
print(run_gravitational_inner_me(aks,z))
zs2=[(6,6),(4,6),(6,4),(4,4)]
zs2=[complex(*x) for x in zs2]
print(naive_gravitational(qs,zs2,z))

#step 1: Apply theorem 2.1 to finest boxes about their centers
#step 2: Add and shift for larger boxes (about their centers)
#step 3: Calculate local expansion of everything in the interaction list first (Lemma 2.4).
#        Now add that to the tilde Psi already calculated for parent
#        Calculate tilde Psi for children by shifting using Lemma 2.5
#step 4: Step 3 for final level
#step 5: Evaluate Psi for layer n
#step 6: Manually calculate local points
#step 7: Add step 5 and 6



'''n,N=3,150
class Box(): pass
levels=[[Box(0,(0,0))]]+[[] for _ in range(n)]
for level in range(n):


points=[(random.random()-0.5,random.random()-0.5) for _ in range(N)]
for x,y in points:
  

#Step 1
for box in levels[-1]:
  
'''
