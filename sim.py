import cmath, math, random

p = 20
N = 1000
n=int(math.log(N)/math.log(10))

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


'''
zs=[(1,1),(-1,1),(1,-1),(-1,-1)]
qs=[1,1,1,1]
z=(10,10)

zs=[complex(*x) for x in zs]
z=complex(*z)

Q,aks = create_gravitational_me(qs,zs)
print(run_gravitational_me(Q,aks,z))
print(naive_gravitational(qs,zs,z))
'''


#Lemma 2.3
def shift_outer_me(a0, aks, z0):
  return [sum(aks[k-1]*pow(z0,l-k)*math.comb(l-1,k-1) for k in range(1,l+1)) - a0*pow(z0,l)/l for l in range(1,p+1)]
'''
aks=shift_outer_me(Q,aks,5+5j)
print(run_gravitational_me(Q,aks,z))
zs2=[(6,6),(4,6),(6,4),(4,4)]
zs2=[complex(*x) for x in zs2]
# we have aks == create_gravitational_me(qs,zs2)
print(naive_gravitational(qs,zs2,z))
'''


#Lemma 2.4
def shift_inner_me(a0,aks,z0):
  return [sum(aks[k-1]/pow(z0,k)*pow(-1,k) for k in range(1,p+1)) + a0*cmath.log(-z0)] +\
         [(1/pow(z0,l) * sum(aks[k-1]/pow(z0,k)*math.comb(l+k-1,k-1)*pow(-1,k) for k in range(1,p+1))) - a0/(l*pow(z0,l))
             for l in range(1,p+1)]
def run_potential_inner_me(bs,z):
  return sum(bs[l]*pow(z,l) for l in range(p+1))
def run_gravitational_inner_me(bs,z):
  return sum(l*bs[l]*pow(z,l-1) for l in range(p+1))
'''
_,aks = create_gravitational_me(qs,zs)
aks=shift_inner_me(Q,aks,5+5j)
z=1+0j
print(run_gravitational_inner_me(aks,z))
zs2=[(6,6),(4,6),(6,4),(4,4)]
zs2=[complex(*x) for x in zs2]
print(naive_gravitational(qs,zs2,z))
'''

#Lemma 2.5
def shift_power_series_center(aks,z0):
  return [sum(aks[k]*math.comb(k,l)*pow(-z0,k-l) for k in range(l,p+1)) for l in range(p+1)]
'''
aks=shift_power_series_center(aks,-1-1j)
print(run_gravitational_inner_me(aks,-2-2j))
_,aks = create_gravitational_me(qs,zs2)
print(naive_gravitational(qs,zs2,-1-1j))
'''

class Box():
  def __init__(self,level,i,j):
    self.level,self.i,self.j,self.particles=level,i,j,[]
    self.center=complex((i+0.5)*pow(2,-level)-0.5,(j+0.5)*pow(2,-level)-0.5)
    self.Q=0
    self.Phi=[0 for _ in range(p+1)]
    self.Psi=[0 for _ in range(p+1)]
    self.Psi2=[0 for _ in range(p+1)]
  def __str__(self):
    return str(self.i)+','+str(self.j)

boxes=[]
for level in range(n+1):
  boxes.append([])
  for i in range(2**level):
    boxes[-1].append([])
    for j in range(2**level):
      boxes[-1][-1].append(Box(level,i,j))

particles=[(random.random()-0.5,random.random()-0.5) for _ in range(N)]
for x,y in particles:
  i,j=int((x+0.5)*2**n), int((y+0.5)*2**n)
  boxes[-1][i][j].particles.append(complex(x,y))

def get_interaction_list(level,i,j):
  r=[]
  for di1 in range(-1,2):
    for dj1 in range(-1,2):
      for di2 in range(2):
        for dj2 in range(2):
          i2,j2=(i//2+di1)*2+di2 ,(j//2+dj1)*2+dj2
          if i2 in range(2**level) and j2 in range(2**level) and (abs(i2-i)>1 or abs(j2-j)>1):
            r.append((i2,j2))
  return r

#Step 1
print('Step 1')
for row in boxes[-1]:
  for box in row:
    box.Q, box.Phi = create_gravitational_me([1 for _ in range(len(box.particles))],
                                             [x-box.center for x in box.particles])

'''
tmp=boxes[-1][-1][-1]
c=-1-1j
#print(tmp.center,c-tmp.center)
print(run_gravitational_me(tmp.Q,tmp.Phi,c-tmp.center))
#print([x+tmp.center for x in tmp.particles])
print(naive_gravitational([1]*len(tmp.particles),tmp.particles,c))
'''

#Step 2
print('Step 2')
for level in range(n-1,-1,-1):
  for row in boxes[level]:
    for box in row:
      for di in range(2):
        for dj in range(2):
          child = boxes[level+1][box.i*2+di][box.j*2+dj]
          shifted=shift_outer_me(child.Q, child.Phi, child.center-box.center)
          box.Phi=[x+y for x,y in zip(box.Phi,shifted)]
          box.Q+=child.Q

'''
tmp=boxes[-2][-1][-1]
c=-1-1j
print(run_gravitational_me(tmp.Q,tmp.Phi,c-tmp.center))
tmp2=boxes[-1][-1][-1].particles+boxes[-1][-1][-2].particles+boxes[-1][-2][-1].particles+boxes[-1][-2][-2].particles
print(naive_gravitational([1]*len(tmp2),tmp2,c))
'''



#Step 3
print('Step 3')
for level in range(1,n):
  for row in boxes[level]:
    for box in row:
      for i,j in get_interaction_list(level,box.i,box.j):
        box2=boxes[level][i][j]
        shifted = shift_inner_me(box2.Q,box2.Phi,box2.center-box.center)
        box.Psi=[x+y for x,y in zip(box.Psi,shifted)]
      box.Psi=[x+y for x,y in zip(box.Psi,box.Psi2)]
  for row in boxes[level]:
    for box in row:
      for di in range(2):
        for dj in range(2):
          child = boxes[level+1][box.i*2+di][box.j*2+dj]
          shifted = shift_power_series_center(box.Psi, child.center-box.center) #not sure
          child.Psi2=shifted

#Step 4
print('Step 4')
for row in boxes[n]:
  for box in row:
    for i,j in get_interaction_list(n,box.i,box.j):
      box2=boxes[n][i][j]
      shifted = shift_inner_me(box2.Q,box2.Phi,box2.center-box.center)
      box.Psi=[x+y for x,y in zip(box.Psi,shifted)]
    box.Psi=[x+y for x,y in zip(box.Psi,box.Psi2)]

#Step 5
print('Step 5')
for row in boxes[n]:
  for box in row:
    neighbor_particles=set()
    for i in range(box.i-1,box.i+2):
      for j in range(box.j-1,box.j+2):
        if i in range(2**n) and j in range(2**n):
          neighbor_particles.update(boxes[n][i][j].particles)
    box.forces=[]
    for x in box.particles:
      neighbor_particles.remove(x)
      box.forces.append(naive_gravitational([1]*len(neighbor_particles),neighbor_particles,x))
      neighbor_particles.add(x)

#Step 6 + 7
print('Step 6 + 7')
particles_final=[]
forces_final=[]
for row in boxes[n]:
  for box in row:
    for i,x in enumerate(box.particles):
      box.forces[i]+=run_gravitational_inner_me(box.Psi,x-box.center)
      particles_final.append(x)
      forces_final.append(box.forces[i])


print(len(particles_final),len(forces_final))
forces_manual=[]
particles_set=set(particles_final)
for particle in particles_final:
  particles_set.remove(particle)
  forces_manual.append(naive_gravitational([1]*len(particles_set),particles_set,particle))
  particles_set.add(particle)

for _ in range(5):
  i=random.randint(0,N-1)
  print(forces_final[i],'=?',forces_manual[i])




