import cmath, math, random, time, matplotlib.pyplot as plt

def naive_gravitational(qs,zs,z):
  return sum(qi/(z-zi) for qi,zi in zip(qs,zs))

#Theorem 2.1
def create_gravitational_me(qs,zs):
  return (sum(qs), [sum(-qi*pow(zi,k)/k for (qi,zi) in zip(qs,zs)) for k in range(1,p+1)])
def run_gravitational_me(Q,aks,z):
  return Q/z+sum(-k*aks[k-1]/pow(z,k+1) for k in range(1,p+1))

#Lemma 2.3
def shift_outer_me(a0, aks, z0):
  return [sum(aks[k-1]*pow(z0,l-k)*math.comb(l-1,k-1) for k in range(1,l+1)) - a0*pow(z0,l)/l for l in range(1,p+1)]

#Lemma 2.4
def shift_inner_me(a0,aks,z0):
  return [sum(aks[k-1]/pow(z0,k)*pow(-1,k) for k in range(1,p+1)) + a0*cmath.log(-z0)] +\
         [(1/pow(z0,l) * sum(aks[k-1]/pow(z0,k)*math.comb(l+k-1,k-1)*pow(-1,k) for k in range(1,p+1))) - a0/(l*pow(z0,l))
             for l in range(1,p+1)]
def run_gravitational_inner_me(bs,z):
  return sum(l*bs[l]*pow(z,l-1) for l in range(p+1))

#Lemma 2.5
def shift_power_series_center(aks,z0):
  return [sum(aks[k]*math.comb(k,l)*pow(-z0,k-l) for k in range(l,p+1)) for l in range(p+1)]

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

def FMM(n,N,p):
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
  for row in boxes[-1]:
    for box in row:
      box.Q, box.Phi = create_gravitational_me([1 for _ in range(len(box.particles))],
                                               [x-box.center for x in box.particles])

  #Step 2
  for level in range(n-1,-1,-1):
    for row in boxes[level]:
      for box in row:
        for di in range(2):
          for dj in range(2):
            child = boxes[level+1][box.i*2+di][box.j*2+dj]
            shifted=shift_outer_me(child.Q, child.Phi, child.center-box.center)
            box.Phi=[x+y for x,y in zip(box.Phi,shifted)]
            box.Q+=child.Q

  #Step 3
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
            shifted = shift_power_series_center(box.Psi, box.center-child.center)
            child.Psi2=shifted

  #Step 4
  for row in boxes[n]:
    for box in row:
      for i,j in get_interaction_list(n,box.i,box.j):
        box2=boxes[n][i][j]
        shifted = shift_inner_me(box2.Q,box2.Phi,box2.center-box.center)
        box.Psi=[x+y for x,y in zip(box.Psi,shifted)]
      box.Psi=[x+y for x,y in zip(box.Psi,box.Psi2)]

  #Step 5
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
  particles_final=[]
  forces_final=[]
  for row in boxes[n]:
    for box in row:
      for i,x in enumerate(box.particles):
        box.forces[i]+=run_gravitational_inner_me(box.Psi,x-box.center)
        particles_final.append(x)
        forces_final.append(box.forces[i])
  return particles_final,forces_final


def naive_forces(particles_final):
  forces_manual=[]
  particles_set=set(particles_final)
  for particle in particles_final:
    particles_set.remove(particle)
    forces_manual.append(naive_gravitational([1]*len(particles_set),particles_set,particle))
    particles_set.add(particle)
  return forces_manual


p = 30
n = 4
Ns,FMM_timings,naive_timings=[],[],[]
for N in range(1000,15001,1000):
  Ns.append(N)
  t1=time.time()
  particles,forces_FMM = FMM(n,N,p)
  t2=time.time()
  forces_manual = naive_forces(particles)
  t3=time.time()
  FMM_timings.append(t2-t1)
  naive_timings.append(t3-t2)
  print('N:',N)
  print('FMM Time:',t2-t1)
  print('Naive Time:',t3-t2)
  for _ in range(5):
    i=random.randint(0,N-1)
    print('  err:',abs(forces_FMM[i]-forces_manual[i]))
plt.plot(Ns,naive_timings,label="Naive")
plt.plot(Ns,FMM_timings,label="FMM")
plt.xlabel("N")
plt.ylabel("Runtime (s)")
plt.legend()
plt.show()



