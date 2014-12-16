from numpy import dot, sqrt, array_equal, random, zeros, asarray, ones
from generators import generators_A, generators_D, generators_E, generators_F, generators_H

def orbit(X,n,v=[]):
    if X == 'A': G = generators_A(n)
    elif X == 'D': G = generators_D(n)
    elif X =='E': G = generators_E(n)
    elif X == 'F' : G = generators_F(n)
    elif X == 'H' : G = generators_H(n)
    print(G)
    
    R = all_reflections(G)
    if not v:
        v = random_vector(n,R)
    assert len(v) == n+1,"The dimension of v must be consistent with the dimension n (recall we are in homogeneous coordinates)."
    assert not is_on_hyperplane(v,R), "The vector v passed as argument lies on a hyperplane"
    
    orbit_v = []
    orbit_v.append(v)
    for x in orbit_v:
        for g in G:
            if not array_equal(x,g):
                im_x = reflection(x,g)
                if not contains(orbit_v,im_x): 
                    orbit_v.append(im_x)
    
    # Matrix representation of the hyperplanes and the orbit                    
    matrix_R = zeros((len(R),n+1))
    matrix_orbit = zeros((len(orbit_v),n+1))
    for i in xrange(len(R)): matrix_R[i,:] = R[i]
    for i in xrange(len(orbit_v)): matrix_orbit[i,:] = orbit_v[i]
    print('The %s hyperplanes of reflections of this Coxeter group are (notice that the rotations are not computed by this program) : \n %s' %(len(R),matrix_R))
    print('\nThe orbit of the vector %s is : \n %s' % (v,matrix_orbit))
    return (len(R), len(orbit_v))

def all_reflections(X):
    R = list(X)
    for r in R:
        for x in X:
            if not array_equal(x,r):
                v = reflection(r,x)
                if not contains(R,v): R.append(v)
    return R

def reflection(p,r):
    v = p - 2*r*dot(p,r)/dot(r,r)
    return v

def contains(S,v):
    v = simplification(v)
    for vector in S:
        vector = simplification(vector)
        if array_equal(vector,v) or array_equal(vector,-v): return True
    return False

def simplification(vector):
    if sqrt(dot(vector,vector)) != 0: return vector / sqrt(dot(vector,vector))
    return vector

def random_vector(n,R):
    while True:
        v = random.randint(2*n, size = n+1)
        v[0] = 0.
        if not is_on_hyperplane(v,R): break
    return v
        
def is_on_hyperplane(v, R):
    for r in R:
        if dot(v,r) == 0: 
            return True
    return False    