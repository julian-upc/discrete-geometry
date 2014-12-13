from numpy import zeros, sqrt, array

def generators_A(n):
    A = []
    for i in xrange(n):
        v = zeros(n+1)
        v[i+1] = 1.0
        if i != n-1 : v[i+2] = -1.0
        A.append(v)
    return A

def generators_D(n):   
    D = []
    for i in xrange(n-1):
        v = zeros(n+1)
        v[i+1] = 1.0
        v[i+2] = -1.0
        D.append(v)
    v = zeros(n+1)
    v[n-1] = 1.0
    v[n] = 1.0
    D.append(v)
    return D

def generators_E(n):
    assert n==8, 'The group E for n %s is not defined ! Try again with n=8.' % n  
    E = []
    for i in xrange(n-2):
        v = zeros(n+1)
        v[i+1] = 1.0
        v[i+2] = -1.0
        E.append(v)
    E.append(array((0.,0.,0.,0.,0.,0.,1.,1.,0.)))
    E.append(array((0.,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)))      
    return E

def generators_F(n):
    assert n==4, 'The group F for n %s is not defined ! Try again with n=4.' % n  
    return [array((0.,1.,-1.,0.,0.)), array((0.,0.,1.,-1.,0.)), array((0.,0.,0.,1.,0.)), array((0.,-0.5,-0.5,-0.5,-0.5))] 

def generators_H(n):
    assert n==4 or n==3, 'The group H for n %s is not defined ! Try again with n=3 or n=4.' % n
    t = 0.5+0.5*sqrt(5)
    if n==3: return [array((0.,2.,0.,0.)), array((0.,-t,1.0/t,-1.)), array((0.,0.,0.,2.))]
    return [array((0.,-t,1.0/t,1.0/t,1.0/t)), array((0.,-1.,1.,0.,0.)), array((0.,0.,-1.,1.,0.)), array((0.,0.,-1.,1.,0.))] 
