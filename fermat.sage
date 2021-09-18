# The Hodge conjecture for Fermat varieries x_0^m + .... + x_{n+2}^m of degree m is True if the following claim is true (although it can also be true if the following is false, e.g m=25):
# Let (a_1,...,a_{m-1};y) be the set of solutions of the system of equations:
# \sum_i (i*k mod m)*a_i = m*y for k in (Z/mZ)^*, i.e k invertible
# They form a semigroup under addition. We say (a_1,...,a_{m-1};y) has length y, p=(a_1,...,a_{m-1};y) is indecomposable if is not a sum of two other elements disjoint from p. By Gordon's lemma, there are finitely many indecomposables.
# We say p is quasi-decomposable if p + u = t + q, where u has length 1 and t,q not equal p
# CLAIM: For fixed m, if each indecomposable is quasi-decomposable then the Hodge conjecture is True for every fermat of degree m, regardless of the dimension n.
from itertools import product
from fractions import gcd
import sys
import numpy as np
#from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
def how_many_indec(m):
    p = MixedIntegerLinearProgram(base_ring=QQ)
    w = p.new_variable(integer=True, nonnegative=True)
    for k in range(1,m):
        if gcd(k,m) == 1:
            l=0
            for i in range(1,m):
                l += ((i*k) % m)*w[i-1]
            l += -m*w[m-1]
            #print l
            p.add_constraint(l == 0)
    p.add_constraint(w[m-1] >= 1)
    indec = p.polyhedron(backend='normaliz').integral_points_generators()[0]
    indec_less = [ x for x in indec if x[-1]> 2]
    return len(indec_less)

def poly_sol(m):
    p = MixedIntegerLinearProgram(base_ring=QQ)
    w = p.new_variable(integer=True, nonnegative=True)
    for k in range(1,m):
        if gcd(k,m) == 1:
            l=0
            for i in range(1,m):
                l += ((i*k) % m)*w[i-1]
            l += -m*w[m-1]
            p.add_constraint(l == 0)
    p.add_constraint(w[m-1] >= 1)
    return p.polyhedron(backend='normaliz')

def lengthOne(m):
    p = MixedIntegerLinearProgram(base_ring=QQ)
    w = p.new_variable(integer=True, nonnegative=True)
    for k in range(1,m):
        if gcd(k,m) == 1:
            l=0
            for i in range(1,m):
                l += ((i*k) % m)*w[i-1]
            l += -m*w[m-1]
            p.add_constraint(l == 0)
    p.add_constraint(w[m-1] == 1)
    return p.polyhedron(backend='normaliz').integral_points()

def get_indec(m):
    p = MixedIntegerLinearProgram(base_ring=QQ)
    w = p.new_variable(integer=True, nonnegative=True)
    #print 'x is %d and m is %d' % (x,m)
    for k in range(1,m):
        if gcd(k,m) == 1:
            l=0
            for i in range(1,m):
                l += ((i*k) % m)*w[i-1]
            l += -m*w[m-1]
            #print l
            p.add_constraint(l == 0)
    p.add_constraint(w[m-1] >= 1)
    return p.polyhedron(backend='normaliz').integral_points_generators()[0]


def get_standard(m,primes):
    result = []
    for p in primes:
        d = m/p
        if p == 2:
            for i in range(1,m):
                if (p*i) % m != 0 :#(d/gcd(i,d))>2:#(p*i) % m != 0 and 2*((p*i) % m) != m:##
                    temp = [i,(i+d) % m,(m-2*i) % m,d]
                    #print temp
                    std = []
                    for e in range(1,m):
                        std.append(temp.count(e))
                    std.append(2)
                    if tuple(std) not in result:
                        result.append(tuple(std))
        else:
            for i in range(1,m):
                if (p*i) % m != 0: #and 2*((p*i) % m) != m:#d/gcd(i,d)>2:
                    #print i
                    temp = [0]*(p+1)
                    for k in range(p):
                        temp[k]= (i+k*d) % m
                    temp[p]=(m-p*i) % m
                    #print temp
                    std = []
                    for e in range(1,m):
                        std.append(temp.count(e))
                    std.append((p+1)/2)
                    #print gcd(i,d)
                    #print tuple(std)
                    #print '--'
                    if p%2 == 1:
                        if tuple(std) not in result:
                            result.append(tuple(std))
                    else:
                        if tuple(std) not in result:
                            result.append(tuple(2*x for x in std))
    return result

def reverse_to(y,m):
    r=[]
    n = len(y)-2
    for e in range(1,m):
        r.append(y.count(e))
    r = r + [n/2 + 1]
    return tuple(r)
#for p in product(range(1,33),range(1,33),range(1,33)):
    #print reverse_to((a, -a % 33, b, -b % 33,c, -c % 33),33)

def convert_to_u(x,m):
    last = x[-1]
    n = 2*(last-1)
    r = []
    for k in range(m-1):
        if x[k] != 0:
            r = r + [k+1]*x[k]
    return tuple(r)

def get_linears(m):
    print 'getting evens...'
    evens = get_evens(m)
    result = []
    for ev in evens:
        print "="*10
        print ev
        print '\n'
        count=0
        for p in product(range(1,m),repeat=(ev+2)/2):
            #print p
            l=[]
            for k in p:
                l=l+[k,-k%m]
            result.append(reverse_to(l,m))
            samples=(m-1)**((ev+2)/2)
            count+=1
            sys.stdout.write("Progress: %.2f%%   \r" % (float(100*count)/samples))
            sys.stdout.flush()
        print '\n'
    return list(set(result))

def get_evens(m):
    l=[]
    indy=get_indec(m)
    for n in indy:
        if 2*(n[-1]-1) not in l and n[-1]>2:
            l.append(2*(n[-1]-1))
    print('evens: ',l)
    return l

#for e in get_linears(21):
#    print convert_to_u(e,21)
#exit()

def get_points_length_less_m(x,m):
    p = MixedIntegerLinearProgram(base_ring=QQ)
    w = p.new_variable(integer=True, nonnegative=True)
    #print 'x is %d and m is %d' % (x,m)
    for k in range(1,m):
        if gcd(k,m) == 1:
            l=0
            for i in range(1,m):
                l += ((i*k) % m)*w[i-1]
            l += -m*w[m-1]
            #print l
            p.add_constraint(l == 0)
    p.add_constraint(w[m-1] >= 1)
    p.add_constraint(w[m-1] <= x)
    return p.polyhedron(backend='normaliz').integral_points()

arr = []

def get_indec_less(m,prm):
    p = poly_sol(m)
    print 'getting indecomposable elements for |m= %d| ...' % m
    indec = p.integral_points_generators()[0]
    #length_one = [ x for x in indec if x[-1]==1]
    #print 'there are %d length one' % len(length_one)
    standards = [list(x) for x in get_standard(m,prm)]
    print 'there are %d STANDARDS ELEMENTS' % len(standards)
    indec_less = [ x for x in indec if x[-1]>= 3 and list(x) not in standards]
    print 'there are %d indec of length>=3' % len(indec_less)
    return indec_less
    
def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            if i not in factors:
                factors.append(i)
    if n > 1:
        if n not in factors:
            factors.append(n)
    #print('finished computing primes:',factors)
    return factors

'''with open('somefile.txt', 'a') as the_file:
    for k in range(49,56):
        fac = prime_factors(k)
        if 1:#len(fac)>1: #k % 3 != 0 and k % 2 != 0:# and  and k % 5 != 0 
            #print k
            ind = get_indec_less(k,fac)
            lasts=[]
            if len(ind)>0:
                for x in ind:
                    #print convert_to_u(x,k)
                    if x[-1] not in lasts:
                        lasts.append(x[-1])
            #print '%d : %s indecomposables' % (k,list(lasts))
            #print 'len: %d' % len(ind)
            print k,max(lasts)
            print '\n'
            the_file.write('(%d,%d)' % (k,max(lasts)))
exit()'''
m = 21
primes = prime_factors(m)
quasi = []
dict_ = {}
indec_less = get_indec_less(m,primes)
lasts_ =[]
length_one = lengthOne(m)
for el in indec_less:
    last = el[-1]
    print '---'
    print el
    print 'position: %d' % indec_less.index(el)
    count=0
    if last not in dict_:
        possible = get_points_length_less_m(last,m)
        dict_[last] = possible
        for el2,el3,el4 in product(length_one,possible,possible):
            if el + el2 == el3 + el4 and (el != el3 and el != el4):
                print 'I am quasi'
                quasi.append(el)
                break
            samples=len(possible)*len(possible)*len(length_one)
            count+=1
            sys.stdout.write("Progress: %.2f%%   \r" % (float(100*count)/samples))
            sys.stdout.flush()
        if el not in quasi:
            print 'This element is not quasi'
            print 'The HC CAN NOT be predicted for degree %d using this method, there are only %d quasi of %d' % (m,len(quasi),len(indec_less))
            break
    else:
        for el2,el3,el4 in product(length_one,dict_[last],dict_[last]):
            if el + el2 == el3 + el4 and (el != el3 and el != el4):
                print 'I am quasi'
                quasi.append(el)
                break
            samples=len(dict_[last])*len(dict_[last])*len(length_one)
            count+=1
            sys.stdout.write("Progress: %.2f%%   \r" % (float(100*count)/samples))
            sys.stdout.flush()
        if el not in quasi:
            print 'This element is not quasi'
            print 'The HC CAN NOT be predicted for degree %d using this method, there are only %d quasi of %d' % (m,len(quasi),len(indec_less))
            break
print 'The HC is TRUE for degree %d fermats' % m