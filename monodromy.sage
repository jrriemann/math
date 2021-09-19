#For a family of elliptic curves over [0,1,infinity], this code graphs the result of the monodromy action on the base space.
import numpy as np
var('z x')
theta = np.linspace(0*np.pi, 2*np.pi, 40)
#pt =  2/(3*np.sqrt(3)) #(-1+I*np.sqrt(3))/2
for pt in [0,1]:#[2/(3*np.sqrt(3)),-2/(3*np.sqrt(3)),0]:
    zeds =  pt+0.1*(np.cos(theta) + I*np.sin(theta)) 
    data1=[]
    data2 = []
    data3 = []
    sol = solve(z*(x^3-1)+z^2, x) #solve(z*(x-1)*(x+1)*x+z^2, x) #
    for v in zeds.tolist():
        #print v
        data1.append(sol[0].rhs().substitute(z=CDF(v)))
        data2.append(sol[1].rhs().substitute(z=CDF(v)))
        data3.append(sol[2].rhs().substitute(z=CDF(v)))
    points([CDF(p) for p in data1],color='red') + points([CDF(p) for p in data2],color='green') + points([CDF(p) for p in data3],color='blue')
