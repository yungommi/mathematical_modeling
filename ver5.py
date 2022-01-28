# open data
lines=open('ev.txt','r')
e=list(lines)
del(e[0])

a=[]
c=[]
E=[]
for line in e:
    l=line.split()
    a.append(float(l[0]))
    c.append(float(l[1]))
    E.append(float(l[2]))


import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

# 1. make a contour plot E(a,c)
npoint=200
e_a=np.linspace(min(a),max(a), npoint)
e_c=np.linspace(min(c),max(c), npoint)
e_a, e_c = np.meshgrid(e_a,e_c)
E_cubic=scipy.interpolate.griddata((a,c), E, (e_a, e_c), method='cubic') #cubic spline

plt.subplot(221)
plt.title('Fig. 1')
plt.pcolormesh(e_a,e_c, E_cubic, cmap=plt.get_cmap('rainbow'))
plt.xlabel('a ($\AA$)')
plt.ylabel('c ($\AA$)')
plt.colorbar()
plt.contour(e_a, e_c, E_cubic, 10, linewidths=0.5, color='k')


# 2.
import math
a1=np.array(a)
c1=np.array(c)
V=(3*math.sqrt(3)/2)*a1**2*c1 # caculate V
E2=[]
VE=sorted([[V[i],E[i]] for i in range(len(E))]) # sorting E by V
n=3
    
while n<len(V):
        E1=[]
        for i in range(n-3,n):
            E1.append([VE[i][1],VE[i][0]]) # selct points having similar volume
            E1=sorted(E1) # select lowest energy 
        
        E2.append(E1[0]) # E2= selected data
        n=n+3
#selected data
new_V=[E2[i][1] for i in range(len(E2))]
new_E=[E2[i][0] for i in range(len(E2))]


plt.subplot(222)
plt.plot(V,E,'d',color='r',marker='d',label='all given data')
plt.plot(new_V,new_E,'d',color='b', marker='o', label='selected data')
plt.title('Fig. 2')
plt.xlabel('Volume ($\AA\ ^3$)')
plt.ylabel('Energy (eV)')
plt.legend()



#3. polynomial fitting
import scipy.optimize
def rand_func(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d # cubic polynomial fitting

#plot Fig.3
new_V=np.array(new_V) 
new_E=np.array(new_E)
plt.subplot(223)
plt.plot(new_V,new_E,'d',color='r', marker='d',label='selected points') #selected data
popt, pcov= scipy.optimize.curve_fit(rand_func,new_V,new_E)
perr=np.sqrt(np.diag(pcov))
plt.title('Fig. 3')
plt.xlabel('Volume ($\AA\ ^3$)')
plt.ylabel('Energy (eV)')
plt.plot(new_V, rand_func(new_V, *popt),color = 'b', label='fit: a=%5.3f$\pm$%3.2f, b=%5.3f$\pm$%3.2f, c=%5.3f$\pm$%3.2f, d=%5.3f$\pm$%3.2f' %(popt[0],perr[0],popt[1],perr[1],popt[2],perr[2],popt[3],perr[3]))
plt.legend()

a=popt[0]; b=popt[1]; c=popt[2]; d=popt[3]
print '3.\nThe fitted polynomial is y= %5.5f x^3 + %5.5f x^2 + %5.5f x + %5.5f' %(a,b,c,d) # print fitted polynomial

# calculation of error
error= new_E- a*new_V**3 - b*new_V**2 - c*new_V -d
Sr=sum(error**2)
meanE=sum(new_E)/len(new_E)
St=sum((new_E-meanE)**2)
r2=(St-Sr)/St
print 'error of fitted polynomial: Sr=%5.3f, St=%5.3f, r^2=%5.3f'%(Sr,St,r2)

# Bisection
print '\n< Bisection method >'
def r_f(x):
    return 3*a*x**2+2*b*x+c # first derivative of fitted parabola

tol = 0.01
x_low = 110
x_high = 150

print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %('Iteration', 'x_low', 'x_high', 'x_r', 'f(x_r)', 'error(%)')
print "____________________________________________________________________"

i = 1
err = 100

while abs(err) >= tol:
    x_r = (x_low + x_high)/2.
    f_xr = r_f(x_r)
 
    o_i = i #iteration
    o_x_low = x_low
    o_x_high = x_high
    o_x_r = x_r
    o_f_xr = f_xr

    if i==1:
        pass # not to be able to caculate approximate error in 1st iteration
    else:
        err=(previous_x_r-x_r)/previous_x_r*100 # Ea
       
    
    if r_f(x_low)*f_xr < 0:
        x_high = x_r
    else:
        x_low = x_r
    i = i+1
    previous_x_r=x_r
    
  
    print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|"\
          %(o_i, o_x_low, o_x_high, o_x_r, o_f_xr, err)
V0=o_x_r; E0= rand_func(o_x_r,popt[0],popt[1],popt[2],popt[3])
print '\nV0 is %10.8s (Angstrom^3) and E0 is %10.8s (eV)'%(V0,E0)



# open method - Newton-Raphson method
print '\n< Newton-Raphson method >'

print "|%9.9s|%10.8s|%10.8s|%12.10s|"\
      %('Iteration', 'x_i', 'x_i+1', 'error(%)')
print "______________________________________________"

def first_d(x):
    return 6*a*x+2*b # first derivative of the function(first derivative of fitted cubic polynomial)

x_i=155
i = 1
err = 100

while abs(err) >= tol:
    
    previous_x=x_i
    x_i1=x_i-r_f(x_i)/first_d(x_i) 
    o_i = i
    o_xi=x_i
    
    if i==1:
        pass # not to be able to caculate approximate error in 1st iteration
    else:
        err=(previous_x-x_i1)/previous_x*100 # Ea
    i+=1   
    x_i=x_i1
  
    print "|%9.9s|%10.8s|%10.8s|%12.7s|"\
          %(o_i, o_xi,x_i1, err)

print '\nV0 is %10.8s (Angstrom^3) and E0 is %10.8s (eV)'%(x_i1,rand_func(x_i1,popt[0],popt[1],popt[2],popt[3]))

#4
from scipy.optimize import leastsq
import numpy as np

def Murnaghan(V, parameters): #3rd order EOS
    E0=parameters[0]; B0=parameters[1]; b0=parameters[2]; V0 = parameters[3]
    V2=(V0/V)**(2.0/3.0)
    return E0+9.0/16.0*B0*V0*((V2-1.0)**3.0*b0-(V2-1.0)**2.0*(4.0*V2-6.0))
    
def residuals(pars, y, x):
    #we will minimize this function's value
    return y - Murnaghan(x, pars)

x0 = [ E0, 0.54, 2, V0] #initial guess of parameters [E0, B0, B'0, V0]

plsq = leastsq(residuals, x0, args=(new_E, new_V))
print '\n4.\nThe coefficients of fitted function are\nE0=%5.5f(eV) B0=%5.5f(GPa) B\'0=%5.5f, V0=%5.5f(Angstrom^3)'%(plsq[0][0],plsq[0][1],plsq[0][2],plsq[0][3])

#plot the fitted curve on top
x = np.linspace(new_V.min(), new_V.max(), 50)
y = Murnaghan(x, plsq[0])
plt.subplot(224)
plt.plot(new_V,new_E, 'ro', label='selected points') #plot selected points
plt.plot(x,y, 'k-',label='third-order Birch-Murnaghan EOS fitting')
plt.title('Fig. 4')
plt.xlabel('Volume ($\AA\ ^3$)')
plt.ylabel('Energy (eV)')
plt.legend()

plt.show()

