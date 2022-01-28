import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.optimize import leastsq

#Reading the given data (ev.txt)
a = []; c = []; V = []; E = []; i = 1
for line in open('ev.txt'):
    if i ==1 :
        pass
    else :
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
        V.append(float(float(a[i-2])**2 * (np.sqrt(3))/2 * float(c[i-2])))
    i = i + 1

#Problem 1
a_plot = np.linspace(a[0],a[-1],63)
c_plot = np.linspace(c[0],c[-1],63)
a_plot,c_plot = np.meshgrid(a_plot,c_plot)

#3d plot
cubic = scipy.interpolate.griddata((a,c),E,(a_plot,c_plot),method='cubic')
plt.subplot(221)
plt.title('Fig.1')
plt.pcolormesh(a_plot,c_plot,cubic,cmap=plt.get_cmap('rainbow'))
plt.colorbar()
plt.contour(a_plot,c_plot,cubic,63,linewidths=0.5,colors='k')


#Problem 2(find the lowest energy)
V_E = []
for i in range(0,len(V)):
    V_E.append([V[i],E[i]])
V_E = sorted(V_E) #Rearrange the Volume data in ascending order.

lowestE = [] # Make the data in this sequence;[Volume, Energy] 
i = 0
while i <= 60 :   # If, E_n-1 < E_n < E_n+1 choose E_n and save its volume data in lowestE
    if V_E[i][1] > V_E[i+1][1] and V_E[i+1][1] < V_E[i+2][1]:
        lowestE.append(V_E[i+1])
    else:
        pass
    i = i + 1
    
x_V = []; y_E = [] #Split the volume and energy data in lowestE=[]
for i in range(0,len(lowestE)):
    x_V.append(lowestE[i][0]) #Volume   
    y_E.append(lowestE[i][1]) #Energy

plt.subplot(222) #The graph of the problem 2
plt.title('FIg.2')
plt.plot(V,E,'o',label = 'all given data')
plt.plot(x_V,y_E,'d',color = 'r', label = 'selected data')
plt.xlabel('Volume(A^3)')
plt.ylabel('Total Energy(eV)')
plt.legend()


# Problem 3
 #curve-fitting in terms of cubic function (To satisfy the condtion of the problem 3 ['iterations should be printed more than 5 lines until the error is lower than 1%])
coef,residual,_,_,_ = np.polyfit(x_V,y_E,3,full = True) 
x_plot = np.linspace(x_V[0],x_V[-1]) # Make the x range (Refer to the 'lowest E')
y_plot = np.polyval(coef,x_plot)
error = float(residual)*100
print 'fit: E = %f^3 + %fV^2 + %fV + %f' %(coef[0],coef[1],coef[2],coef[3])
print 'the error of fitting is %f%s' %(error,'%')

# Make the function to find the root of the equation
def E0(V):  #The graph of the problem 2
    return (coef[0]*V**3 + coef[1]*V**2 + coef[2]*V + coef[3])
def dE(V):  #dE/dV
    return (coef[0]*3*V**2 + coef[1]*2*V + coef[2])
def ddE(V):  #d^2E/dV^2
    return (coef[0]*3*2*V + coef[1]*2)

# Using bracket method, bisection method, To find the point of centraflexure, find the root of the 1st derivative of equation 
V_low = x_V[0]; V_high = x_V[-1] #Make the starting point (Use the first data and last data from the x_V)
tol = 0.08 # Make the error less than 1%
#Bisection method
print''
print '--Using bisection method(braketing the root)--'
print"|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|" %('iteration','V_low','V_high','V0','E0','error(%)')
print"----------------------------------------------------------------"
i = 1; err_b = 100
while abs(err_b) > tol:
    V_r = (V_low + V_high)/2
    if i == 1:
        print"|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|" %(i,V_low,V_high,V_r,E0(V_r),'-')
    else:
        if dE(V_low)*dE(V_r)<0: # By intermediate value theorem, reset the V-high or V_low (BIsection method)
            V_high = V_r
        else:
            V_low = V_r
        newV_r = (V_low + V_high)/2
        err_b = ((newV_r - V_r)/newV_r)*100
        print"|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|" %(i,V_low,V_high,newV_r,E0(V_r),abs(err_b))      
    i = i+1
print 'V0 is %f, E0 is %f' %(newV_r,E0(newV_r))

#Using open method, Newton-raphson method
print''
print'--Using Newton-Raphson method(open method)--'
print"|%9.9s|%10.8s|%10.8s|%10.8s|%12.10s|" %('iteration','V_initial','V0','E0','error(%)')
print"----------------------------------------------------------------"

V = x_V[-1]; i = 1; err_0 = 100
while abs(err_0) > tol:
    newV = V - dE(V)/ddE(V) # The equation of the Newton-raphson method
    err_0 = ((newV - V)/newV)*100
    print"|%9.9s|%10.8s|%10.8s|%10.8s|%12.10s|" %(i,V,newV,E0(newV),abs(err_0))
    V = newV; i = i + 1        
print 'V0 = %f, E0 = %f' %(V,E0(V))
# We can find out the Newton-raphson method's iteration is less than bisection method. This refer to in some case except for the unique equations, the Newton-raphson is much more useful to find the root we want.
#The graph of the problem 3
plt.subplot(223)
plt.title('Fig.3')
plt.plot(x_plot,y_plot,label = 'fitted curve')
plt.plot(x_V,y_E,'d',color = 'r', label = 'Selected data')
plt.xlabel('Volume(A^3)')
plt.xlim(x_V[0],x_V[-1])
plt.ylabel('Energy(eV)')
plt.plot(newV_r,E0(newV_r),'o',color = 'y',label = 'ground state(bisection method)')
plt.plot(V,E0(V),'<',color = 'g', label = 'Newton method')
plt.legend(loc = 'best')

#Problem 4
co= np.polyfit(x_V,y_E,2)

def En0(V):  #The curve from the Problem 2 ( Make the parabolic equation to fit third-order Borch-Murnaghan equation)
    return (co[0]*V**2 + co[1]*V + co[2])
def dEn(V):  #dE/dV
    return (co[0]*2*V + co[1])
def ddEn(V):  #d^2E/dV^2
    return (co[0]*2)
#Using Open method ( Newton-Raphson mehod), find the ground state volume
V = x_V[-1]; i = 1; err_0 = 100
while abs(err_0) > tol:
    newV = V - dEn(V)/ddEn(V)
    err_0 = ((newV - V)/newV)*100
    V = newV; i = i + 1      

V0 = V; E0 = En0(V0) #The ground state data by the open method ( Newton-Raphson method)
 # The coefficent of the qudratic term of the curve from the Problem 2
B0 = V0 * 2 * co[0]; BP = 2 * co[0]

def rand_func(V1,a,b,c):  # Birch-Murnaghan equation ( Problem 4)
    Ev = E0+(9.0/16.0)*B0*V0*((((V0/V1)**(2.0/3.0)-1.0)**3.0)*BP-(((V0/V1)**(2.0/3.0)-1.0)**2.0)*(4.0*(V0/V1)**(2.0/3.0)-6.0))
    return Ev

vdata = np.array(x_V)
edata = np.array(y_E)

z = rand_func(vdata,co[0],co[1],co[2])
popt,pcov = scipy.optimize.curve_fit(rand_func,vdata,edata)

#Print the ground state and bulk modulus
print 'The ground state total energy is %f [eV], the ground state volume is %f [A^3],\
the bulk modulus is %f [GPa], and the derivative of bulk modulus is %f.'\
%(E0,V0,B0,BP)


#The graph of the Problem 4
plt.subplot(224)
plt.title('Fig.4')
plt.plot(vdata,rand_func(vdata,*popt),'-r',label='Birch-Murnaghan')
plt.plot(vdata,rand_func(vdata,*popt),'d')
plt.xlabel('Volume(A^3)')
plt.ylabel('Energy(eV)')
plt.legend()
plt.show()
