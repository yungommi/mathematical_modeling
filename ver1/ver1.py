#Problem 1

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

#Read Data from ev.txt and store them into a, c, E.
a = []; c = []; E = []; i = 1
for line in open("ev.txt"):
    if i == 1:
        pass
    else:
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i = i + 1
    
#Count the points of a and c by using len().     
npoint_a = len(a)
npoint_c = len(c)

#Setting up a regular grid of spline points.
g_a = np.linspace(min(a), max(a), npoint_a)
g_c = np.linspace(min(c), max(c), npoint_c)
g_a, g_c = np.meshgrid(g_a, g_c)

#Interpolate with cubic spline method by using scipy.interpolate.
E_cubic = scipy.interpolate.griddata((a, c), E, (g_a, g_c), method='cubic')

#Plot the 3D contour with cubic spline.
plt.subplot(221)
plt.pcolormesh(g_a, g_c, E_cubic, cmap = plt.get_cmap('rainbow'))
plt.colorbar()
plt.contour(g_a, g_c, E_cubic,linewidths = 0.5, color='k')
plt.title('Fig. 1')

#Label x and y axis.
plt.xlabel('a (Angstrom)')
plt.ylabel('c (Angstrom)')
plt.legend()

#Problem 2

import math

#Make the a,c, E lists to array and calculate the volume of primitive unit cell of hcp. 
a_v=np.array(a)
c_v=np.array(c)
E_v=np.array(E)
v=(math.sqrt(3)*a_v**2*c_v)/2


#Store the volume data and the energy at each volume by list in X.
X=[]
for j in range(len(v)):
   X.append([v[j],E[j]])

#Sort the X from the lowest volume to the highest volume.
v_sort=sorted(X)

v_low=[]; E_low=[]
#Group four data points of the given data set which have similar volume.
for i in range(len(v)/4):
#Reverse the list of data points by using [::-1] and choose the data point which have the lowest energy comparing other three points having similar volume.    
    lowest_E_of_similar_v=min((v_sort[4*i][::-1],v_sort[4*i+1][::-1],v_sort[4*i+2][::-1],v_sort[4*i+3][::-1]))
#Then, seperate the volume and energy at the point of the lowest energy and store them into each list in order.
    v_low.append(lowest_E_of_similar_v[1])
    E_low.append(lowest_E_of_similar_v[0])

#Plot the given data by color of red and the selected data by color of black.
plt.subplot(222)
plt.plot(v, E_v, 'd', color='r', label='all given data')
plt.plot(v_low, E_low, 'd', color='k', label='selected data')
plt.xlabel('Volume (Angstrom^3)')
plt.ylabel('Energy (eV)')
plt.title('Fig. 2')
plt.legend()

#Problem 3

print ("\nProblem 3")

#We chose a third-order polynomial to fit the selected points in problem 2.
#y = a0 + a1*x + a2+x^2 + a3*x^3 + e
#By the least squares, the coefficients yielding the minimum sum of the squares of the residuals are obtained by setting the partial derivatives equal to zero.
#As we chose a third-order polynomial, we need (4,4) matrix and (4,1) matrix.
#The summations required to develop the (4,4) matrix and (4,1) matrix are computed like below.

n=len(v_low)
sum_of_v=0; sum_of_v_2=0; sum_of_v_3=0; sum_of_v_4=0; sum_of_v_5=0; sum_of_v_6=0

for i in range(n):
    sum_of_v = sum_of_v + v_low[i]
    #Sum of v
    sum_of_v_2 = sum_of_v_2 + v_low[i]**2
    #Sum of v^2
    sum_of_v_3 = sum_of_v_3 + v_low[i]**3
    #Sum of v^3
    sum_of_v_4 = sum_of_v_4 + v_low[i]**4
    #Sum of v^4
    sum_of_v_5 = sum_of_v_5 + v_low[i]**5
    #Sum of v^5
    sum_of_v_6 = sum_of_v_6 + v_low[i]**6
    #Sum of v^6
sum_v=[n, sum_of_v, sum_of_v_2, sum_of_v_3, sum_of_v_4, sum_of_v_5, sum_of_v_6]

#Make v_M which is a (4,4) matrix.  
v_M = np.zeros((4,4))
i = 0
while i < 4 :
    v_M[i,0] = sum_v[i]
    v_M[i,1] = sum_v[i+1]
    v_M[i,2] = sum_v[i+2]
    v_M[i,3] = sum_v[i+3]                                        
    i = i+1


sum_of_E=0; sum_of_vE=0; sum_of_v_2E=0; sum_of_v_3E=0;

for i in range(n):
    sum_of_E = sum_of_E + E_low[i]
    #Sum of E
    sum_of_vE = sum_of_vE + E_low[i] * v_low[i]
    #Sum of v*E
    sum_of_v_2E = sum_of_v_2E + E_low[i] * (v_low[i]**2)
    #Sum of v^2*E
    sum_of_v_3E = sum_of_v_3E + E_low[i] * (v_low[i]**3)
    #Sum of v^3*E

#Make E_M which is a (4,1) matrix.
E_M = np.zeros((4,1))
E_M[0,0] = sum_of_E
E_M[1,0] = sum_of_vE
E_M[2,0] = sum_of_v_2E
E_M[3,0] = sum_of_v_3E

#Inverse the v_M matrix and multiply the v_inv and E_M matrix to get coefficients of the polynomial.   
v_inv = np.linalg.inv(v_M)
coeff = np.matmul(v_inv, E_M)

#Print the equation of the fitted third-order polynomial.
print ("\n")
print ("The coeffiencts of the third order polynomial: [a0 = %.8f , a1 = %.8f , a2 = %.8f , a3 = %.8f ]"\
      %(coeff[0],coeff[1],coeff[2],coeff[3])
print("\nThe polynomial fitting result is: y= %.8f + %.8f x + %.8f x^2 + %8f x^3"
      %(coeff[0],coeff[1],coeff[2],coeff[3])

    
#Plot the curve of the fitted third-order polynomial.
x=np.linspace(30,60,60)
y=coeff[0] + coeff[1]*x + coeff[2]*x**2 + coeff[3]*x**3
plt.subplot(223)
plt.plot(x, y, '--', color='r', label='polynomial fitting')
plt.plot(v_low, E_low, 'd', color='k', label='selected data')
plt.xlabel('Volume (Angstrom^3)')
plt.ylabel('Energy (eV)')
plt.title('Fig. 3')
plt.legend()

#Define the function of the fitted third-order polynomial.
def f(x):
    return float(coeff[0] + coeff[1]*x + coeff[2]*(x**2) + coeff[3]*(x**3))

#Print the percentage error of each selected data.
print "\n"
print "|%10.10s|%25.25s|%20.20s|%20.20s|%20.20s|"\
      %('Data No.','percentage error (%)','Fitted energy (eV)','Original energy (eV)','Volume (Angstrom^3)')
print "------------------------------------------------------------------------------------------------------------------------"
error_of_data=0
for i in range(n):
    error_of_data=abs(f(v_low[i])-E_low[i])*100
    print "|%10.8s|%25.8s|%20.8s|%20.8s|%20.8s| "\
          %(i+1,error_of_data, f(v_low[i]), E_low[i], v_low[i])


#At first, find the ground state volume and energy by golden-section search.
print "\n"
print "By golden-section search"
print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.9s|%10.8s|%10.8s|%12.10s|"\
      %('iteration','x_low','f(x_low)','x_2','f(x_2)','x_1','f(x_1)','x_high','f(x_high)','x_opt','f(x_opt)','error(%)')
print "------------------------------------------------------------------------------------------------------------------------------------------"

#Set up the constant d.
d=(1+math.sqrt(5))/2

#Set up an initial value of of variable.
x_low = 20
x_high = 60
err = 100
i=1

while abs(err)>1:
    d=(math.sqrt(5)-1)/2*(x_high-x_low)
    x1 = x_low + d
    y1 = f(x1)
    x2 = x_high - d
    y2 = f(x2)    
    y_low = f(x_low)
    y_high = f(x_high)
    x_opt = 0
    if y2 < y1:
        x_opt = x2
    else :
        x_opt = x1
    y_opt = f(x_opt)
    err=(3-math.sqrt(5))/2*abs((x_high-x_low)/x_opt)*100
    print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %(i,x_low,y_low,x2,y2,x1,y1,x_high,y_high,x_opt,y_opt,abs(err))
    #If f(x2)<f(x1), then f(x2) is the minimum and x1 becomes the new x_high for the next round.
    #Also, old x2 will be new x1.
    if y2 < y1:
        x_high = x1
        x1 = x2
    #If f(x1)<f(x2), then f(x1) is the minimum and x2 becomes the new x_low for the next round.
    #Also, old x1 will be new x2.
    else :
        x_low = x2
        x2 = x1
    i=i+1
i=i-1

print "\nAfter %i iteration, the ground state volume is predicted as %7.5f (Angstrom^3) with an error of %5.4f"\
      %(i,x_opt,err)
print "\nAfter %i iteration, the ground state energy is predicted as %7.5f (eV) with an error of %5.4f"\
      %(i,y_opt,err)

#Second, find the ground state volume and energy by Newton method of optimization.
print "\n"
print "By Newton method of optimization"
print "|%9.9s|%10.8s|%10.9s|%10.9s|%10.8s|%10.8s|%12.10s|"\
      %('iteration','x_old','f1(x_old)','f2(x_old)','x_new','f(x_new)','error(%)')
print "-----------------------------------------------------------------------------------------------------------"

#Set up an initial value of of variable.
x_old = 0
er = 100
i=1

#Define the first derivative function of f(x).
def f1(x)  :
    return float(coeff[1] + 2*coeff[2]*x + 3*coeff[3]*(x**2))

#Define the second derivative function of f(x).
def f2(x) :
    return float(2*coeff[2] + 6*coeff[3] * x)

while abs(er)>1:
    x_new=x_old - (f1(x_old))/(f2(x_old))
    er=(x_new-x_old)/x_new*100
    print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %(i,x_old,f1(x_old),f2(x_old),x_new,f(x_new),abs(er))
    x_old=x_new
    i=i+1
i=i-1

print "\nAfter %i iteration, the ground state volume is predicted as %7.5f (Angstrom^3) with an error of %5.4f"\
      %(i,x_old,abs(er))
print "\nAfter %i iteration, the ground state energy is predicted as %7.5f (eV) with an error of %5.4f"\
      %(i,f(x_old),abs(er))


#Problem 4

print "\nProblem 4"
#scipy.optimize.leastsq library is fitting a desired function from the given data by minimizing the sum of the squares of the residuals.
from scipy.optimize import leastsq

#Make the lists to array for calculation.
E_l=np.array(E_low)
v_l=np.array(v_low)

#Define the Birch-Murnaghan function with parameter list and a given V.
def murna(para,V):
    E0,V0,B0,Bp = para
    U=np.power(V0/V,2./3.)
    S=U-1.0
    X = E0 + 9.0/16.0*B0*V0*(Bp*S**3.0-(4.0*U-6.0)*S**2.0)
    return X

#Calculate the error.
def objects(para,E,V):
    error = E - murna(para,V)
    return error

#To find the values of parameter, we need to determine the initial values of parameter randomly.
#As we have to find ground state total energy and the ground state volume, determine E0 and V0 as minimum of the selected data.
x0 = [min(E_l),min(v_l),1,1]


parameters,flag= leastsq(objects, x0, args=(E_l, v_l))

#Print the coefficients of the desired curve.
print "\nCoefficients are [E0 = %.8f (eV)l V0 = %.8f (Anstrom^3)l B0 = %.8f (GPa)l B0' = %8f] "\
      %(parameters[0],parameters[1],parameters[2],parameters[3])

#Plot the graph with the fitted curve and given data.
V = np.linspace(30,60,60)
EOS = murna(parameters,V)

plt.subplot(224)
plt.plot(V,EOS,'--',color='y',label='Third-order Birch-Murnaghan EOS')
plt.plot(v_low, E_low, 'd', color='k', label='selected data')
plt.xlabel('Volume (Angstrom^3)')
plt.ylabel('Energy (eV)')
plt.title('Fig. 4')
plt.legend()
plt.tight_layout()
plt.show()
