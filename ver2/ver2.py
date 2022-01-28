#import methods to use
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.interpolate
import scipy.optimize

#before solve the problem set some factors
#read data as list(a and c are different lattice parameter, E is for energy of each data)
a=[]; c=[]; E=[]; i=1
#by using iterator 'for', read Data from 'ev.txt' and store the data in each list(a, c, E)
# i designates line of the data
for line in open("ev.txt"):
    if i==1:
        pass
    else:
        a.append(float(line.split()[0]))    #use append method to bring each data in each list
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i=i+1

#set array for each list
a_arr=np.array(a); c_arr=np.array(c); E_array=np.array(E)



#------Problem 1

#data for a contour
#setting up a regular grid of spline points.
#npoint is to divide into 500 space with 'a array' and 'c array'
npoint= 500
g_a=np.linspace(a_arr.min(), a_arr.max(), npoint)
g_c=np.linspace(c_arr.min(), c_arr.max(), npoint)
#for contour map, use meshgrid method
g_a, g_c = np.meshgrid(g_a, g_c)

#define cubic spline as 'G_cubic',and interpolate with cubic spline method by using scipy.interpolate
G_cubic=scipy.interpolate.griddata((a,c), E, (g_a, g_c), method='cubic')

#3D contour with cubic spline
#make 4 section, Fig.1 located at left up position
plt.subplot(221)
#label title
plt.title('Fig.1')
#to illustrate the height of cubic spline, use color differences
plt.pcolormesh(g_a, g_c, G_cubic, cmap = plt.get_cmap('rainbow'))
plt.colorbar()
plt.contour(g_a, g_c, G_cubic, 10, linewidths = 0.5, colors = 'k')
#label x,y axes
plt.xlabel('a (Angstrom)')
plt.ylabel('c (Angstrom)')





#------Problem 2

#define volume of HCP using a and c array
#we use the unit cell of HCP as parallelogram(the 1/3 of hexagon) so calculate the volume
#(the area of hexagon * height of the hexagonal prism)/3
V=list((math.sqrt(3)*a_arr**2*c_arr)/2)
#make a list of V and E
V_E=[]
for i in range(len(V)):
    V_E.append([V[i], E[i]])
#use sorted to sort the V_E list from the lowest volume to the highest volume
sortV=sorted(V_E)

#Choose 20 group(in the pdf, we can see the passage that 'Choose more than 10 points from the given data set, ~)
#make lists for lowest volume and lowest energy of each group 
low_V=[]
low_E=[]
#we already choose 20 group, we can group 3 data points of the given data that have similar volume(ev. txt has 63 data for each factors) 
for i in range(len(a)/3):
#reverse the list of data points by using [::-1] (we need the sorted data of lowest energy data so inverse the list because the volume is first factor)
#and choose the data point that have the lowest energy comparing other points having similar volume.
    low_E_V=min((sortV[3*i][::-1], sortV[3*i+1][::-1], sortV[3*i+2][::-1]))
#seperate volume and energy and store them in each list
    low_V.append(low_E_V[1])
    low_E.append(low_E_V[0])

#Fig.2 located at right up position
plt.subplot(222)
#label title
plt.title('Fig.2')
#plot the given data by color of blue and selected data red
plt.plot(V, E, 'd', color='b', label='all given data')
plt.plot(low_V, low_E, 'd', color='r', label='selected data')
#label x,y axes
plt.xlabel("Volume(Angstrom^3)")
plt.ylabel("Energy(eV)")
#show legend
plt.legend()




#------Problem 3

#volume energy funcion
#we choose third-order polynomial that we think has the least error polynomial regression
#there are reasons that we set the polynomial into cubic polynomial
#1. when use linear equation or quadratic polynomial, it is not accurate to fit the curve because in case of quadraic, it is symmetric but this curve is not(the dots, decreasing left side is more steep than right side)   
#2. when use linear equation or quadratic, it has limitation on some methods like newton method (limitation on using derivative)
#3. even if higher polynomial may be accrurate, in this case there are many unknown valuse which is not suitable for the polynomial fitting and can make error bigger
#so we choose third-order polynomial to fit the selected points in problem 2 ( cubic: most accurate and simple polynomials)

#define third-order polynomial
def VEf(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

#scipy.optimize.curve_fit by using VEf function to get coefficient and error(fitting the dots of volume and energy selected in problem 2, can get a, b, c, d value)
popt, pcov=scipy.optimize.curve_fit(VEf, low_V, low_E)
#standard deviation errors during fitting (by using value of 'pcov', standard deviation errors for each a, b, c, d)  
perr=np.sqrt(np.diag(pcov))
#set the factors to plot
x=np.linspace(min(V), max(V), 100)
y=VEf(x, *popt)
#Fig.3 located at left down position
plt.subplot(223)
#label title
plt.title('Fig.3')
#plot the curve of third-order polynomial and label the value of coefficient(a, b, c, d)
plt.plot(x, y, color='r', label='cubic fitting : a=%.6f$\pm$%.6f, b=%.6f$\pm$%.6f, \n\t\t c=%.6f$\pm$%.6f, d=%.6f$\pm$%.6f' \
         %(popt[0],perr[0],popt[1],perr[1],popt[2],perr[2],popt[3],perr[3]))
#plot the dots selected in problem 2
plt.plot(low_V, low_E, 'd', color='k')
#label x,y axes
plt.xlabel("Volume(Angstrom^3)")
plt.ylabel("Energy(eV)")
#show legend
plt.legend()

print('Problem 3')
print('')
print('third-order polynomial is %7.5f x^3 + %7.5f x^2 + %7.5f x + %7.5f' %(popt[0],popt[1],popt[2],popt[3]))
print('')
#Find the ground state volume V0 and E0
#use bracket method(bisection method) and open method(newton method) to get roots of thd differential of the polynomial
#before start, we set three types of function
def fx(x):
    return float(popt[0]*x**3+popt[1]*x**2+popt[2]*x+popt[3])
def fx1(x):
    return float(3*popt[0]*x**2+2*popt[1]*x+popt[2])   #differential of volme energy function
def fx2(x):
    return float(6*popt[0]*x+2*popt[1]) #second differential of volume energy function

#---bracket method (bisection method)
#choose bisection method because it is the simplest method to calculate iterations
#set the values
#tolerance = 1%
tol=1
#in bisection method, we set two points to start 
V_low=min(V)
V_high=max(V)
#p_r is to calculate percentage of relative error
p_r=0
print("bracket method (bisection method)")
print("|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %('iteration', 'x_low', 'x_high', 'x_r', 'fx_r', 'error(%)'))
print("---------------------------------------------------------------------")

#iteration 'i' is set on 0
i=0; err=100
#iteraions should be printed more than 5 lines and the error is lower than 1%
while err > tol or ((err < tol) and i < 5):
    i=i+1
    V_r=(V_low+V_high)/2.  #a system of measuring in bisection method to calculate Xr
    err=abs((p_r-V_r)/V_r)*100   #way to calculate the ea(relative error) in the textbook
    if err == 100:
        #in first iteration, there is no ea, so the space of error is blank
        print("|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|"\
              %(i, V_low, V_high, V_r, fx(V_r), ' '))
    else:
        #after first iteration, relative error will be happen
        print("|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|"\
              %(i, V_low, V_high, V_r, fx(V_r), err))
#in bisection method, this is the way to find interval that include the root
    if fx1(V_low)*fx1(V_r)<0:    # if the root is between V_low and V_r
        V_high = V_r             # V_r will be the value that has same sign one of the two
    else:
        V_low = V_r    #if the root is between V_high and V_r
    p_r=V_r     #change p_r to V_r to calculate next relative error 

#final estimated value
print('')
print("After %i iteration, the ground state Volume is estimated as %7.5f (A^3) with Energy %7.8f eV with error %.2f (%%)"\
      %(i, V_r, fx(V_r), err))
print('')



#---open method (newton method)
#choose open method which is the simplest method to calculate iterations
print("open method (newton method)")
print("|%9.9s|%16.15s|%15.12s|%13.9s|"\
      %('iteration', 'Volume (Å^3)', 'Energy (ev)', 'error (%)'))

print('-----------------------------------------------------------')

#set up an initial value
x_new=max(V)
x_old=0
i=0
err2=100

#iteraions should be printed more than 5 lines and the error is lower than 1%
while err2>tol or ((err2<tol) and i<5):
    i=i+1
    err2=abs((x_old-x_new)/x_new)*100   #way to calculate the ea(relative error) in the textbook
    if err2==100:
        #in first iteration, there is no ea, so the space of error is blank
        print("|%9.9s|%15.13s|%15.8s|%13.5s|"\
              %(i, x_new, fx(x_new), ' '))
    else:
        #after first iteration, relative error will be happen
        print("|%9.9s|%15.13s|%15.8s|%13.5s|"\
              %(i, x_new, fx(x_new), err2))
    x_old=x_new        #change x_old to x_new to calculate next relative error
    if fx2(x_new)==0:
        print('false value')  #prevent that denominator is 0 
        break
    else:
        x_new=x_new-(fx1(x_new)/fx2(x_new))  #in newton method, this is the way to find the root

#final estimated value
print('')
print("After %i iteration, the ground state Volume is estimated as %7.5f (Å^3)"\
      %(i, x_new))
print("the ground state energy is estimated as %7.8f (eV) with error = %.2f (%%)"\
      %(fx(x_new), err2))
print('')



#------Problem 4

#use third-order Birch-Murnaghan equation
#define the Birch-Murnaghan EOS
def BMeq(V, E0, B0, V0, B1):
      return E0+(9/16.)*B0*V0*(((V0/V)**(2/3.)-1)**3*B1-((V0/V)**(2/3.)-1)**2*(4*(V0/V)**(2/3.)-6))

#set the value that will be used in scipy.optimize.curve_fit
x_vals=[fx(x_new), 1, x_new, 1]
#scipy.optimize.curve_fit by using third-order Birch-Murnaghan EOS to get coefficient and error
#set the boundary to use x_vals
popt, pcov=scipy.optimize.curve_fit(BMeq, low_V, low_E, x_vals)
perr=np.sqrt(np.diag(pcov))

#print the result after optimization of Birch-Murnaghan EOS
print('')
print('Problem 4')
print('The ground state total energy : %.5f eV' %(popt[0]))
print('The ground state volume : %.5f Å^3' %(popt[2]))
print('The bulk modulus : %.5f GPa' %(popt[1]))
print('The derivative of bulk modulus : %.5f' %(popt[3]))

#Fig.4 located at right down position
plt.subplot(224)
#label the title
plt.title('Fig.4')
#plot the graph to use Birch-Murnagha EOS
plt.plot(low_V, BMeq(low_V, *popt), '-', label='Birch-Murnaghan EOS')
plt.plot(low_V, low_E, 'bo', label='Selected data')
#label x,y axes
plt.xlabel("Volume(Angstrom^3)")
plt.ylabel("Energy(eV)")
#show legend
plt.legend()
plt.show()
