#import numpy, matplotlib.pyplot, scipy.interpolate to use various methods
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

#read data as list(a and c are for length between atoms, E is for energy of each data)
a = []; c = []; E = []; i =1
#make list for the given data (use open to bring 'ev.txt' file)
#i is to designate line of the data
for line in open("ev.txt"):
    if i == 1:
        pass
    else:
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i += 1
#with append method you can bring each data in order
#set array for each list
a_array = np.array(a); c_array = np.array(c); E_array = np.array(E)



#----------------------------------------------------------------------------------------
#1---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
npoint = 500
#npoint is to divide into 500 space with 'a array' and 'c array', which is named as 'g_a' and 'g_c' 
g_a = np.linspace(a_array.min(), a_array.max(), npoint)
g_c = np.linspace(c_array.min(), c_array.max(), npoint)
#to plot a contour map, you should use meshgrid method
g_a, g_c = np.meshgrid(g_a, g_c)

#define cubic spline as 'G_cubic'. Use scipy.interpolate.griddata to use spline method 
G_cubic = scipy.interpolate.griddata( (a, c), E, (g_a, g_c), method = 'cubic')

#show Fig. 1 at left up position
plt.subplot(221)
#labeling title and x,y axes
plt.title("Fig. 1")
plt.xlabel("lattice parameter a($\AA$)")
plt.ylabel("lattice parameter c($\AA$)")
#to illustrate the height of cubic spline on plane, use color differeces with using 'rainbow'
plt.pcolormesh(g_a, g_c, G_cubic, cmap = plt.get_cmap('rainbow'))
#show colorbar
plt.colorbar()
#plot contour
plt.contour(g_a, g_c, G_cubic, 20, linewidths = 0.5, colors = 'k')



#----------------------------------------------------------------------------------------
#2---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#define volume(the area of hexagon * height of the hexagonal prism)
V = list( 3 * ( 3**0.5 ) / 2 * a_array**2 * c_array)

#set volume in order of their sizes
selectedV = sorted(V)
#make a list for sorted volume's E in order of size of the volume
E_selectedV = [ E[ V.index(i) ] for i in selectedV]

#get the lowest energy comparing other points having similar volume(difference <= 1Angstrong^3)
while True:
    len_selectedV = len(selectedV)
    #use 'for' to review selectedV from the lowest value
    for i in range( len(selectedV) - 1 ):
        if selectedV[i+1] - selectedV[i] > 1:
            pass
            #this is the case when difference of adjacent volumes is bigger than 1
        else:
            #if difference of adjacent volumes isn't bigger than 1, erase volume with bigger energy and its energy; then break(because len of selectedV is change)
            if E_selectedV[i+1] > E_selectedV[i]:
                selectedV.pop(i+1)
                E_selectedV.pop(i+1)
            else:
                selectedV.pop(i)
                E_selectedV.pop(i)
            break
    #use len_selectedV to find if there is change in selectedV
    if len_selectedV == len(selectedV):
        #break 'while' when there is no change in selectedV
        break

#show Fig. 2 at right up position
plt.subplot(222)
#labeling title and x,y axes
plt.title("Fig. 2")
plt.xlabel("volume ($\AA$$^3$)")
plt.ylabel("energy (eV)")
#plot dots for the lowest energy comparing other points having similar volumes with blue circle
plt.plot(selectedV, E_selectedV, 'bo', label = 'the lowest energy comparing other points having similar volume')
#plot values of all volumes and energies with red dot
plt.plot(V, E, 'r.', label = 'all given data')
#show legend
plt.legend()



#----------------------------------------------------------------------------------------
#3---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#set list into array of selectedV and E_selectedV(these two list are shown in problem 2)
selectedV_array = np.array(selectedV); E_selectedV_array = np.array(E_selectedV)

#to use curve_fit import scipy.optimize
import scipy.optimize

#volume energy function(third-degree polynomial)
#set the polynomial into cubic polynomial for these reasons
#1. higher polynomial may be more accurate, but those have many unknown values to consider which is not suitable for polynomial fitting
#2. if you use quadratic polynomial(and linear equation), it is not accurate enough to fit the curve because if you consider the dots, the decreasing left side is steeper than increasing right side(quadratic is symmetric but this curve isn't)
#3. when you use quadratic polynomial, you cannot use newton method for the derivative of quadratic polynomial(linear equation)
#4. 4 or higher degree polynomial has possibility to make error bigger
#Therfore cubic would be most accurate polynomial in simple polynomials
def VEfunc(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

#get coefficient and error of volume energy function
#you can get a, b, c, d value with scipy.optimize.curve_fit, fitting the dots of volume and energy chosen in problem 2
popt, pcov = scipy.optimize.curve_fit(VEfunc, selectedV_array, E_selectedV_array)
#standard deviation errors for each a, b, c, d value using 'pcov' value(covariance)
perr = np.sqrt(np.diag(pcov))

print "=================================================================================="
print "problem 3."
#with caculated a, b, c, d values you can get the cubic polynomial; show all numbers of value(not using %.xf to show only to x decimal); u means using unicode(find unicode from www.fileformat.info)
print "E(V) = ", popt[0], u"( \u00B1", perr[0], u") * x\u00B3+ ", popt[1], u"( \u00B1", perr[1], u")  * x\u00B2+ ", popt[2], u"( \u00B1", perr[2], ")  * x + ", popt[3], u"( \u00B1", perr[3], ")"
print "----------------------------------------------------------------------------------"

#show Fig. 3 at left down position
plt.subplot(223)
#labeling title and x,y axes
plt.title("Fig. 3")
plt.xlabel("volume ($\AA$$^3$)")
plt.ylabel("energy (eV)")
#plot the curve of cubic polynomial and label it with significant figures that can show the value of coefficient(a, b, c, d)
plt.plot(selectedV_array, VEfunc(selectedV_array, *popt), 'r-', label = 'fit with cubic : a=%.6f$\pm$%.6f, b=%.6f$\pm$%.6f,\n\t\t c=%.6f$\pm$%.6f, d=%.6f$\pm$%.6f' %(popt[0],perr[0],popt[1],perr[1],popt[2],perr[2], popt[3], perr[3]))
#plot the dots chosen in problem 2
plt.plot(selectedV_array, E_selectedV_array, 'bo')
#show legend
plt.legend()

#caculate the ground state volume V_0 and energy E_0
#caculate with bracket method(bisection method) and open method(newton method) by getting roots of the differential of cubic polynomial
#differential of volume energy function(VEfunc)
def difVEfunc(x):
    return 3 * popt[0] * x**2 + 2 * popt[1] * (x) + popt[2]

#bracket method(bisection method)-----------------------------------------------------------------------
tol = 1     #tolerance = 1%
#set two points to start the bisection method; the dot behind is to designate value as float
V_low = 120.
V_high = 140.
#set previous V_r to caculate percentage relative error
preV_r = 0

print "bracket method(bisection method)\n"
#setting significant numbers of each value and space of each value
print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|" %('Iteration', 'V_low', 'V_high', 'V_r', 'E\'(V_r)', 'error (%)')
print "+---------+----------+----------+----------+----------+------------+"

i = 0
err = 100
#error should be under 1% and the lines should be at least 5 lines
while err > tol or ( ( err < tol ) and i < 5):
    i += 1
    #V_r is the estimated value; the dot behind is to designate value as float
    V_r = ( V_low + V_high ) / 2.
    #err is the percent relative error
    err = abs( ( preV_r - V_r ) / V_r ) * 100
    #because there is no percent relative error for the first iteration, make it a blank
    if err == 100:
        print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|" %(i, V_low, V_high, V_r, difVEfunc(V_r), ' ')
    else:
        print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|" %(i, V_low, V_high, V_r, difVEfunc(V_r), err)
    #progressing bisection method, deciding whether the root is in interval (V_low, V_r) or (V_r, V_high)
    if difVEfunc(V_low) * difVEfunc(V_r) <0:
        #if root is in interval (V_low, V_r)
        V_high = V_r
    else:
        V_low = V_r
    #designate previous value(preV_r) to caculate errors  
    preV_r = V_r

#print the final estimated value
print u"\nAfter %i iternation, the E\u2080 is predicted as %.5f (eV) with V\u2080 %.3f (\u00c5\u00B3)" %(i, VEfunc(V_r, *popt), V_r)
print "----------------------------------------------------------------------------------"

#open method(newton method)--------------------------------------------------------------------------
tol = 1     #tolerance = 1%
#set the starting point of newton method, the dot behind is to designate value as float
V_r = 150.
#set previous V_r to caculate percentage relative error
preV_r = 0

#second differential of volume energy function(VEfunc)
def dif2VEfunc(x):
    return 6 * popt[0] * x + 2 * popt[1]

print "open method(Newton method)\n"
print "|%9.9s|%10.8s|%10.8s|%12.10s|" %('Iteration', 'V_r', 'E\'(V_r)', 'error (%)')
print "+---------+----------+----------+------------+"

i=0
err=100
#error should be under 1% and the lines should be at least 5 lines
while err > tol or ( ( err < tol ) and i < 5):
    i += 1
    dif2E_Vr = dif2VEfunc(V_r)
    difE_Vr = difVEfunc(V_r)
     #err is the percetage relative error
    err = abs( ( preV_r - V_r ) / V_r ) * 100
     #because there is no percent relative error for the first iteration, make it a blank
    if err==100:
        print "|%9.9s|%10.8s|%10.8s|%12.7s|" %(i, V_r, difE_Vr, ' ')
    else:
        print "|%9.9s|%10.8s|%10.8s|%12.7s|" %(i, V_r, difE_Vr, err)
        
    preV_r = V_r
    #the equation for the newton method is x_n+1 = x_n - ( f(x_n) / f'(x_n) ) and therfore f'(x_n) should not be zero. for this equation, dif2E_Vr should not be zero
    if dif2E_Vr==0:
        print 'Can\'t solve'
        break
    else:
        #the equation for the newton method
        V_r = V_r - ( ( difE_Vr ) / ( dif2E_Vr ) )

#print the final estimated value
print u"\nAfter %i iternation, the E\u2080 is predicted as %.5f (eV) with V\u2080 %.5f (\u00c5\u00B3)" %(i, VEfunc(V_r, *popt), V_r)



#----------------------------------------------------------------------------------------
#4---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#third-order Birch-murnaghan equation of state (BMfunc) (in the README file) (B_1 means B'_0)
def BMfunc(V, E_0, B_0, V_0, B_1):
    return E_0 + 9 / 16. * B_0 * V_0 * ( ( ( V_0 / V )**( 2 / 3. ) - 1 )**3 * B_1 - ( ( V_0 / V )**( 2 / 3. ) - 1 )**2 * ( 4 * ( V_0 / V )**( 2/3. ) - 6 ) )

#get coefficient and error of third-order Birch-Murnaghan EOS(BMfunc)
#caculate the V_0, E_0, B_0, B_1(B'_0) value with curve fit; set boundary : E_0 (-inf, 0), B_0 (-inf, inf), V_0 (120, 140), B'_0 (-inf, inf)
popt, pcov = scipy.optimize.curve_fit(BMfunc, selectedV_array, E_selectedV_array, bounds = ( [-np.inf, -np.inf, 120, - np.inf], [0, np.inf, 140, np.inf] ) )  #, maxfev=100000
perr = np.sqrt(np.diag(pcov))

print "=================================================================================="
print "problem 4."
print u"\tE\u2080 = %.5f \u00B1 %.5f" %(popt[0], perr[0])
print u"\tV\u2080 = %.5f \u00B1 %.5f" %(popt[2], perr[2])
#popt[1] is B_0 value in the BMfunc
print u"\tB\u2080 = %.5f \u00B1 %.5f" %(popt[1], perr[1])
#popt[3] is B_1 value in the BMfunc which is B'_0
print u"\tB'\u2080 = %.5f \u00B1 %.5f" %(popt[3], perr[3])


#show Fig. 4 at right down position
plt.subplot(224)
#labeling title and x,y axes
plt.title("Fig. 4")
plt.xlabel("volume ($\AA$$^3$)")
plt.ylabel("energy (eV)")
#plot the graph for Birch-Murnagha EOS
plt.plot(selectedV_array, BMfunc(selectedV_array, *popt), 'r-', label='Birch-Murnaghan EOS(fitting curve)')
plt.plot(selectedV, E_selectedV, 'bo', label = 'selected points')
#show legend
plt.legend()
plt.show()
