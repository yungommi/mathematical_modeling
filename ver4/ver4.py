import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import scipy.interpolate
import scipy.optimize
import math
#making figures
fig = plt.figure(figsize = (16,16))
#reading ev.txt and splicing a,c,e
Arr = np.array(map(lambda x: map(float,x), [line.split() for line in open('ev.txt')][1:]))
a, c, e = list(Arr[:,0]), list(Arr[:,1]), list(Arr[:,2])

#1
g_a, g_c = np.meshgrid(np.linspace(min(a), max(a), 200),np.linspace(min(c), max(c), 200))
#using cubic spline method
G_cubic = scipy.interpolate.griddata((a,c),e, (g_a,g_c), method = 'cubic')
ax = fig.add_subplot(2,2,1)
mesh=ax.pcolormesh(g_a,g_c, G_cubic,cmap = plt.get_cmap('rainbow'))
#making contour plot
ax.contour(g_a,g_c,G_cubic,10, linewidths =1.4,color = 'k')
plt.colorbar(mesh)
ax.set_title('Fig1');ax.set_xlabel('a($\AA$)');ax.set_ylabel('c($\AA$)')

#2
#dictionary to find V by E and find E by V
V_Edict = {}; E_Vdict = {}
#function to calculate V with given a and c and matching V - E
def putintoDict(List):
    a,c,e = List
    V= sqrt(3)*(a**2)*c/2
    V_Edict[V] = e
    E_Vdict[e] = V
#sorting list of V and slicing into list of length 6 (the number of a's is bigge than 60 so we can get more than 10points)
map(putintoDict, Arr)
Vlist = sorted(V_Edict.keys())
V = [Vlist[i:i+6] for i in range(0, len(Vlist), 6)]
#select the minimum V in sliced list of V 
minE = [min(map(lambda x: V_Edict[x],n))  for n in V]
VofminE = [E_Vdict[x] for x in minE]
#plotting the selected V's and corresponding E's and all V's and E's
ax = fig.add_subplot(2,2,2)
ax.scatter(V_Edict.keys(), V_Edict.values(),s =12, label = 'all given data', color = 'r')
ax.scatter(VofminE, minE, s = 12, label = 'selected data',color = 'b', marker = 'd')
ax.set_title('Fig.2');ax.set_xlabel('V(($\AA$)^3)');ax.set_ylabel('E(eV)');ax.legend()

#3
#polyfitted the selected points in #2 in order of 2 beacuse the plot in #2 looks like second-order graph
polycoeff = np.polyfit(VofminE, minE, 2,full = True)
linspace = np.linspace(min(VofminE), max(VofminE), 100)
poly_plot = np.polyval(polycoeff[0], linspace)
#By setting full=True in np.polyfit, we can get residuals (sum of residual error of the points)
residuals = polycoeff[1][0]
#plotting fitted polynomial and included selected points as marker
ax = fig.add_subplot(2,2,3)
ax.plot(linspace,poly_plot, label = 'fitted polynomial equation')
ax.scatter(VofminE, minE,label = 'selected  points' ,color = 'r', marker = 'd')
eq = '%.4fx^2 %.4fx +%.4f'%(polycoeff[0][0], polycoeff[0][1], polycoeff[0][2])
#putting text including coeff of polynomials and residual error
ax.text(0.95, 0.70, 'Eq: '+ eq, horizontalalignment='right',transform = ax.transAxes,
        verticalalignment='bottom',color='black', fontsize=9)
ax.text(0.95, 0.64, 'Residuals: '+ str(residuals),horizontalalignment='right',transform = ax.transAxes,
        verticalalignment='bottom',color='black', fontsize=9)
ax.set_title('Fig.3');ax.set_xlabel('V(($\AA$)^3)');ax.set_ylabel('E(eV)');ax.legend()
#defining functions doing golden section search and gradient descent
def golden_section(func, xl, xu):
    print("Golden Sectio Starting point xl = %d(minimum V) , xu = %d(maximum V)" %(xl,xu))
    from scipy.constants import golden
    ea = 10 ;noiteration = 0 ;xopt = 0
    while ea>1:
        d = (golden-1)*(xu-xl)
        x1,x2 = xl+d, xu-d
        noiteration +=1
        xu,xl,xopt = [xu,x2,x1] if func(x1)<func(x2) else [x1,xl,x2]
        ea = (2-golden)*(math.fabs((xu-xl)/(xopt+0.0)))*100
        print("Golden Section Iteration %d, ea: %.4f, x = %.4f" %(noiteration, ea, xopt))
    return xopt
def gradient_descent(a1,a0,x):
    ea ,xi,steplength, noiteration = 10000,x,20,0
    print("Gradient descent starting point x_start = %d (maximum V)" % max(Vlist))
    while ea>1:
        x = xi
        xi = x- steplength*(a1*x+a0)
        ea = math.fabs((xi-x)/(xi+0.0))*100
        steplength-=0.01 ; noiteration+=1
        print("Gradient Descent Iteration %d, ea: %.4f, x = %.4f" %(noiteration, ea, xi))
    return xi
#finding V0 using gradient descent and golden section 
V0_open = gradient_descent(2*polycoeff[0][0], polycoeff[0][1],max(Vlist))
V0 = golden_section(lambda x: np.polyval(polycoeff[0], x),min(Vlist), max(Vlist))
#find the E0 of V0 calculated by golden section search by putting V0 into fitted polynomial
E0 = np.polyval(polycoeff[0], V0)

#4
#defining function : third order Birch_Murnaghan eqation
def fittingEq(x,tpl):
    v0,e0,b0,b00 = tpl
    def VV(V):
        return np.sign(v0/V)*np.power(np.fabs(v0/V),2.0/3)    
    return e0 + (9.0/16)*b0*v0*(np.power((VV(x)-1),3)*b00 - np.power((VV(x)-1),2)*(4*VV(x)-6)    ) 
#setting range of x and defining error function
xlim = np.linspace(min(Vlist), max(Vlist), 200)
ef = lambda tpl,x,y:fittingEq(x,tpl)-y
#used optimize.leastsq which operates by making ef(error function) to zero because optimize-curve fit is much costly than this func
#coeff = initial guess of fittingEq's parameters(v0,e0,b0,b0`)
#putted b0,b0` as initial guess of 1 because I had no initial information of b0,b0`
#putted V0,E0 as initial guess because these are the values obtained by 2nd order polyfit
coeff = [V0,E0, 1, 1]
popt,pcov = scipy.optimize.leastsq(ef,coeff ,args = (VofminE, minE))
#plotting the third order Birch_Murnaghan eqation
ax = fig.add_subplot(2,2,4)
ax.plot(Vlist, fittingEq(Vlist,popt), label = 'Third-order Birch_Murnaghan')
ax.scatter(VofminE, minE, s = 12, color = 'r', marker = 'd', label = 'selected points')
ax.legend();ax.set_title('Fig.4');ax.set_xlabel('V(($\AA$)^3)');ax.set_ylabel('E(eV)')
#printing the values (V0,E0,B0,B0`) obtained by optimization
print("Result: V0(A^3): %f, E0(eV): %f, B0(Gpa): %f, B0`: %f"  %(popt[0], popt[1], popt[2], popt[3]))
#showing plot and adjusting the space between subplots
fig.subplots_adjust(hspace = 0.3, wspace = 0.3)
plt.show()
