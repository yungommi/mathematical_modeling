import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

# Set up a points from 'ev.txt'

a=[];c=[];E=[]
i=1
for line in open('ev.txt'):
    if i==1:
        pass
    else:
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i = i + 1


#Fig.1

npoint=200
g_x1 = np.linspace(min(a), max(a), npoint)
g_x2 = np.linspace(min(c), max(c), npoint)
g_x1, g_x2 = np.meshgrid(g_x1, g_x2)


# Interpolate 
G_cubic = scipy.interpolate.griddata((a,c), E, (g_x1, g_x2), method='cubic')


# Contour with cubic spline method
plt.subplot(221)
plt.title('Fig.1')
plt.pcolormesh(g_x1, g_x2, G_cubic, cmap = plt.get_cmap('rainbow'))
plt.colorbar()
plt.contour(g_x1, g_x2, G_cubic, 10, linewidths = 0.5, colors = 'k')
plt.xlabel('a (Angstrom)')
plt.ylabel('c (Angstrom)')




#Fig.2

V=[]

for j in range(len(a)):
    V.append(3*(3**0.5)*a[j]*a[j]*c[j]/2)

All_given_data=[]
Selected_E=[]
Selected_V=[]
for k in range(len(a)):
    All_given_data.append([V[k],E[k]])
    
All_given_data=sorted(All_given_data)

for l in range(len(a)/5):
    if min(All_given_data[5*l][1],All_given_data[5*l+1][1],All_given_data[5*l+2][1],All_given_data[5*l+3][1],All_given_data[5*l+4][1])==All_given_data[5*l][1]:
        Selected_E.append(All_given_data[5*l][1])
        Selected_V.append(All_given_data[5*l][0])
                            
    elif min(All_given_data[5*l][1],All_given_data[5*l+1][1],All_given_data[5*l+2][1],All_given_data[5*l+3][1],All_given_data[5*l+4][1])==All_given_data[5*l+1][1]:
        Selected_E.append(All_given_data[5*l+1][1])
        Selected_V.append(All_given_data[5*l+1][0])
                             
    elif min(All_given_data[5*l][1],All_given_data[5*l+1][1],All_given_data[5*l+2][1],All_given_data[5*l+3][1],All_given_data[5*l+4][1])==All_given_data[5*l+2][1]:
        Selected_E.append(All_given_data[5*l+2][1])
        Selected_V.append(All_given_data[5*l+2][0])
                             
    elif min(All_given_data[5*l][1],All_given_data[5*l+1][1],All_given_data[5*l+2][1],All_given_data[5*l+3][1],All_given_data[5*l+4][1])==All_given_data[5*l+3][1]:
        Selected_E.append(All_given_data[5*l+3][1])
        Selected_V.append(All_given_data[5*l+3][0])
                             
    else :
        Selected_E.append(All_given_data[5*l+4][1])
        Selected_V.append(All_given_data[5*l+4][0])  
    
plt.subplot(222)
plt.title('Fig.2')
plt.plot(V,E,'d', color='g', label = "All given data")
plt.plot(Selected_V,Selected_E,'d',color='r', label ="Selected data")
plt.xlabel("Volume")
plt.ylabel("eV")
plt.legend()

#Fig3

n=len(Selected_V)

sum_of_v=0; sum_of_v_2=0; sum_of_v_3=0; sum_of_v_4=0; sum_of_v_5 = 0; sum_of_v_6 = 0

for i in range(n):
    sum_of_v = sum_of_v + Selected_V[i]
    sum_of_v_2 = sum_of_v_2 + Selected_V[i]**2
    sum_of_v_3 = sum_of_v_3 + Selected_V[i]**3
    sum_of_v_4 = sum_of_v_4 + Selected_V[i]**4
    sum_of_v_5 = sum_of_v_5 + Selected_V[i]**5
    sum_of_v_6 = sum_of_v_6 + Selected_V[i]**6
sum_v=[n, sum_of_v, sum_of_v_2, sum_of_v_3, sum_of_v_4, sum_of_v_5, sum_of_v_6]

sum_of_E=0; sum_of_vE=0; sum_of_v_2E=0; sum_of_v_3E = 0

for i in range(n):
    sum_of_E = sum_of_E + Selected_E[i]
    sum_of_vE = sum_of_vE + Selected_E[i] * Selected_V[i]
    sum_of_v_2E = sum_of_v_2E + Selected_E[i] * (Selected_V[i]**2)
    sum_of_v_3E = sum_of_v_3E + Selected_E[i] * (Selected_V[i]**3)
E_M=np.zeros((4,1))
E_M[0,0] = sum_of_E
E_M[1,0] = sum_of_vE
E_M[2,0] = sum_of_v_2E
E_M[3,0] = sum_of_v_3E
    
v_M = np.zeros((4,4))
i = 0
while i < 4 :
    v_M[i,0] = sum_v[i]
    v_M[i,1] = sum_v[i+1]
    v_M[i,2] = sum_v[i+2]
    v_M[i,3] = sum_v[i+3]
    i = i+1

v_inv = np.linalg.inv(v_M)
coeff = np.matmul(v_inv, E_M)
coeff_array = np.array(coeff)

x=np.linspace(100,200,200)
y = coeff_array[0] + coeff_array[1]*x + coeff_array[2]*(x**2) + coeff_array[3]*(x**3)

print "# Figure.3"
print "\nThe polynomial fitting result is: y= %.7f + %.7f x1 + %.7f x2 + %.7f x3"\
      %(coeff[0],coeff[1],coeff[2],coeff[3])

plt.subplot(223)
plt.title('Fig.3')
plt.plot(Selected_V,Selected_E,'d',color='r', label ="Selected data")
plt.plot(x,y,'-',color='k', label = 'curve')
plt.ylim([-16.5,-14.9])
plt.xlabel("Volume")
plt.ylabel("eV")
plt.legend()


#Error
E_mean = sum(Selected_E)/len(Selected_E)
St=0
for i in range(n):
    stz = (Selected_E[i]-E_mean)**2
    St = St + stz

Sr = 0

for i in range(n):
    srz = (Selected_E[i]-(coeff_array[0] + coeff_array[1]*Selected_V[i] + coeff_array[2]*(Selected_V[i]**2)+coeff_array[3]*(Selected_V[i]**3)))**2
    Sr = Sr +srz

R_square = (St-Sr)/St
print 'R_square is', R_square
print ""

# bisection
target_V = 131.099485
tol = 0.001
def Y(x):
    return coeff_array[0] + coeff_array[1]*x + coeff_array[2]*(x**2) +coeff_array[3]*(x**3) 
def Yprime(x):
    return coeff_array[1] + 2*coeff_array[2]*x +3*coeff_array[3]*(x**2) 
x_low = 125.
x_high = 135.
x_r = (x_low +x_high)/2
print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %('iteration','x_low','x_high','x_r','f(x_r)','error (%)')
print "--------------------------------------------------------------------"

i = 1.
err = 100.

while abs(err) > tol :
    f_xr = Yprime(x_r)
    if (Yprime(x_low) - Yprime(target_V))*(f_xr - Yprime(target_V)) < 0:
        x_high = x_r
    else :
        x_low = x_r
    x_r = (x_low+x_high)/2
    err = ((target_V - x_r)/target_V)*100    
    print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|"\
%(i,x_low,x_high,x_r,f_xr,err)
    plt.plot(x_r,f_xr,'d',label = i)
    i = i+1

print "\nAfter %i iteration, the volume is predicted as %7.5f (Å3) with Energy %7.8f eV"\
      %(i,x_r,Y(x_r))
print "\n"


# newton method
tol = 0.01
i=1
print "|%9.9s|%10.8s|%10.8s|%12.10s|"\
      %('iteration','lastX','nextX','error (%)')
print "----------------------------------------------"
def derivative(f, x, h):
      return (f(x+h) - f(x-h)) / (2.0*h)
def Yprime(x):
	return coeff_array[1] + coeff_array[2]*x*2 + coeff_array[3]*3*(x**2)
def solve(f, x0, h):
    lastX = x0
    nextX = lastX + 10* h
    err = 100
    i=1
    while (abs(err) > tol):  
        newY = f(nextX)                     
        err = (nextX-lastX)/nextX*100
        lastX = nextX
        nextX = lastX - newY / derivative(f, lastX, h)
        print "|%9.9s|%10.8s|%10.8s|%12.10s|"\
      %(i,lastX,nextX,err)
        i=i+1
    return nextX


xFound = solve(Yprime, 150, 0.01)
print "\nthe ground state volume(V) =", xFound
print "the ground state energy(eV) =", Y(xFound)
print ""

#Fig.4

#GSE = E0, GSV = V0, GSB = B0, GSBD = B0'

E0 = -16.2348
V0 = 130.96155
B0 = 0.65426579
B0P = 3.40847192

GSE = float(E0)
GSV = float(V0)
print "\n# Figure.4"


def f(x,a,b,c,d):
    return a+9*b*c/16*((((V0/x)**(2.0/3.0))-1)*(((V0/x)**(2.0/3.0))-1)*(((V0/x)**(2.0/3.0))-1)*d+(((V0/x)**(2.0/3.0))-1)*(((TC[1]/x)**(2.0/3.0))-1)*(6-4*((TC[1]/x)**(2.0/3.0))))

TC=[E0,V0,B0,B0P]
y=f(x,TC[0],TC[1],TC[2],TC[3])
popt,pcov=scipy.optimize.curve_fit(f,Selected_V,Selected_E)


print ""
print 'The ground state total energy is :',popt[0],'eV'
print 'The ground state volume is :',popt[1],'Å3'
print 'The bulk modulus is :',popt[2],'GPa'
print 'The derivative of bulk modulus is :',popt[3]


plt.subplot(224)
plt.title('Fig.4')
plt.plot(x,f(x,*popt),'-', color='b', label = 'EOS')
plt.plot(Selected_V,Selected_E,'d',color='r', label ="Selected data")
plt.xlabel("Volume")
plt.ylabel("EOS")
plt.legend()
plt.show()
