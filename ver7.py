import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate


## Question 1
#Need group for gather data in text. so make empty group
a=[]; c=[]; E=[]; V=[] ; i=1 ; dV=[] ; dE=[]
#open data in text
#use if to repeat and set 'i=1' not to gather first line in text. first line was not a data.
#use append for add datas in group.

for line in open("ev.txt"):
    if i == 1:
        pass
    else:
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i=i+1

npoint = 200
#define a and c range for graph.
a_plot = np.linspace(min(a), max(a),npoint)
c_plot = np.linspace(min(c), max(c),npoint)
#a and c are not related, they only effects to E. so use meshgrid.
a_plot, c_plot = np.meshgrid(a_plot, c_plot)

## Interpolate (cubic graph used)
G_cubic = scipy.interpolate.griddata( (a, c), E, (a_plot, c_plot), method = 'cubic')

#to plot 1,2,3,4 in one side, use subplot.
plt.subplot(221)
plt.title('Fig. 1')
plt.xlabel("a($\AA$)", fontsize = 7)
plt.ylabel("c($\AA$)", fontsize = 7)
#make contour line to see E
plt.pcolormesh(a_plot, c_plot, G_cubic, cmap = plt.get_cmap('rainbow'))
plt.colorbar()
plt.contour(a_plot, c_plot, G_cubic, 10, linewidths = 0.5, colors = 'k')


## Question 2
#define data V, using a and c datas in the text
#use formula to define hcp volume that base is parallelogram
from numpy import sqrt
j=0
while j<63:
    V.append(2*sqrt(3)/4*float(a[j])**2*float(c[j]))
    j = j + 1

#Gather similar volumes in same groups and find minumum E in each group. 
#make group. Similar volumes will be in same groups.
V0=[]; V1=[]; V2=[]; V3=[]; V4=[]; V5=[]; V6=[]; V7=[]; V8=[]; V9=[];
E0=[]; E1=[]; E2=[]; E3=[]; E4=[]; E5=[]; E6=[]; E7=[]; E8=[]; E9=[];
z=0

#divide ranges to gather similar volumes.V and E is related so define Vx and Ex together.
while z<63:
    if 35<V[z]<36.8:
        V0.append(V[z])
        E0.append(E[z])
    elif 36.8<V[z]<38.6:
        V1.append(V[z])
        E1.append(E[z])
    elif 38.6<V[z]<40.4:
        V2.append(V[z])
        E2.append(E[z])
    elif 40.4<V[z]<42.2:
        V3.append(V[z])
        E3.append(E[z])
    elif 42.2<V[z]<44:
        V4.append(V[z])
        E4.append(E[z])
    elif 44<V[z]<45.8:
        V5.append(V[z])
        E5.append(E[z])
    elif 45.8<V[z]<47.6:
        V6.append(V[z])
        E6.append(E[z])
    elif 47.6<V[z]<49.4:
        V7.append(V[z])
        E7.append(E[z])
    elif 49.4<V[z]<51.2:
        V8.append(V[z])
        E8.append(E[z])
    elif 51.2<V[z]<53:
        V9.append(V[z])
        E9.append(E[z])
    else:
        pass
    z=z+1

#define minimum value in a group. and find all minimum (V,E) in every group. 
#compare values in a group and repeat with using while. when find minumum value, use break to finish repeating .
q0=0
while q0<len(E0)-1:
    if E0[q0] == min(E0):
        q0=q0
        break
    else:
        pass
    q0=q0+1

q1=0
while q1<len(E1)-1:
    if E1[q1] == min(E1):
        q1=q1
        break
    else:
        pass
    q1=q1+1

q2=0
while q2<len(E2)-1:
    if E2[q2] == min(E2):
        q2=q2
        break
    else:
        pass
    q2=q2+1

q3=0
while q3<len(E3)-1:
    if E3[q3] == min(E3):
        q3=q3
        break
    else:
        pass
    q3=q3+1

q4=0
while q4<len(E4)-1:
    if E4[q4] == min(E4):
        q4=q4
        break
    else:
        pass
    q4=q4+1

q5=0
while q5<len(E5)-1:
    if E5[q5] == min(E5):
        q5=q5
        break
    else:
        pass
    q5=q5+1

q6=0
while q6<len(E6)-1:
    if E6[q6] == min(E6):
        q6=q6
        break
    else:
        pass
    q6=q6+1

q7=0
while q7<len(E7)-1:
    if E7[q7] == min(E7):
        q7=q7
        break
    else:
        pass
    q7=q7+1

q8=0
while q8<len(E8)-1:
    if E8[q8] == min(E8):
        q8=q8
        break
    else:
        pass
    q8=q8+1

q9=0
while q9<len(E9)-1:
    if E9[q9] == min(E9):
        q9=q9
        break
    else:
        pass
    q9=q9+1

#plot all given datas 
plt.subplot(222)

plt.plot(V,E,'d',color='r', label = 'all given data')


#plot selected minimum datas in different color.
plt.plot(V0[q0],E0[q0],'d',color='b')
plt.plot(V1[q1],E1[q1],'d',color='b')
plt.plot(V2[q2],E2[q2],'d',color='b')
plt.plot(V3[q3],E3[q3],'d',color='b')
plt.plot(V4[q4],E4[q4],'d',color='b')
plt.plot(V5[q5],E5[q5],'d',color='b')
plt.plot(V6[q6],E6[q6],'d',color='b')
plt.plot(V7[q7],E7[q7],'d',color='b')
plt.plot(V8[q8],E8[q8],'d',color='b')
plt.plot(V9[q9],E9[q9],'d',color='b', label = 'selected')
plt.title('Fig.2')
plt.xlabel('Volume($\AA^3$)', fontsize = 7)
plt.ylabel('Energy(eV)', fontsize = 7)
plt.legend()

## Question 3

# declare the initial value
x = [V0[q0], V1[q1], V2[q2], V3[q3], V4[q4], V5[q5], V6[q6], V7[q7], V8[q8], V9[q9]]
y_true = [E0[q0], E1[q1], E2[q2], E3[q3], E4[q4], E5[q5], E6[q6], E7[q7], E8[q8], E9[q9]]
# calculate the coefficent
h_11 =(y_true[4] - y_true[0])/(x[4]-x[0])
h_12 =(y_true[9] - y_true[4])/(x[9]-x[4])
h_2=(h_12 - h_11)/(x[9]-x[0])

# make the data to plot graph
x_plot = np.linspace(35, 53, 101)
y_plot_inter   = y_true[0] + h_11*(x_plot-x[0])
y_plot_inter_2 = y_true[0] + h_11*(x_plot-x[0]) + h_2*(x_plot-x[0])*(x_plot - x[4])

# plot the graph
plt.subplot(223)
plt.plot(x,y_true,'d', color = 'r', label = 'selected')
plt.plot(x_plot,y_plot_inter_2,'--',color = 'b', label='quadratic interpolation')
plt.ylim([-16.4,-15.8])
plt.xlabel('Volume($\AA^3$)', fontsize = 7)
plt.ylabel('Energy(eV)', fontsize = 7)
plt.title('Fig.3')
plt.legend()




## Bisection Method
# approach : The situation is not moderate to bisection method. So, we coordinate the criteria to choose the new x value.
# the criteria is that. when the f(V_low)<f(V_high) select new x value as V_high. if not select new x value as V_low 
# We think that is possible because the fittting curve is fitted by quadratic interpolation. So the fitting curve is symmetry.
def fny_plot_inter_2(x_plot):
  
    return y_true[0] + h_11*(x_plot-x[0]) + h_2*(x_plot-x[0])*(x_plot - x[4])


# we choose the target value refering to the graph of prob.#2
target_V = 44.6
target_E = fny_plot_inter_2(44.6)



V_low = 35
V_high = 50


print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.10s|"\
      %('Iteration','V_low','V_high','V_r','f(V_r)','error(%)')
print"--------------------------------------------------------------------"

i = 1
err=100

while err >= 0.002 :
        V_r = (V_low + V_high)/2.0
        f_Vr = fny_plot_inter_2(V_r)
        

        err = abs((f_Vr - target_E) / target_E * 100)

        print "|%9.9s|%10.8s|%10.8s|%10.8s|%10.8s|%12.7s|"\
        %(i,V_low,V_high,V_r,f_Vr,err)
        if fny_plot_inter_2(V_low) < fny_plot_inter_2(V_high) :
                    V_high = V_r
        
        else :
                    V_low = V_r

        i=i+1

print"\nAfter %i iternation, the ground state volume is predicted as %7.4f($\AA^3$) with the ground state energy %7.4f(eV)"\
%(i-1, V_r, f_Vr)

fny11_plot_inter_2=y_true[0] + h_11*(x_plot-x[0]) + h_2*(x_plot-x[0])*(x_plot - x[4]) - target_E

##Open Method
#approach : we make the moderate form to using open method by (fitting curve - target E)


#define the derivative of fiiting curve
def fny1_plot_inter_2(x_plot):

    return 2*h_2*x_plot + h_11 - h_2*(x[4] + x[0])

#define the target function by (fitting curve - target E)
def fny2_plot_inter_2(x_plot):

    return y_true[0] + h_11*(x_plot-x[0]) + h_2*(x_plot-x[0])*(x_plot - x[4]) - target_E

# we choose the initiate value refering to the graph of prob.#2
V_initiate = 42

print "|%9.9s|%10.8s|%10.8s|%12.10s|"\
      %('Iteration','V_i','f(i)','error(%)')
print"----------------------------------------------"

i=1
err=100

while err >= 0.002:
        V_r = V_initiate
        f_Vr = fny_plot_inter_2(V_r)
        f_Vrt = fny2_plot_inter_2(V_r)


        err = abs((f_Vr - target_E) / target_E * 100)

        print "|%9.9s|%10.8s|%10.8s|%12.7s|"\
        %(i,V_r,f_Vr,err)
        if  f_Vrt> 0 :
                    i=i+1
                    V_initiate = V_initiate - f_Vrt/fny1_plot_inter_2(V_r)    
        else :
                    break

print"\nAfter %i iternation, the ground state volume is predicted as %7.4f($\AA^3$) with the ground state energy %7.4f(eV)"\
%(i-1, V_r, f_Vr)


## Question 4

from scipy.optimize import leastsq
vols = np.array(x)
energies = np.array(y_true)
def Birch_Murnaghan(parameters, vol):
    E0, B0, BP, V0 = parameters
    E = E0 + B0*V0*(((V0/vol)**(2.0/3.0)-1.0)**3.0*BP-((V0/vol)**(2.0/3.0)-1.0)**2.0*((V0/vol)**(2.0/3.0)*4.0-6.0))*9.0/16.0
    return E

def objective(pars, y, z):
    #minimize the function
    err =  y - Birch_Murnaghan(pars, z)
    return err

x0 = [ -16.0, 0.54, 2.0, 44] #initial guess of parameters

plsq = leastsq(objective, x0, args=(energies, vols))
print('')
print 'Fitted parameters = {0}'.format(plsq[0])

plt.plot(vols,energies, 'ro')

# define z and y
z = np.linspace(min(vols), max(vols), 50)
y = Birch_Murnaghan(plsq[0], z)

# plot curve
plt.subplot(224)
plt.title('Fig.4')
plt.plot(x,y_true,'d', color = 'b', label = 'selected')
plt.plot(z, y, 'k-', label='Curve fit')
plt.xlabel('Volume($\AA^3$)', fontsize = 7)
plt.ylabel('Energy(eV)', fontsize = 7)
plt.legend()
plt.show()
