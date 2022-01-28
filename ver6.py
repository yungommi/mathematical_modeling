import numpy 
import scipy.interpolate
import matplotlib.pyplot as plt

#Problem1

#making list about a, c, E with data ev.txt

a = []
c = []
E = []
i = 1

for line in open("ev.txt"):
    if i == 1:
        pass
    else:
        a.append(float(line.split()[0]))
        c.append(float(line.split()[1]))
        E.append(float(line.split()[2]))
    i = i + 1

#count the number of data     
point_a = len(a)
point_c = len(c)

#make spline point
axis_a = numpy.linspace(min(a), max(a), point_a)
axis_c = numpy.linspace(min(c), max(c), point_c)
axis_a, axis_c = numpy.meshgrid(axis_a, axis_c)

#using scipy.interpolate for cubic spline method
_E_ = scipy.interpolate.griddata((a, c), E, (axis_a, axis_c), method='cubic')

#plotting 3D contour figure 
plt.subplot(221)
plt.pcolormesh(axis_a, axis_c, _E_, cmap = plt.get_cmap('rainbow'))
plt.contour(axis_a, axis_c, _E_,linewidths = 0.5, color='k')
plt.colorbar()


#label axes and title

plt.title('Fig. 1')
plt.xlabel('a (Angstrom)')
plt.ylabel('c (Angstrom)')
plt.legend()


#Problem2


#Problem3


#Problem4


