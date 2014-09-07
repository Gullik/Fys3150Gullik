from numpy import genfromtxt
from pylab import *

my_data = genfromtxt('Results.txt', delimiter=',')

plot(my_data[:,0], my_data[:,1], label = 'Calculated solution')
plot(my_data[:,2], my_data[:,3], label = 'Exact solution')

legend(loc='upper right')

print my_data


xlabel('Steps')
ylabel('Magnitude')
#ylim(0, 1)
title('Plot')
grid(True)
savefig("test.png")
show()


