from numpy import genfromtxt
from pylab import *



# Plotting of the solutions to the first 
eigenvectors_0_5 = genfromtxt('Eigenvectors_w_r_0.5.csv', delimiter=',')
eigenvectors_1 = genfromtxt('Eigenvectors_w_r_1.csv', delimiter=',')
eigenvectors_5 = genfromtxt('Eigenvectors_w_r_5.csv', delimiter=',')

N = len(eigenvectors_0_5)
xx = linspace(0,20,N)

FigSolution = figure()
plot(xx, eigenvectors_0_5[:,0], label = 'w_r = 0.5')
hold('On')
plot(xx, eigenvectors_1[:,0], label = 'w_r = 1')
plot(xx, eigenvectors_5[:,0], label = 'w_r = 5')

legend(loc='upper right')
xlabel('rho')
ylabel('Wavefunction')
#ylim(0, 1)
title('Plot of the groundstate for all with ' + str(N + 2) + ' steps')
grid(True)
savefig("Plot of the groundstate with" + str(N +2) + ".png")
show()
hold('Off')


plot(xx, eigenvectors_0_5[:,0], label = 'Groundstate')
hold('On')
plot(xx, eigenvectors_0_5[:,1], label = 'First exited state')
plot(xx, eigenvectors_0_5[:,2], label = 'Second exited state')


legend(loc='upper right')
xlabel('rho')
ylabel('Wavefunction')
title('The three first states with w_r = 0.5 ')
grid(True)
savefig("3stateswr05.png")
show()
hold('Off')


plot(xx, eigenvectors_1[:,0], label = 'Groundstate')
hold('On')
plot(xx, eigenvectors_1[:,1], label = 'First exited state')
plot(xx, eigenvectors_1[:,2], label = 'Second exited state')


legend(loc='upper right')
xlabel('rho')
ylabel('Wavefunction')
title('The three first states with w_r = 1 ')
grid(True)
savefig("3stateswr1.png")
show()
hold('Off')


plot(xx, eigenvectors_5[:,0], label = 'Groundstate')
hold('On')
plot(xx, eigenvectors_5[:,1], label = 'First exited state')
plot(xx, eigenvectors_5[:,2], label = 'Second exited state')


legend(loc='upper right')
xlabel('rho')
ylabel('Wavefunction')
title('The three first states with w_r = 5 ')
grid(True)
savefig("3stateswr5.png")
show()
hold('Off')
