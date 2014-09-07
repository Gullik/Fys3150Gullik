from numpy import genfromtxt
from pylab import *


# Plotting of the solutions to the first 
Solutions = genfromtxt('Results.csv', delimiter=',')
FigSolution = figure()
plot(Solutions[:,0], Solutions[:,1], label = 'Calculated solution')
plot(Solutions[:,2], Solutions[:,3], label = 'Exact solution')

legend(loc='upper right')

N = len(Solutions)

legend(loc='upper right')
xlabel('x')
ylabel('y')
#ylim(0, 1)
title('Plot with ' + str(N) + ' steps')
grid(True)
savefig("Plot_N_" + str(N) + ".png")
show()

# Plotting the relative error
ErrorTable = genfromtxt('ErrorTable.csv', delimiter=',')

FigRelErr = figure()
plot(log10(1/(ErrorTable[:,0]+1)), ErrorTable[:,1], label = 'Relative Error')

title('Plot of the max relative error dependant on the stepsize')
legend(loc='upper right')
xlabel('log10(1/(N+1))')
ylabel('Max Relative Error')
grid(True)
savefig("RelErrorPlot.png")

show()





