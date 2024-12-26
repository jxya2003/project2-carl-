from numpy import * 
import matplotlib.pyplot as plt 

title = "Pendulum Motion for Linear Air Resistance RK4 alpha = 2.2"
showNumercialIntegration = True
showSmallAngleAprox = False

#Parameters go here
step = 0.01
t_initial = 0 #Starting values for the solution
theta_initial = pi-0.01#Starting position (for the system)
omega_initial = 0 #Starting velocity (for the system)
u_initial = [theta_initial, omega_initial]
t_final = 15

#Define the f(t,u) vector in du/dt=f(t,u)
def f(t, u):
    theta = u[0]
    omega = u[1]
    nu = 0.1
    g = 9.80665
    l = 1
    
    return [omega, -nu*omega -(g/l * sin(theta))]
    #return [omega, -(g/l * sin(theta))]

def smallAngleAproximation(t):
    g = 9.80665
    length = 1
    return theta_initial*cos(t*sqrt(g/length))
    
def getSmallAngleValues():
    thetaValuesSmallAngle = []
    for t in t_values:
        thetaValuesSmallAngle.append(smallAngleAproximation(t))
    
    
    return thetaValuesSmallAngle


def findZero(thetaValues, start):
    interator = start
    
    while True:
        if (interator == len(thetaValues)):
            return -1
        if thetaValues[interator] < 0.09 and thetaValues[interator] > -0.09:
            return interator
        interator += 1
    

def getPeriod(thetaValues):
    firstZero = findZero(thetaValues, 0)
    nextZero = findZero(thetaValues, firstZero + 100)
    thirdZero = findZero(thetaValues, nextZero + 100)
    return t_values[thirdZero] - t_values[firstZero]


def RK4MethodStep(t, u):
    k1 = f(t, u)
    k2 = f(t + step/2, add(u, multiply(step/2, k1)))
    k3 = f(t + step/2, add(u, multiply(step/2, k2)))
    k4 = f(t + step, add(u, multiply(step, k3)))
    return add(u, multiply(step/6, add(add(add(k1, multiply(2, k2)), multiply(2, k3)), k4)))

#Initialize lists to store the values of the points
t_values = [t_initial]
u_values = [u_initial]
while t_values[-1] < t_final: #Loop until you reach or pass the end point
    u_values.append(RK4MethodStep(t_values[-1], u_values[-1])) #Add the next u values (using RK4 Method)
    t_values.append(t_values[-1] + step) #Add the next t value
    
#Extract the data from the u values list of lists
theta_values = transpose(u_values)[0]
omega_values = transpose(u_values)[1]

theta_appoximate = getSmallAngleValues()

#Graph the output
plt.figure()
plt.title(title)
plt.xlabel("Time (s)")
plt.ylabel("Angle (Rad)")
print(getPeriod(theta_values))
if showNumercialIntegration:
    plt.plot(t_values, theta_values, 'bo')
if showSmallAngleAprox:
    plt.plot(t_values, theta_appoximate, 'ro')
    
plt.grid(True)
plt.show()
