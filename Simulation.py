
import math
from math import pi as pi
from math import sin as sin
from math import cos as cos
from math import sqrt as sqrt
import csv

#data file

writer = csv.writer(open('results.csv', 'w', newline=''))
writer.writerow(["X", "Y"])
# define parameters
h1 = 0.6
m_cw = 2
m_b = 0.3375
m_p = 0.025
l1 = 0.2
l2 = 0.8
l3 = 0.2
l = l1 + l2
l_cm = l2 - l / 2
Ib = (1 / 12) * m_b * l**2 + m_b * l_cm**2
alphar = 0.95*pi
# define constants
Cd = 0.42
g = 9.81
r = 0.02
A = float(pi * r ** 2)
rho = 1.225
dt1 = 0.0001
dt2 = 0.01
dtheta = 0.1

# Record X,Y for max range
# stage one
# define initial conditions
theta = float(math.acos((h1-l3)/l2))
thetad = float(0)
alpha = float(0)
alphad = float(0)
thetadd = float(0)
alphadd = float(0)

# define functions
def differential_equations_1():
    global thetadd  # sets out the second time derivative of theta
    global alphadd  # sets out the second time derivative of alpha
    thetadd = float(((m_cw * l1 - m_b * l_cm - m_p * l3) * g * sin(theta) + m_p * l2 * (
                g * sin(alpha) * cos(theta - alpha) - l3 * alphad ** 2 * sin(theta - alpha) - l2 * thetad ** 2 * cos(theta - alpha) * sin(theta - alpha))) / (m_cw * l1 ** 2 - Ib + m_p * l2 ** 2 * (sin(theta - alpha)) ** 2))
    alphadd = float((l2 / l3)*(thetad ** 2 * sin(theta - alpha) - thetadd * cos(theta - alpha)))
    return


def getx():
    x = -(l2 * sin(theta) + l3 * sin(alpha))  # return the x coordinate of the projectile for a given theta and alpha
    return x


def gety():
    y = (-l2 * cos(theta) - l3 * cos(alpha))  # return the y coordinate of the projectile for a given theta and alpha
    return y




while alpha < (alphar):
    writer.writerow([str(getx()), str(gety())])
    differential_equations_1()
    theta = theta + thetad * dt1
    thetad = thetad + thetadd * dt1
    alpha = alpha + alphad * dt1
    alphad = alphad + alphadd * dt1

# stage two
# initial conditions
v0 = sqrt(l2**2*thetad**2 + 2*l2*l3*thetad*alphad*cos(theta - alpha)+l3**2*alphad**2)
phi = float(pi-alpha)
x = getx()
y = gety()
xd = v0*cos(phi)
yd = v0*sin(phi)
xdd = float(0)
ydd = float(0)




# define functions
def differential_equations_2():
    global xdd
    global ydd
    xdd = -((rho/(2*m_p))*sqrt(xd ** 2 + yd ** 2)*Cd*A*xd)
    ydd = -g-((rho/(2*m_p))*sqrt(xd ** 2 + yd ** 2)*Cd*A*yd)
    return

while (y>-h1):
    writer.writerow([str(x), str(y)])
    differential_equations_2()
    x = x + xd * dt1
    xd = xd + xdd * dt1
    y = y + yd * dt1
    yd = yd + ydd * dt1

print('max range is: ' + str(x))

# Record theta(0) vs range
theta0 = float(math.acos((h1-l3)/l2))
writer.writerow(["Theta(0)", "Range"])
# stage one
# define initial conditions
while theta0<(pi/2):
    theta = theta0
    thetad = float(0)
    alpha = float(0)
    alphad = float(0)
    thetadd = float(0)
    alphadd = float(0)






    while alpha < (alphar):

        differential_equations_1()
        theta = theta + thetad * dt2
        thetad = thetad + thetadd * dt2
        alpha = alpha + alphad * dt2
        alphad = alphad + alphadd * dt2

    # stage two
    # initial conditions
    v0 = sqrt(l2**2*thetad**2 + 2*l2*l3*thetad*alphad*cos(theta - alpha)+l3**2*alphad**2)
    phi = float(pi-alpha)
    x = getx()
    y = gety()
    xd = v0*cos(phi)
    yd = v0*sin(phi)
    xdd = float(0)
    ydd = float(0)

    # define functions

    while (y>-h1):

        differential_equations_2()
        x = x + xd * dt2
        xd = xd + xdd * dt2
        y = y + yd * dt2
        yd = yd + ydd * dt2

    writer.writerow([str(theta0), str(x)])
    theta0 += dtheta