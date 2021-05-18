'''
solving gps equation using Newton-Raphson Method (Multivariate)
that is explained well here : http://fourier.eng.hmc.edu/e176/lectures/ch2/node6.html
'''

import numpy as np #import numpy library to be abe to use matrices operations

#instialsing the variables that i well use
x1 = 25700
y1 = 6800
x2 = 23030
y2 = 13300
x3 = 17100
y3 = 20375
v = 300
tl1 = 105.10200
ts1 = 36.532944917384
tl2 = 105.50200
ts2 = 37.8525643607006
tl3 = 108.00200
ts3 = 38.6962494012455
#the matrix of answer and i gave it 0 as initial value
xye=np.array([[0],[0],[0]])

#find The first derivative with respect to x.
def y(y,yn):
    return 2*(y-yn)
#find The first derivative with respect to y.
def x(x,xn):
    return 2*(x-xn)
#find The first derivative with respect to t.
def t(e,tln,tsn,v):
     return -2*(e+tln-tsn)*v*v
#equation system f(x,y,e)
def d(x,xn,y,yn,e,tln,tsn,v):
    return (x-xn)*(x-xn)+(y-yn)*(y-yn)-(e+tln-tsn)*(e+tln-tsn)*v*v

#(x-xn)(x-xn)+(y-yn)(y-yn)-(e+tln-tsn)*(e+tln-tsn)*v**2


#finding the inverse of jacobian matrix
def Jacobian_inverse():
    #the elements of jacoban matrix 
    e11=x(xye[0,0],x1)
    e21=x(xye[0,0],x2)
    e31=x(xye[0,0],x3)
    e12=y(xye[1,0],y1)
    e22=y(xye[1,0],y2)
    e32=y(xye[1,0],y3)
    e13=t(xye[2,0],tl1,ts1,v)
    e23=t(xye[2,0],tl2,ts2,v)
    e33=t(xye[2,0],tl3,ts3,v)
    #puting the elements in the matrix
    A=[[e11,e12,e13], 
    [e21,e22,e23],
    [e31,e32,e33]]
    return np.array(np.linalg.inv(A))


#finding eqution results matrix
def main_eqution():
    f1=d(xye[0,0], x1, xye[1,0], y1, xye[2,0], tl1, ts1, v)
    f2=d(xye[0,0], x2, xye[1,0], y2, xye[2,0], tl2, ts2, v)
    f3=d(xye[0,0], x3, xye[1,0], y3, xye[2,0], tl3, ts3, v)
    return np.array([[f1],[f2],[f3]])


#Newton-Raphson method iterations
for i in range(4):
    #we repeat the following process to get higher accuracy
    B= Jacobian_inverse()   
    fun=main_eqution()
    #the matrix of answer=the matrix of answer-inverse of jacobian matrix*results matrix
    xye=np.subtract(xye,B.dot(fun))

#printing the solution
print('After 3 iterations ')       
print('x: ',round(xye[0,0],10))
print('y: ',round(xye[1,0],10))
print('time delta: ',round(xye[2,0],10))