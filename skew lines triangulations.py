from math import sin,cos
import numpy as np

def triangulate_1(x1,y1,z1,pitch1,yaw1,x2,y2,z2,pitch2,yaw2):

    #x1, y1, z1, l1, m1, n1 = find_lenpara(glat,glng,galt,lat1,lng1,alt1,yaw1,pitch1)
    #x2, y2, z2, l2, m2, n2 = find_lenpara(glat,glng,galt,lat2,lng2,alt2,yaw2,pitch2)

    #x1,y1,z1 = 30.801,-3.4658,3
    #x2,y2,z2 = 10.8971,30.580998,2.5

    l1,m1,n1 = find_dcs(yaw1,pitch1)
    l2,m2,n2 = find_dcs(yaw2,pitch2)

    h1 = (l1 * l1) + (m1 * m1) + (n1 * n1)
    h2 = (l2 * l2) + (m2 * m2) + (n2 * n2)
    #k10 = (pow(x1,2)*(pow(n1,2)+pow(m1,2)))+pow(y1,2)*(pow(l1,2)+pow(n1,2))+pow(z1,2)*(pow(m1,2)+pow(l1,2))-2*n1*m1*y1*z1-2*n1*l1*x1*z1-2*m1*l1*x1*y1)/h1+(pow(x2,2)*(pow(n2,2)+pow(m2,2))+pow(y2,2)*(pow(l2,2)+pow(n2,2))+pow(z2,2)*(pow(m2,2)+pow(l2,2))-2*n2*m2*y2*z2-2*n2*l2*x2*z2-2*m2*l2*x2*y2)/h2
    #k10 = (((pow(x1,2)*(pow(n1,2)+pow(m1,2)))+(pow(y1,2)*(pow(n1,2)+pow(l1,2)))+(pow(z1,2)*(pow(m1,2)+pow(l1,2)))-(2*n1*m1*y1*z1)-(2*n1*l1*x1*z1)-(2*m1*l1*x1*y1))/h1)+(((pow(x2,2)*(pow(n2,2)+pow(m2,2)))+(pow(y2,2)*(pow(n2,2)+pow(l2,2)))+(pow(z2,2)*(pow(m2,2)+pow(l2,2)))-(2*n2*m2*y2*z2)-(2*n2*l2*x2*z2)-(2*m2*l2*x2*y2))/h2)
    k4 = ((-2*x1*n1*n1 - 2*x1*m1*m1 + 2*n1*l1*z1 + 2*m1*l1*y1)/h1) + ((-2*x2*n2*n2 - 2*x2*m2*m2 + 2*n2*l2*z2 + 2*m2*l2*y2)/h2)
    k5 = ((-2*y1*l1*l1 - 2*y1*n1*n1 + 2*n1*m1*z1 + 2*m1*l1*x1)/h1) + ((-2*y2*l2*l2 - 2*y2*n2*n2 + 2*n2*m2*z2 + 2*m2*l2*x2)/h2)
    k6 = ((-2*z1*m1*m1 - 2*z1*l1*l1 + 2*n1*m1*y1 + 2*n1*l1*x1)/h1) + ((-2*z2*m2*m2 - 2*z2*l2*l2 + 2*n2*m2*y2 + 2*n2*l2*x2)/h2)

    k7 = (2*m1*l1/h1) + (2*m2*l2/h2)
    k8 = (2*n1*l1/h1) + (2*n2*l2/h2)
    k9 = (2*n1*m1/h1) + (2*n2*m2/h2)

    k1 = (((n1*n1) + (m1*m1))/h1) + (((n2*n2) + (m2*m2))/h2)
    k2 = (((l1*l1) + (n1*n1))/h1) + (((n2*n2) + (l2*l2))/h2)
    k3 = (((m1*m1) + (l1*l1))/h1) + (((m2*m2) + (l2*l2))/h2)

    A = np.array([[2*k1,-k7,-k8],[-k7,2*k2,-k9],[-k8,-k9,2*k3]])
    B = np.array([-k4,-k5,-k6])
    A_inv = np.linalg.inv(A)
    P = A_inv.dot(B)
    
    '''print([k1,k2,k3],[k4,k5,k6],[k7,k8,k9])
    print('A: ',A)
    print('B: ',B)
    print('A_inv:',A_inv)
    print('P: ',P)'''
    return P[0],P[1],P[2]

def triangulate_2(x1,y1,z1,pitch1,yaw1,x2,y2,z2,pitch2,yaw2):

    #x1,y1,z1 = 30.801,-3.4658,3
    #x2,y2,z2 = 10.8971,30.580998,2.5

    l1,m1,n1 = find_dcs(yaw1,pitch1)
    l2,m2,n2 = find_dcs(yaw2,pitch2)

    lamda1 = (-z1/n1)
    lamda2 = (-z2/n2)

    Px,Py = (x1+(lamda1*l1)),(y1+(lamda1*m1))
    Qx,Qy = (x2+(lamda2*l2)),(y2+(lamda2*m2))

    return ((Px+Qx)/2),((Py+Qy)/2),0

def find_dcs(yaw,pitch):
    #convert to radians
    yaw = (np.pi/180)*yaw
    pitch = (np.pi/180)*pitch
    
    sinp = sin(pitch)
    l = sinp*cos(yaw)
    m = sinp*sin(yaw)
    n = cos(pitch)
    
    return l,m,n

print(triangulate_2(30.801,-3.4658,3,100,61.9,10.8971,30.580998,2.5,95,17.53))
