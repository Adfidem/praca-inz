# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 13:55:29 2021

@author: Szogunator
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import optimize as opt
import math
import sys
from scipy.integrate import quad
from scipy.integrate import odeint
from astropy.time import Time
from astropy.constants import G  #Gravitational consatant
from astropy.constants import M_earth  #Gravitational consatant
from astropy.constants import M_sun  #Gravitational consatant
from matplotlib.patches import Ellipse

def period(mu,semi_major_axis_length):
    #print((2*math.pi/math.sqrt(mu))*(semi_major_axis_length**(3/2))/60/60/24,"days")
    return 2*math.pi/math.sqrt(mu)*semi_major_axis_length**(3/2)#seconds!
def days_to_seconds(time):
    return time*24*60*60
def seconds_to_days(time):
    return time/60/60/24

class vector:
    #https://tomaszgolan.github.io/js-python/wyklady/js-python_w10/#wektor-iloczyn-skalarny
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
    def __str__(self):
        return "[{}, {}, {}]".format(self.x, self.y, self.z)
    def __mul__(self, w): # operator mnożenia, iloczyn skalarny
        return self.x*w.x + self.y*w.y + self.z*w.z
    def vector_mul(self, w):#iloczyn wektorowy
        return vector(self.y*w.z-self.z*w.y, self.z*w.x-self.x*w.z, self.x*w.y-self.y*w.x)
    def vector_add(self, w):
        return vector(self.x+w.x, self.y+w.y, self.z+w.z)
    def vector_subtract(self, w):
        #print(type(w.x),type(self.x))
        return vector(self.x-w.x, self.y-w.y, self.z-w.z)
    def vector_mul_scalar(self, s):
        return vector(self.x*s, self.y*s, self.z*s)
    def vector_div_scalar(self, s):
        return vector(self.x/s, self.y/s, self.z/s)
    def angle_between_vectors(self, w):#https://www.omnicalculator.com/math/angle-between-two-vectors#angle-between-two-3d-vectors-example
        return math.acos((self.x*w.x+self.y*w.y+self.z*w.z)/(math.sqrt(self.x**2+self.y**2+self.z**2)*math.sqrt(w.x**2+w.y**2+w.z**2)))
    def vector_length(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    def to_np_vector(self):
        return np.array([[self.x, self.y, self.z]]).T
#jak uzywac iloczynu wektororowego
#A = vector(1,2,3)
#B = vector(4,5,6)
#print(A.vector_mul(B))
#jak uzywac skalarnego
#print(A*B)


class param:
    k = vector(0,0,1)
    def __init__(self, name):
        self.name = name
    def body(self, mass):
        self.mass = mass #kg
    def mu(self):
        return self.mass*G.value/1000000000 #km^3⋅s^–2 
    def orbit(self, apoapsis, periapsis, eccentricity, inclination, arg_of_periapsis, orbited_body, mean_anomaly):#, epoch_MJD
        self.periapsis  = periapsis #km
        self.apoapsis = apoapsis #km
        self.eccentricity = eccentricity #[-]
        self.inclination = math.pi/2-math.radians(inclination)#input in degrees
        self.orbited_body = orbited_body#class param of the body being orbited
        self.mean_anomaly = math.radians(mean_anomaly)#input in degrees
        #self.epoch_MJD = epoch_MJD#input modified Julian Date - I don't think this is neaded any more
        #self.ascending_node = ascending_node#input in degrees
        self.arg_of_periapsis = math.radians(arg_of_periapsis)#input in degrees
          
    def semi_major_axis_length(self):#a
        return (self.periapsis+self.apoapsis)/2
    def semi_major_axis_unit_vector(self):
        answer = vector(self.semi_major_axis_length()*self.sin(self.inclination)*self.cos(self.arg_of_periapsis), self.semi_minor_axis_length()*self.sin(self.inclination)*self.sin(self.arg_of_periapsis), self.semi_minor_axis_length()*self.cos(self.inclination))
        return answer.vector_div_scalar(self.semi_major_axis_length())
    def semi_minor_axis_length(self):#b or c
        return (self.periapsis+self.apoapsis)/2*math.sqrt(1-self.eccentricity*self.eccentricity)
    def semi_minor_axis_unit_vector(self):
        answer = vector(self.semi_major_axis_length()*self.sin(self.inclination)*self.cos(self.arg_of_periapsis+math.pi/2), self.semi_minor_axis_length()*self.sin(self.inclination)*self.sin(self.arg_of_periapsis+math.pi/2), self.semi_minor_axis_length()*self.cos(self.inclination))
        return  answer.vector_div_scalar(self.semi_minor_axis_length())
    #cos i sin zapewniaja sprawdzaja czy kat nie jest zbyt blisko pi/2 i pi
    def cos(self, angle_rad):
        if abs(angle_rad) > math.pi/2-0.00000001 and abs(angle_rad) < math.pi/2+0.00000001:
            return 0
        else:
            return math.cos(angle_rad)
    def sin(self, angle_rad):
        if abs(angle_rad) > math.pi-0.00000001 and abs(angle_rad) < math.pi+0.00000001:
            return 0
        else:
            return math.sin(angle_rad)
    
    def set_time(self, isot_time):
        self.MJD = isot_time.mjd
    
    def ellipse_function(self,X,A,B):
        return B*math.sqrt(1-((X-A)**2)/(A**2))
    def area_under_ellipse(self,x):
        A = self.semi_major_axis_length()
        B = self.semi_minor_axis_length()
        return quad(self.ellipse_function,0, x, args=(A,B))[0]
    
    def azimuth_angle(self):#returns angle in rad + arg_of_periapsis (so true anomaly + arg_of_periapsis)
        time_current = days_to_seconds(self.MJD)#days_to_seconds(self.Julian_Day())
        half_ellipse_area = math.pi*self.semi_major_axis_length()*self.semi_minor_axis_length()/2
        P = period(self.orbited_body.mu(), self.semi_major_axis_length())
        time0 = P/2#time of half a orbit
        unit_time = time_current % P #definiuje pozycje na orbicie (narazie zakladajac ze poczotek czasu to kiedy obiekt byl na peryapsie)
        #print(seconds_to_days(unit_time), "unit time")
        if unit_time == P/2:#check extreems (apoapsis)
            return math.pi
        if unit_time == 0:#priapsis
            return 0
        
        which_half_of_orbit = 1
        if unit_time > P/2:
            which_half_of_orbit = -1
            unit_time = unit_time - P/2
            
        swept_area = half_ellipse_area*unit_time/time0
        X_max = self.semi_major_axis_length()*2
        #print(X_max, "X_max")
        rotation_point = self.semi_major_axis_length()-math.sqrt(self.semi_major_axis_length()**2-self.semi_minor_axis_length()**2)
        #print(rotation_point, "rotation")
        #print(half_ellipse_area, "half_ellipse_area")
        
        error = np.zeros(3) 
        error[0] = self.semi_major_axis_length()#guessed x
        error[1] = 2#scale of guess
        error[2] = 0 #save it right (1) or left (2) of rotation point
        accuracy = 1 #km along x axis
        x_previous = 0
        #print(error[0], "initial x")
        
        
        
        while 1==1:
            area = self.area_under_ellipse(error[0])
            
            #print(area,"area________________new loop pass, with entering x at", error[0])
            if error[1] > 1e+300:
                #max 1e+300
                #print("exited loop through error[1] bar exeption, with error at %f" % error[0])
                #print("and error[1] at", error[1])
                break
            
            triangle_area = np.zeros(2)
            if error[0] > rotation_point:
                #print("entered triangle 1")#subtract triangle, right of rotation point
                triangle_area[0] = 0.5*(error[0]-rotation_point)*(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length()))
                triangle_area[1] = -1
                error[2] = 1
                #print(triangle_area[0], "triangle_area", error[0]-rotation_point, "triangle base")
            else:
                #print("entered triangle 2")#add triangle, left of rotation point
                triangle_area[0] = 0.5*(rotation_point-error[0])*(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length()))
                triangle_area[1] = 1
                error[2] = 2
                #print(triangle_area[0], "triangle_area", error[0]-rotation_point, "triangle base")
            
            if (area+(triangle_area[0]*triangle_area[1]))-swept_area >= 0:
                '''print("______________________")
                print((area+(triangle_area[0]*triangle_area[1])), "area current")
                print("%f" % swept_area, "swept_area")
                print((area+(triangle_area[0]*triangle_area[1]))-swept_area , "differance")
                print("entered option 1 : left of x")'''
                error[1] = error[1]*2
                error[0] = error[0] - X_max/error[1]
                #print(error[0]/error[1], "error[0]/error[1]")
            elif (area+(triangle_area[0]*triangle_area[1]))-swept_area < 0:
                '''print("______________________")
                print((area+(triangle_area[0]*triangle_area[1])), "area current")
                print("%f" % swept_area, "swept_area")
                print((area+(triangle_area[0]*triangle_area[1]))-swept_area , "differance")
                print("entered option 2 right of x")'''
                error[1] = error[1]*2
                error[0] = error[0] + X_max/error[1]
            else:
                sys.exit("azimuth_angle loop failed to meet continue or exit criteria")
            
            if abs(x_previous - error[0]) < accuracy:#check chanege of x, exit if change less than accuracy
                #print("exited loop through x accuracy, with error at", error)
                break
            x_previous = error[0]
            '''print("continued loop with error[0] at %f" % error[0])
            print("and error[1] at", error[1])
            print(error[0]/(self.semi_major_axis_length()*2)*100 , "% x - of boundry line______________")'''
        
        
        if error[2] == 1:
            if which_half_of_orbit == 1:#first half
                #print("----------- exit angle 1:right of rotation_point, first half of orbit")
                angle = math.pi-math.atan(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length())/(error[0]-rotation_point))
                return angle
            elif which_half_of_orbit == -1:#second half
                #print("----------- exit angle 2:right of rotation_point, second half of orbit")
                angle = 2*math.pi-math.atan(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length())/(error[0]-rotation_point))
                return angle
            else:
                sys.exit("azimuth_angle failed to determin whether angle is + ro - on the y axis")
        elif error[2] == 2:
            if which_half_of_orbit == 1:
                #print("----------- exit angle 3:left of rotation_point, first half of orbit")
                angle = math.atan(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length())/(rotation_point-error[0]))
                return angle
            elif which_half_of_orbit == -1:
                #print("----------- exit angle 4:left of rotation_point, second half of orbit")
                angle = math.pi+math.atan(self.ellipse_function(error[0],self.semi_major_axis_length(),self.semi_minor_axis_length())/(rotation_point-error[0]))
                return angle 
            else:
                sys.exit("azimuth_angle failed to determine whether angle is + ro - on the y axis")
        else:
            sys.exit("azimuth_angle failed to determine whether angle is left or right of rotation_point")    
    

    
     
    def frame_of_reference_correction(self):#move x axis from center of ellipse to foci
        #print(math.sqrt(self.semi_major_axis_length()*self.semi_major_axis_length()-self.semi_minor_axis_length()*self.semi_minor_axis_length()),"frame x correction")
        return math.sqrt(self.semi_major_axis_length()*self.semi_major_axis_length()-self.semi_minor_axis_length()*self.semi_minor_axis_length())
    def position_vector(self):
        calculated_azimuth_angle_and_arg_of_periapsis = self.azimuth_angle() + self.arg_of_periapsis
        return vector(self.semi_major_axis_length()*self.sin(self.inclination)*self.cos(calculated_azimuth_angle_and_arg_of_periapsis)-self.frame_of_reference_correction(), self.semi_minor_axis_length()*self.sin(self.inclination)*self.sin(calculated_azimuth_angle_and_arg_of_periapsis), self.semi_minor_axis_length()*self.cos(self.inclination))
    def distance_from_orbited_body(self):
        r = self.position_vector()
        return math.sqrt(r*r)
    def speed(self):
        return math.sqrt(self.orbited_body.mu()*(2/self.distance_from_orbited_body()-1/self.semi_major_axis_length()))
        
    def velocity_vector(self):
        #https://math.stackexchange.com/questions/655853/ellipse-tangents-in-3d
        c = vector(-self.frame_of_reference_correction(),0,0)#center of ellipse, changed to "-" correction
        paramiter = self.azimuth_angle()
        u = self.semi_major_axis_unit_vector().vector_mul_scalar(self.semi_major_axis_length()*math.sin(paramiter))
        v = self.semi_minor_axis_unit_vector().vector_mul_scalar(self.semi_minor_axis_length()*math.cos(paramiter))
        tangent_unit_vector = c.vector_subtract(u.vector_add(v))
        return tangent_unit_vector.vector_mul_scalar(self.speed())
        #return vector(0, math.sqrt(2*Sun.mu()*self.apoapsis/(self.periapsis*(self.apoapsis+self.periapsis))) ,0)#correct mu (and everything else xd)
    def specific_angular_momentum(self):
        return self.position_vector().vector_mul(self.velocity_vector())#trzeba pomyslec !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #return self.velocity_vector().vector_mul(self.position_vector())
    def position_of_ascending_node(self):
        k = vector(0,0,1)
        return k.vector_mul(self.specific_angular_momentum())
    def angle_to_ascending_node(self):
        return self.position_vector().angle_between_vectors(self.position_of_ascending_node())

def draw_elipse(semi_minor_axis_length, semi_major_axis_length):
    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=(0,0), width = semi_minor_axis_length, height = semi_major_axis_length, edgecolor='r', fc='None', lw=2)
    ax.set_xlim(-semi_major_axis_length, semi_major_axis_length)
    ax.set_ylim(-semi_major_axis_length, semi_major_axis_length)
    ax.add_patch(ellipse)

def get_2D_elipse_parameters_from_3D_points(origin, destination):
    #unify frame of reference, center of ellipse at (0,0,0)
    origin_position = origin.position_vector().vector_add(vector(origin.frame_of_reference_correction(),0,0))
    destination_position = destination.position_vector().vector_add(vector(destination.frame_of_reference_correction(),0,0))
    
    #return to 2D ellipse / convert to local frame of reference
    #https://math.stackexchange.com/questions/1167717/transform-a-plane-to-the-xy-plane
    reference_plane_vector = vector(0,0,1)
    plane_vector = origin_position.vector_mul(destination_position)
    #plane_unit_vector = plane_vector.vector_div_scalar(plane_vector.vector_length())
    
    costau = math.cos(plane_vector.angle_between_vectors(reference_plane_vector))
    sintau = math.sin(plane_vector.angle_between_vectors(reference_plane_vector))
    rotation_vector = plane_vector.vector_mul(reference_plane_vector)
    unit_rotation_vector = rotation_vector.vector_div_scalar(rotation_vector.vector_length())
    
    rotation_matrix = np.array([
            [costau+unit_rotation_vector.x**2*(1-costau), unit_rotation_vector.x*unit_rotation_vector.y*(1-costau), unit_rotation_vector.y*sintau], 
            [unit_rotation_vector.x*unit_rotation_vector.y*(1-costau), costau+unit_rotation_vector.y**2*(1-costau),-unit_rotation_vector.x*sintau],
            [-unit_rotation_vector.y*sintau, unit_rotation_vector.x*sintau, costau]
    ])
    
    destination_position_2D = rotation_matrix @ destination_position.to_np_vector()
    origin_position_2D = rotation_matrix @ origin_position.to_np_vector()
    
    x1 = origin_position_2D[[0]]
    y1 = origin_position_2D[[1]]
    x2 = destination_position_2D[[0]]
    y2 = destination_position_2D[[1]]
    
    if(origin_position_2D[[2]] > 1 or destination_position_2D[[2]] > 1):
        string = ' '.join(["get_2D_elipse_parameters_from_3D_point failed to nulify 'z' coordinate", np.array2string(origin_position_2D), "derived origin coordinate \n" , np.array2string(destination_position_2D), "derived estinatio coordinate"])
        sys.exit(string)
    b2 = (x2**2*y1**2-x1**2*y2**2)/(x2**2-x1**2)
    a = np.sqrt(x1**2/(1-y1**2/b2))
    b = np.sqrt(b2)
    
    if(b > a):#make shure semi_major_axis_length > semi_minor_axis_length
        j = a
        a = b
        b = j
    
    semi_major_axis_length = a
    semi_minor_axis_length = b
    
    draw_elipse(semi_minor_axis_length, semi_major_axis_length)
    
    answer = np.zeros(2)
    answer[0] = semi_major_axis_length
    answer[1] = semi_minor_axis_length
    print(answer)
    return answer


def transfer_to_massles_body(time_of_departure, time_of_arrival, origin ,destination):#in MJD
    ship = param("ship")
    origin.set_time(time_of_departure)
    destination.set_time(time_of_arrival)
    
    ship.orbit(origin.apoapsis, origin.periapsis, origin.eccentricity, origin.inclination, origin.arg_of_periapsis, origin.orbited_body, origin.mean_anomaly)
    ship.set_time(time_of_departure)
    
    print(origin.frame_of_reference_correction(), destination.frame_of_reference_correction(), "f.o.r.c")
    
    print(destination.semi_major_axis_length())
    axis = get_2D_elipse_parameters_from_3D_points(origin, destination)
    
     
    
    
    '''
    #change k - make shure it't the right ascending node
    #point2 = destination.position_of_ascending_node() #ascending_node
    #point1 = origin.position_vector()
    point1 = vector(3.98,0.2,0)
    point2 = vector(-3.98,0.2,0)
    
    print(point1,point2)
    
    
    #not shure if it should be + or -
    #point1.x = point1.x + origin.frame_of_reference_correction()
    #point2.x = point2.x + destination.frame_of_reference_correction()
    
    print(point1.x, point2.x)
    #1 - x y z
    #2 - p q r
    print("---", point2.x**2*point1.y**2+point2.x**2*point1.z**2, "----", -point1.x**2*point2.y**2-point1.x**2*point2.z**2)
    b2 = (point2.x**2*point1.y**2+point2.x**2*point1.z**2-point1.x**2*point2.y**2-point1.x**2*point2.z**2)/(point2.x**2-point1.x**2)
    print(b2)
    b = math.sqrt(b2)
    a = math.sqrt((point2.x**2*b2)/(b2-point2.y**2-point2.z**2))
    
    print("a", a)
    print("b", b)
    #find_orbit(ship)
    #opt.minimize(find_orbit, )
        
    '''
    return

def find_orbit(ship):#get ellipse from two points
    
    
    
    #answer = ship.velocity_vector()
    return #answer*answer

def transfer_from_massles_body():
    
    return

def planetary_departure():
    
    return
def planetary_rendezvous():
    
    return



def solve():
    return 0


start_of_analysis = param("time start")
#second, minute, hour, day, month, year
start_of_analysis.set_time(Time('2023-01-01T00:00:00', format = 'isot'))
end_of_analysis = param("time end")
end_of_analysis.set_time(Time('2098-01-01T00:00:00', format = 'isot'))

Sun = param("Sun")
Sun.body(M_sun.value)

Earth = param("Earth")
Earth.body(M_earth.value)
epoch_earth = Time(2000.0, format = 'jyear')
#  apoapsis[km], periapsis[km], eccentricity[-], inclination[deg], arg_of_periapsis [deg], orbited_body, mean_anomaly[deg],#### epoch_MJD [MJD]
Earth.orbit(152100000, 147095000, 0.0167086, 0, 0, Sun, 358.617)#365.256 #, epoch_earth.mjd

Earth.set_time(Time('2004-05-12T14:45:30', format = 'isot'))
asteroid_1996FG3 = param("1996FG3")
asteroid_1996FG3.orbit(212728172.12260202, 102474541.42333502, 0.35, 2, 24.08, Sun, 202.32) #, 59600.0
asteroid_1996FG3.set_time(Time('2004-05-12T14:45:30', format = 'isot'))




transfer_to_massles_body(Time('2004-05-12T14:45:30', format = 'isot'),Time('2008-05-12T14:45:30', format = 'isot'), Earth, asteroid_1996FG3)




'''
print(math.degrees(Earth.azimuth_angle()), "deg")
print(Earth.position_vector())
print(Earth.speed(),"speed")
print(math.degrees(asteroid_1996FG3.azimuth_angle()), "deg")
print(asteroid_1996FG3.position_vector())
Earth.velocity_vector()
asteroid_1996FG3.velocity_vector()
print(Earth.position_of_ascending_node(), asteroid_1996FG3.position_of_ascending_node())



def simpleHohman(r1,r2,mu):
    deltaVd = math.sqrt(mu/r1)*(math.sqrt(2*r2/(r1+r2))-1)
    print(deltaVd*1000)
    deltaVa = math.sqrt(mu/r2)*(1-math.sqrt(2*r1/(r1+r2)))
    print(deltaVa*1000)
    return(abs(deltaVa+deltaVd))


print(simpleHohman(6700,10000,Earth.mu())*1000)
'''

"""
step = 10
target = np.zeros((step,20))

ship = np.zeros((step,20))

"""