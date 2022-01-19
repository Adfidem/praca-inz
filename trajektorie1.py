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


G = 6.67430151515151515151515E-11 #Gravitational consatant

def period(mu,semi_major_axis_lengtht):
    #print((2*math.pi/math.sqrt(mu))*(semi_major_axis_lengtht**(3/2))/60/60/24,"days")
    return 2*math.pi/math.sqrt(mu)*semi_major_axis_lengtht**(3/2)#seconds!
def days_to_seconds(time):
    return time*24*60*60
def seconds_to_days(time):
    return time/60/60/24
def J2000_to_Julian_Day(time):
    return (time-2000)*365.25+2451545
def convert_Julian_Day_to_MJD(time):
    return time-2400000.5

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
        return self.mass*G/1000000000 #km^3⋅s^–2 
    def orbit(self, apoapsis, periapsis, eccentricity, inclination, orbited_body, mean_anomaly, epoch_MJD):
        self.periapsis  = periapsis #km
        self.apoapsis = apoapsis #km
        self.eccentricity = eccentricity #[-]
        self.inclination = math.pi/2-math.radians(inclination)#input in degrees
        self.orbited_body = orbited_body#class param of the body being orbited
        self.mean_anomaly = math.radians(mean_anomaly)#input in degrees
        self.epoch_MJD = epoch_MJD#input modified Julian Date
          
    def semi_major_axis_lengtht(self):#a
        #print((self.periapsis+self.apoapsis)/2)
        return (self.periapsis+self.apoapsis)/2
    def semi_minor_axis_lengtht(self):#b or c
        return (self.periapsis+self.apoapsis)/2*math.sqrt(1-self.eccentricity*self.eccentricity)

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
    
    def time(self,second, minute, hour, day, month, year): 
        if second > 59 or second < 0:
            sys.exit("seconds input out of range")
        elif minute > 59 or minute < 0:
            sys.exit("minute input out of range")
        elif hour > 24 or hour < 0:
            sys.exit("hour input out of range")
        elif day > 31 or day < 1:
            sys.exit("day input out of range")
        elif month > 12 or month < 1:
            sys.exit("month input out of range")
        elif year > 2099 or year < 1901:
            sys.exit("year input out of range")
        else:
            self.second = second
            self.minute = minute
            self.hour = hour
            self.day = day
            self.month = month
            self.year = year
    
    def convert_standard_time_notation_to_Julian_Day(self):#Boulet 1991   
        Jo = 367*self.year-int(7/4*(self.year+int((self.month+9)/12)))+int(275*self.month/9)+self.day+1721013.5
        UT = self.hour+self.minute/60+self.second/3600
        #print("jo", Jo)
        #print("UT", UT)
        return Jo + UT/24 #days
    def convert_standard_time_notation_to_MJD(self):#https://scienceworld.wolfram.com/astronomy/ModifiedJulianDate.html
        Jo = 367*self.year-int(7/4*(self.year+int((self.month+9)/12)))+int(275*self.month/9)+self.day+1721013.5-2400000.5
        UT = self.hour+self.minute/60+self.second/3600
        #print("jo", Jo)
        #print("UT", UT)
        return Jo + UT/24 #days
#    def velocity_unit_vector(self, X, Y, Z):
#        A = X/(self.semi_major_axis_lengtht**2)
#        B = Y
#        return 0
#    def set_current_MJD(self, current_MJD):
#        self.current_MJD = current_MJD
    
    
    
    def elipse_function(self,X,A,B):
        return B*math.sqrt(1-((X-A)**2)/(A**2))
    def area_under_elipse(self,x):
        A = self.semi_major_axis_lengtht()
        B = self.semi_minor_axis_lengtht()
        return quad(self.elipse_function,0, x, args=(A,B))[0]
    
    def azimuth_angle(self):#returns angle in rad
        time_current = days_to_seconds(self.convert_standard_time_notation_to_MJD())#days_to_seconds(self.Julian_Day())
        half_elipse_area = math.pi*self.semi_major_axis_lengtht()*self.semi_minor_axis_lengtht()/2
        P = period(self.orbited_body.mu(), self.semi_major_axis_lengtht())
        time0 = P/2#time of half a orbit
        unit_time = time_current % P #definiuje pozycje na orbicie (narazie zakladajac ze poczotek czasu to kiedy obiekt byl na peryapsie)
        #print(seconds_to_days(unit_time), "unit time")
        if unit_time == P/2:#check extreems (apoapsis)
            return math.pi
        if unit_time == 0:#priapsis
            return 0
        
        wich_half_of_orbit = 1
        if unit_time > P/2:
            wich_half_of_orbit = -1
            unit_time = unit_time - P/2
            
        swept_area = half_elipse_area*unit_time/time0
        X_max = self.semi_major_axis_lengtht()*2
        #print(X_max, "X_max")
        rotation_point = self.semi_major_axis_lengtht()-math.sqrt(self.semi_major_axis_lengtht()**2-self.semi_minor_axis_lengtht()**2)
        #print(rotation_point, "rotation")
        #print(half_elipse_area, "half_elipse_area")
        
        error = np.zeros(3) 
        error[0] = self.semi_major_axis_lengtht()#guessed x
        error[1] = 2#scale of guess
        error[2] = 0 #save it right (1) or left (2) of rotation point
        accuracy = 1 #km along x axis
        x_previous = 0
        #print(error[0], "initial x")
        
        
        
        while 1==1:
            area = self.area_under_elipse(error[0])
            
            #print(area,"area________________new loop pass, with entering x at", error[0])
            if error[1] > 1e+300:
                #max 1e+300
                #print("exited loop through error[1] bar exeption, with error at %f" % error[0])
                #print("and error[1] at", error[1])
                break
            
            triangle_area = np.zeros(2)
            if error[0] > rotation_point:
                #print("entered triangle 1")#subtract triangle, right of rotation point
                triangle_area[0] = 0.5*(error[0]-rotation_point)*(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht()))
                triangle_area[1] = -1
                error[2] = 1
                #print(triangle_area[0], "triangle_area", error[0]-rotation_point, "triangle base")
            else:
                #print("entered triangle 2")#add triangle, left of rotation point
                triangle_area[0] = 0.5*(rotation_point-error[0])*(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht()))
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
            print(error[0]/(self.semi_major_axis_lengtht()*2)*100 , "% x - of boundry line______________")'''
        
        
        if error[2] == 1:
            if wich_half_of_orbit == 1:#first half
                print("----------- exit angle 1:right of rotation_point, first half of orbit")
                angle = math.pi-math.atan(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht())/(error[0]-rotation_point))
                return angle
            elif wich_half_of_orbit == -1:#second half
                print("----------- exit angle 2:right of rotation_point, second half of orbit")
                angle = 2*math.pi-math.atan(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht())/(error[0]-rotation_point))
                return angle
            else:
                sys.exit("azimuth_angle failed to determin whether angle is + ro - on the y axis")
        elif error[2] == 2:
            if wich_half_of_orbit == 1:
                print("----------- exit angle 3:left of rotation_point, first half of orbit")
                angle = math.atan(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht())/(rotation_point-error[0]))
                return angle
            elif wich_half_of_orbit == -1:
                print("----------- exit angle 4:left of rotation_point, second half of orbit")
                angle = math.pi+math.atan(self.elipse_function(error[0],self.semi_major_axis_lengtht(),self.semi_minor_axis_lengtht())/(rotation_point-error[0]))
                return angle
            else:
                sys.exit("azimuth_angle failed to determine whether angle is + ro - on the y axis")
        else:
            sys.exit("azimuth_angle failed to determine whether angle is left or right of rotation_point")    
    

    
     
    def frame_of_referance_correction(self):#move x axis from center of elipse to foci
        #print(math.sqrt(self.semi_major_axis_lengtht()*self.semi_major_axis_lengtht()-self.semi_minor_axis_lengtht()*self.semi_minor_axis_lengtht()),"frame x correction")
        return math.sqrt(self.semi_major_axis_lengtht()*self.semi_major_axis_lengtht()-self.semi_minor_axis_lengtht()*self.semi_minor_axis_lengtht())
    def position_vector(self):
        calculated_azimuth_angle = self.azimuth_angle()
        return vector(self.semi_major_axis_lengtht()*self.sin(self.inclination)*self.cos(calculated_azimuth_angle)-self.frame_of_referance_correction(), self.semi_minor_axis_lengtht()*self.sin(self.inclination)*self.sin(calculated_azimuth_angle), self.semi_minor_axis_lengtht()*self.cos(self.inclination))
    def distance_from_orbited_body(self):
        r = self.position_vector()
        return math.sqrt(r*r)
    def speed(self):
        return math.sqrt(self.orbited_body.mu()*(2/self.distance_from_orbited_body()-1/self.semi_major_axis_lengtht()))
    def velocity_vector(self):
        return odeint()
        #return vector(0, math.sqrt(2*Sun.mu()*self.apoapsis/(self.periapsis*(self.apoapsis+self.periapsis))) ,0)#correct mu (and everything else xd)


Sun = param("Sun")
Sun.body(1.9885E+30)


Earth = param("Earth")
Earth.body(5.97237E+24)
#  apoapsis[km], periapsis[km], eccentricity[-], inclination[deg], orbited_body, mean_anomaly[deg], epoch_MJD [MJD]
Earth.orbit(152100000, 147095000, 0.0167086, 0, Sun, 358.617, convert_Julian_Day_to_MJD(J2000_to_Julian_Day(0)))#365.256
Earth.time(30,45,14,12,5,2004)
print(math.degrees(Earth.azimuth_angle()), "deg")
print(Earth.position_vector())
print(Earth.speed(),"speed")

asteroid_1996FG3 = param("1996FG3")
asteroid_1996FG3.orbit(212728172.12260202, 102474541.42333502, 0.35, 2, Sun, 202.32, 59600.0)
asteroid_1996FG3.time(30,45,14,12,5,2004)
print(math.degrees(asteroid_1996FG3.azimuth_angle()), "deg")
print(asteroid_1996FG3.position_vector())



'''
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