# Steven Fuller
# June 9, 2020

#importing libraries to install dependent libraries.
import subprocess
import sys

# Defining a function that will be used to install dependent libraries
# Function will take 1 argument that is the dependent library to install
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Calling function to install each dependent library
install("matplotlib")
install("numpy")
install("python-math")
install("tk")


import matplotlib
import numpy as np
import math
from numpy import *
from tkinter import *
from tkinter.ttk import *
from tkinter.filedialog import *

root = Tk()
root.title("Neutron Transport")
root.geometry("410x150")

data = []
ellipsoid = []

#function used to open files
def open_data_file():
    global data 
    fileName = askopenfilename()
    if fileName is not None:
        data = np.load(fileName)

def open_ellipsoid_file():
    global ellipsoid 
    fileName = askopenfilename()
    if fileName is not None:
        ellipsoid = np.load(fileName)

def simulation(data, ellipsoid):
    #constants
    num_of_ellipse = 3
    num_of_neutrons = 100_000
    g = 9.81
    total_distance = 53.0

    #Slices
    slice_distance = 0.5 #want to get down to 0.01
    total_slices = int((1/2) * (1/slice_distance) * total_distance + 1)

    #counters
    reflect = 0 
    bounce = 0
    incident_neutrons = 0 
    detector_incident_reflected = 0 
    detector_incident_unreflected = 0
    unreflected_sensitivity = 0
    reflected_sensitivity = 0
    x_bounces = 0
    y_bounces = 0
    z_bounces = 0

    #ellipsoid deffinition
    X_0 = 0 
    Y_0 = 0
    Z_0 = 0.5 * total_distance
    R = 1

    #limits the length of the arrays 
    data = data[:num_of_neutrons]

    ellipsoid = ellipsoid[:num_of_ellipse]

    #using with open closes the file even if there is an error
    with open("reflected_file.txt", "w+") as f:
        f.write("j, c, a, reflect, bounce, (x,y,z) bounces, weight, weight, initial (x,y,z) position") 
        f.write("initial (x,y,z) velocity, initial wavelength, total initial velocity, final wavelength")
        f.write("total final velocity, final (x,y,z) velocity, fianl (x,y,z) position \n")

        print("start loops")
        for i in range(num_of_ellipse):

            last_update = 0
            for j,neutron in enumerate(data):   #total lines 7856990

                #Normalizes the flux to 1 MW operating power
                weight = neutron[0]/5
                #adjusts z position to 0 instead of -0.001
                neutron[3] = 0

                #repeated calculations  
                distance_over_time = (total_distance / neutron[6]) 
                slice_over_time = (slice_distance/neutron[6]) 
                total_initial_velocity = sqrt(neutron[4]**2 + neutron[5]**2 + neutron[6]**2)

                # final velocity calculations
                x_final_velocity = neutron[4]
                y_final_velocity = neutron[5] - g * distance_over_time
                z_final_velocity = neutron[6]
                total_final_velocity = sqrt(neutron[4]**2 + (neutron[5] - g * (total_distance / neutron[5]))**2 + neutron[6]**2)

                #final position calculations
                x_final_position = neutron[1] + neutron[4]*distance_over_time
                y_final_position = neutron[2] + neutron[5]*distance_over_time - (0.5*g)*distance_over_time**2
                z_final_position = neutron[3] + neutron[6]*distance_over_time

                #calculate the wavelength in angstroms
                wavelength_initial = 3956 / total_initial_velocity
                wavelength_final = 3956 / total_final_velocity

                incident_neutrons += weight

                if(z_final_position == total_distance and (x_final_position**2 + y_final_position**2) <= 0.0625):
                    #print("direct hit")
                    #neutrons that made it to detector without hitting anything
                    detector_incident_unreflected += weight  
                    #Calculate the n-nbar sensitivity in approximate ILL/year
                    unreflected_sensitivity += ( (weight * distance_over_time**2) / (1500000000.0*1.2))

                else:
                    for k in range(total_slices):
                        if(z_final_position < Z_0 and z_final_position > 0.0 and (x_final_position**2 + y_final_position**2 < (ellipsoid[i,2]**2)*(1-((z_final_position - Z_0)**2)/ellipsoid[i,0]**2))):
                            #final position calculations
                            x_final_position = neutron[1] + neutron[4] * slice_over_time
                            y_final_position = neutron[2] + neutron[5] * slice_over_time - (0.5*g)*(slice_over_time)**2
                            z_final_position = neutron[3] + slice_distance

                            #final velocity calculations
                            x_final_velocity = neutron[4]
                            y_final_velocity = neutron[5] - g * slice_over_time
                            z_final_velocity = neutron[6]
                            total_final_velocity = sqrt( x_final_velocity**2 + y_final_velocity**2 + z_final_velocity**2)

                            #for next step final calculations are saved as initial
                            x_initial_position = x_final_position
                            y_initial_position = y_final_position
                            y_initial_position = z_final_position
                            x_initial_velocity = x_final_velocity
                            y_initial_velocity = y_final_velocity
                            z_initial_velocity = z_final_velocity
                        
                        if(z_final_position < Z_0 and z_final_position > 0.0 and (x_final_position**2 + y_final_position**2 > (ellipsoid[i,2]**2)*(1-((z_final_position - Z_0)**2)/ellipsoid[i,0]**2))):
                            reflect = 1
                            bounce += 1

                            #bounce point
                            x_bounces = x_final_position - x_initial_velocity*(slice_distance/z_initial_velocity)
                            y_bounces = y_final_position - y_initial_velocity*(slice_distance/z_initial_velocity) + (0.5*g)*((slice_distance/z_initial_velocity)**2)
                            z_bounces = z_final_position - slice_distance

                            #save bounce point as initial point for next calculation
                            x_initial_position = x_bounces
                            y_initial_position = y_bounces
                            z_initial_position = z_bounces

                            #find the normilization of the uit normal vecto n-hat
                            norm = sqrt(pow(2*x_bounces/pow(ellipsoid[i,2],2),2) + pow(2*y_bounces/pow(ellipsoid[i,2],2),2) + pow(2*z_bounces/pow(ellipsoid[i,0],2),2))

                            #Calculate the components of the ellipsoid's unit normal vector n-hat
                            norm_x = (2*x_bounces/pow(ellipsoid[i,2],2))/norm
                            norm_y = (2*y_bounces/pow(ellipsoid[i,2],2))/norm
                            norm_z = (2*z_bounces/pow(ellipsoid[i,2],2))/norm

                            dot_product = (x_initial_velocity * norm_x) + (y_initial_velocity * norm_y) + (z_initial_velocity * norm_z)

                            #calculate the new reflected velocity
                            x_final_velocity = x_initial_velocity - 2 * dot_product * norm_x
                            y_final_velocity = y_initial_velocity - 2 * dot_product * norm_y
                            z_final_velocity = z_initial_velocity - 2 * dot_product * norm_z
                            total_final_velocity = sqrt( x_final_velocity**2 + y_final_velocity**2 + z_final_velocity**2)

                            #set initial velocity to final velocity
                            x_initial_velocity = x_final_velocity
                            y_initial_velocity = y_final_velocity
                            z_initial_velocity = z_final_velocity

                            #calculate new final position
                            x_final_position = x_initial_position + x_initial_velocity * (slice_distance/z_initial_velocity)
                            y_final_position = y_initial_position + y_initial_velocity * (slice_distance/z_initial_velocity) - (0.5*g) * (slice_distance/z_initial_velocity)**2
                            z_final_position = z_initial_position + z_initial_velocity * (slice_distance/z_initial_velocity)

                            #set initial position to final position
                            x_initial_position = x_final_position
                            y_initial_position = y_final_position
                            z_initial_position = z_final_position
                    
                            if(z_final_position[j] >= Z_0):
                                #final velocity calculations
                                x_final_velocity = x_initial_velocity
                                y_final_velocity = y_initial_velocity - g * (Z_0 / z_initial_velocity)
                                z_final_velocity = z_initial_velocity

                                #total velocity initial and final
                                total_initial_velocity = sqrt( x_initial_velocity**2 + y_initial_velocity**2 + z_initial_velocity**2)
                                total_final_velocity = sqrt( x_final_velocity**2 + y_final_velocity**2 + z_final_velocity**2)

                                #final position calculations
                                x_final_position = x_initial_position + x_initial_velocity * (Z_0 / z_initial_velocity)
                                y_final_position = y_initial_position + y_initial_velocity * (Z_0 / z_initial_velocity) - (0.5*g) * (Z_0 / z_initial_velocity)**2
                                z_final_position = z_initial_position + z_initial_velocity * (Z_0 / z_initial_velocity)

                                #wavelength calculations
                                initial_wavelength = 3956 / total_initial_velocity
                                final_wavelength = 3956 / total_final_velocity

                        #if neutrons have reached the end of the beam line
                        if(z_final_position >= total_distance -0.01 or z_final_position <= total_distance + 0.01):
                            if(x_final_position**2 + y_final_position**2 <= 0.0625):
                                detector_incident_reflected += weight
                                #calculate the n-nbar sensitivity
                                reflected_sensitivity += ((weight * ((total_distance - z_bounces)/z_final_velocity)**2) / (1500000000.0*1.2))
                            
                            else:
                                #print("I did not hit the detector \n")
                                pass 
                            
                #print out data to file
                output = (j,  ellipsoid[i][2], ellipsoid[i][0], reflect, bounce,
                    x_bounces, y_bounces, z_bounces,
                    weight, neutron[1], neutron[2], neutron[3], 
                    neutron[4], neutron[5], neutron[6],
                    wavelength_initial, total_initial_velocity,
                    wavelength_final, total_final_velocity,
                    x_final_velocity, y_final_velocity, z_final_velocity,
                    x_final_position, y_final_position, z_final_position)
                
                f.write(str(output))
                f.write('\n') 
                    
                reflect = 0
                bounce = 0
                x_bounces = 0
                y_bounces = 0
                z_bounces = 0

                value = (j+1) / num_of_neutrons * 100
                if value - last_update >= 10:
                    last_update = value
                    neutron_progress['value'] = value
                    root.update()

            print("The ellipsoid configuration is (a,b,c,f) = (", ellipsoid[i][0], ", ", ellipsoid[i][1], ", ", ellipsoid[i][2], ", ", ellipsoid[i][3], ")\n")
            print("The FINAL TOTAL value of the INCIDENT flux is: ", incident_neutrons, "n/s \n")
            print("The pre-value of the UNREFLECTED detector flux is: ", detector_incident_unreflected, "n/s \n")
            print("The FINAL TOTAL value of the UNREFLECTED+REFLECTED detector flux is: ", detector_incident_reflected + detector_incident_unreflected, "n/s \n")
            print("The UNREFLECTED n-nbar <Nt^2> sensitivity is :", unreflected_sensitivity, "\n")
            print("The UNREFLECTED+REFLECTED n-nbar <Nt^2> sensitivity is :", reflected_sensitivity + unreflected_sensitivity, "\n")

            incident_neutrons = 0
            detector_incident_reflected = 0
            detector_incident_unreflected = 0
            unreflected_sensitivity = 0
            reflected_sensitivity = 0
            
            value = (i+1) / num_of_ellipse * 100
            ellipsoid_progress['value'] = value
            root.update()   

text = Text(root,height=2, width=30)
text.insert(END, "Files must be of .npy format!")
text.grid(columnspan=3, row=0)

neutron_progress = Progressbar(root, orient=HORIZONTAL, length=250, mode='determinate')
neutron_progress.grid(column=1, row = 3, columnspan=2)
ellipsoid_progress = Progressbar(root, orient=HORIZONTAL, length=250, mode='determinate')
ellipsoid_progress.grid(column=1, row = 4, columnspan=2)

data_load = Button(root, text='Load data file', command=open_data_file, height=2, width=15)
data_load.grid(row=1, column=0)

ellipsoid_load = Button(root, text='Load ellipsoid file', command=open_ellipsoid_file, height=2, width=15)
ellipsoid_load.grid(row=1, column=1)

run = Button(root, text='Run', command=lambda : simulation(data, ellipsoid), height=2, width=7)
run.grid(row=1, column=2)

text2 = Text(root,height=2, width=20)
text2.insert(END, "Neutron Progress: ")
text2.grid(column=0, row=3)

text3 = Text(root,height=2, width=20)
text3.insert(END, "Ellipsoid Progress: ")
text3.grid(column=0, row=4)

root.mainloop()