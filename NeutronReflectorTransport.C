#include <cmath>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <iomanip>

using namespace std;

int main()
{
  const int numofellipse=3; //Number of ellipse configurations to check via MC
  const double pi=3.141592653; //Approximate value of pi
  const double g=9.81; //Acceleration due to gravity in [m]/[s^2]
  const double distance_total=53.0; //Total distance from source to final detector in [m]
  const double slice_distance = 0.5; //0.001; //Step distance to transport along z-axis
  const int total_slices = (1./2.)*(1./slice_distance)*distance_total+1; //Number of iterations for transport of a neutron in [mm] through the half-ellipsoid
  int reflect=0, bounce=0; //Value =0 if normal transport complete, while value =1 if reflection has occurred; and the bounce number accumulator
  double inc_int=0, det_int_reflected=0, det_int_unreflected=0; //Values of incident and detector neutron integrals/accumulators
  double unreflected_sensitivity=0, reflected_sensitivity; //Value of the n-nbar sensitivity as the summation of weight*t_final*t_final
  double int_array[numofellipse][numofellipse][numofellipse][numofellipse]; //Array to store c/a values and sensitivity values
  double get_weight, weight; //To be read in; will not change
  //To be read in from event file; to be used for transport and to check for reflections; spins will be ignored
  double get_x, x_initial, get_y, y_initial, get_z, z_initial, get_vx, vx_initial, get_vy, vy_initial, get_vz, vz_initial, get_time, t_initial, get_spin_x, get_spin_y, get_spin_z;
  double v_initial, v_final, l_initial, l_final; //Initial/final total velocities/wavelengths of the neutron
  //To be read out; will be the final phase space at the detector (all other events will be ignored and not saved)
  double x_final, y_final, z_final, vx_final, vy_final, vz_final, t_final;
  double X, X_0=0, get_a, a, Y, Y_0=0, get_b, b, Z, Z_0=0.5*distance_total, get_c, c, get_f, f, R=1; //For the definition of the ellipsoid
  double n_x, n_y, n_z, norm; //Define the normalized/unit normal vector components to a point on the ellipsoid
  double dotproduct; //Define the dot product of the initial velocity vector and the normalized normal vector
  double X_numbounce=0, Y_numbounce=0, Z_numbounce=0; //The coordinates of the reflection/bounce point
  int i, j, ei, ci, ai; //Iterators

  //Initialize in stream for uniformly and randomly thrown ellipsoid configurations passing the minimum radius geometric cut at z=0
  ifstream ellconfigcin("/Users/stevenfuller/Desktop/NewNeutronWork/WorkingEllipsoidConfigurations.dat");
  
  //Start looping over all ellipsoidal configurations given a value for a and b=c (surface of revolution along z-axis)
  for(ei=0; ei<numofellipse; ei++)
    {
      ellconfigcin >> get_a >> get_b >> get_c >> get_f; //Read through file line by line (configuration by configuration)
      a = get_a; b = get_b; c = get_c; f = get_f; //Set values of major and minor axes

      //Create a string to same the name and location of the output files based on configuration parameters
      string reflected_file = "/Users/stevenfuller/Desktop/NewNeutronWork/reflectedandtransportedESS_c_";
      char buffer [50]; int n;
      n = sprintf(buffer,"%f_a_%f.dat",c,a);
      reflected_file.append(buffer,n);

      //Initialize in/out streams
      ifstream cin("/Users/stevenfuller/Desktop/NewNeutronWork/new_anni_events.dat");
      ofstream reflcout(reflected_file);
      
      //Initialize loop over all neutron histories
      for(j=0; j<7856990/*number of lines in ANNI event file from Torsten*/; j++)
	{
	  //reflect=0;
	  if(j%1000000==0){cout << "The value of the iterator is j = " << j << endl;}
	  //Read in all data from input file
	  cin >> get_weight >> get_x >> get_y >> get_z >> get_vx >> get_vy >> get_vz >> get_time >> get_spin_x >> get_spin_y >> get_spin_z;
	  //Define weight of all 'individual' neutrons/tracks
	  weight = get_weight/5.0; //Normalize the flux to 1 MW of operating power
	  get_z = 0.0; //Start everything at z=0 rather than -1mm
	  //Save their velocities and then transport them to the end of the detector and (eventually) check their positions
	  vx_initial = get_vx; 
	  vy_initial = get_vy; 
	  vz_initial = get_vz;
	  vx_final = get_vx; 
	  vy_final = get_vy - g*(distance_total/get_vz); 
	  vz_final = get_vz;
	  v_initial = sqrt(get_vx*get_vx + get_vy*get_vy + get_vz*get_vz);
	  v_final = sqrt(get_vx*get_vx + (get_vy - g*(distance_total/get_vz))*(get_vy - g*(distance_total/get_vz)) + get_vz*get_vz);
	  x_initial = get_x; 
	  y_initial = get_y; 
	  z_initial = get_z; 
	  x_final = get_x + get_vx*((distance_total)/get_vz);
	  y_final = get_y + get_vy*((distance_total)/get_vz) - (0.5*g)*(((distance_total)/get_vz)*((distance_total)/get_vz));
	  z_final = get_z + get_vz*((distance_total)/get_vz);

	  //Calculate their wavelenth in angstroms
	  l_initial = 3956/v_initial;
	  l_final = 3956/v_final;

	  inc_int += weight; //Add up all incident neutrons coming from event file by weights

	  if((z_final == distance_total) && (x_final*x_final + y_final*y_final <= 0.0625)) //See if neutrons have already hit the detector
	    {
	      det_int_unreflected += weight; //Add up only those neutrons which hit detector of their own accord
	      unreflected_sensitivity += ((weight*(distance_total/get_vz)*(distance_total/get_vz))/(1500000000.0*1.2)); //Calculate the n-nbar sensitivity in approximate ILL/year
	      //Print out all pertinent transported data to .dat file
	      reflcout << setprecision(3) << fixed << showpoint << j << " " << c << " " << a << " " << reflect  << " " << bounce << " " << X_numbounce << " " << Y_numbounce << " " << Z_numbounce << " " << weight << " " << get_x << " " << get_y << " " << get_z << " " << x_final << " " << y_final << " " << z_final << " " << get_vx << " " << get_vy << " " << get_vz << " " << v_initial << " " << l_initial << " " << vx_final  << " " << vy_final  << " " << vz_final  << " " << v_final << " " << l_final << endl;
	    }
	  else //Otherwise, we must simulate to see what the effect a reflection will have on the neutron
	    {
	      //cout << "I'M INSIDE THE ELSE STATEMENT" << endl;
	      //Reset variables to start again
	      x_initial = get_x; 
		  y_initial = get_y; 
		  z_initial = get_z;
	      x_final = get_x; 
		  y_final = get_y; 
		  z_final = get_z;
	      vx_initial = get_vx; 
		  vy_initial = get_vy; 
		  vz_initial = get_vz;
	      vx_final = get_vx; 
		  vy_final = get_vy; 
		  vz_final = get_vz;

	      for(i=0; i<total_slices; i++)
		{
		  //Make sure that the neutron is inside (not on or outside) the ellipsoid
		  if((z_final<Z_0) && (z_final>=0.0) && (x_final*x_final + y_final*y_final < c*c*(1-((z_final-Z_0)*(z_final-Z_0))/(a*a))))
		    {
		      //cout << "LOOK AT ME, I'M INSIDE THE ELLIPSOID and on STEP #" << i << endl;
		      //Save their previous velocities before every iteration and then transport them to the surface of the ellipsoid, checking their positions at every step
		      vx_initial = vx_initial; vy_initial = vy_initial; vz_initial = vz_initial;
		      v_initial = sqrt(vx_initial*vx_initial+vy_initial*vy_initial+vz_initial*vz_initial);
		      //Update their previous positions
		      x_initial = x_initial; y_initial = y_initial; z_initial = z_initial;
		      //Update their final positions after the small step
		      x_final = x_initial + vx_initial*((slice_distance)/vz_initial);
		      y_final = y_initial + vy_initial*((slice_distance)/vz_initial) - (0.5*g)*(((slice_distance)/vz_initial)*((slice_distance)/vz_initial));
		      z_final = z_initial + slice_distance;
		      //Update their final velocities
		      vx_final = vx_initial; 
			  vy_final = vy_initial - g*(slice_distance/vz_initial); 
			  vz_final = vz_initial;
		      //cout << "MY Z-COMPONENT OF VELOCITY IS " << vz_final << endl;
		      v_final = sqrt(vx_final*vx_final + vy_final*vy_final + vz_final*vz_final);
		      //Now, for the next iteration, save their NEW STARTING positions and velocites
		      x_initial = x_final; 
			  y_initial = y_final; 
			  z_initial = z_final;
		      vx_initial = vx_final; vy_initial = vy_final; vz_initial = vz_final;
		    }

		  //If the neutron hits or moves just outside the ellipsoid...
		  if((z_final<Z_0) && (z_final>=0.0) && (x_final*x_final + y_final*y_final > c*c*(1-((z_final-Z_0)*(z_final-Z_0))/(a*a))))
		    {
		      //cout << "LOOK AT ME, I'M OUTSIDE THE ELLIPSOID BUT BEFORE THE CENTER and on STEP #" << i << endl;
		      //cout << "The ORIGINAL bounce point is (" << X_numbounce << "," << Y_numbounce << "," << Z_numbounce << ")" << endl;
		      reflect = 1;
		      bounce += 1;
		      //...go back to the initial value, which we will call the bounce/reflection point, which will then be used to calculate a gradient
		      X_numbounce = x_final - vx_initial*((slice_distance)/vz_initial);
		      Y_numbounce = y_final - vy_initial*((slice_distance)/vz_initial) + (0.5*g)*(((slice_distance)/vz_initial)*((slice_distance)/vz_initial));
		      Z_numbounce = z_final - slice_distance;
		      //cout << "The NEW bounce point is (" << X_numbounce << "," << Y_numbounce << "," << Z_numbounce << ")" << endl;
		      //Save this bounce point as the original starting point for the next iteration
		      x_initial = X_numbounce; 
			  y_initial = Y_numbounce; 
			  z_initial = Z_numbounce;
		      //Now that we have gone back to the point previous to going outside the ellipsoid, we must change the velocity components according to the
		      //law of reflection by computing the component-wise normalized normal vector to the surface: vin_f = vin_o - 2*(sum(vin_o*n_i,{i,1,3}))*n_i
		      //Find the norm(alization) of the unit normal vector n_hat
		      norm = sqrt(pow(2*X_numbounce/pow(c,2),2) + pow(2*Y_numbounce/pow(c,2),2) + pow(2*Z_numbounce/pow(a,2),2));
		      //Calculate the components of the ellipsoid's unit normal vector n_hat
		      n_x = (2*X_numbounce/pow(c,2))/norm; n_y = (2*Y_numbounce/pow(c,2))/norm; n_z = (2*Z_numbounce/pow(a,2))/norm;
		      dotproduct = (vx_initial * n_x) + (vy_initial * n_y) + (vz_initial * n_z); //Calculates the full dot product of vn_vec.n_hat
		      
			  //Calculate the new velocity by components for the reflected neutron at the end of this step
		      vx_final = vx_initial - 2*dotproduct*n_x; 
			  vy_final = vy_initial - 2*dotproduct*n_y; 
			  vz_final = vz_initial - 2*dotproduct*n_z;
		      //cout << "MY Z-COMPONENT OF VELOCITY AFTER THE BOUNCE IS " << vz_final << endl;
		      v_final = sqrt(vx_final*vx_final+vy_final*vy_final+vz_final*vz_final);

		      vx_initial = vx_final; 
			  vy_initial = vy_final; 
			  vz_initial = vz_final;

		      //Calculate the new final position during this step
		      x_final = x_initial + vx_initial*((slice_distance)/vz_initial);
		      y_final = y_initial + vy_initial*((slice_distance)/vz_initial) - (0.5*g)*(((slice_distance)/vz_initial)*((slice_distance)/vz_initial));
		      z_final = z_initial + vz_initial*((slice_distance)/vz_initial);
		      //cout << "MY NEW POSITION IS (" << x_final << "," << y_final << "," << z_final << ")" << endl;
		      //Now, for the next iteration, save their NEW STARTING positions
		      x_initial = x_final; 
			  y_initial = y_final; 
			  z_initial = z_final;
		    }
		}
	      //Now do a full transport to the detector once the neutron reaches beyond the half-way point (the half-ellipsoid's edge)
	      if(z_final>=Z_0)//This should only happen on the very last iteration of the loop
		{
		  //cout << "LOOK AT ME, I'M AFTER THE CENTER AND READY TO BE TRANSPORTED" << endl;
		  //Save their velocities and then transport them to the end of the rest of the distance to the detector
		  vx_initial = vx_initial; 
		  vy_initial = vy_initial; 
		  vz_initial = vz_initial;

		  vx_final = vx_initial; 
		  vy_final = vy_initial - g*(Z_0/vz_initial); 
		  vz_final = vz_final;

		  v_initial = sqrt(pow(vx_initial,2) + pow(vy_initial,2) + pow(vz_initial,2));
		  v_final = sqrt(pow(vx_final,2) + pow(vy_final,2) + pow(vz_final,2));

		  x_initial = x_initial; 
		  y_initial = y_initial; 
		  z_initial = z_initial;

		  x_final = x_initial + vx_initial*((Z_0)/vz_initial);
		  y_final = y_initial + vy_initial*((Z_0)/vz_initial) - (0.5*g)*(((Z_0)/vz_initial)*((Z_0)/vz_initial));
		  z_final = z_initial + vz_initial*((Z_0)/vz_initial);
		  //cout << "MY FINAL POSITION IS (" << x_final << "," << y_final << "," << z_final << ")" << endl;
		  //Calculate their wavelenth in angstroms
		  l_initial = 3956/v_initial;
		  l_final = 3956/v_final;
		  //cout << "This is another test of reflection: X_numbounce = " << X_numbounce << endl;
		}
	      if(z_final == distance_total) //See if neutrons have reached the end of the beamline
		{
		  //cout << "LOOK AT ME, I'VE MADE IT TO THE END OF THE BEAMLINE" << endl;
		  //Add up all neutrons hitting detector region from event file by weights
		  if(x_final*x_final + y_final*y_final <= 0.0625)
		    {
		      det_int_reflected += weight;
		      //cout << "AND I HIT THE DETECTOR!" << endl;
		      reflected_sensitivity += ((weight*((distance_total-Z_numbounce)/vz_final)*((distance_total-Z_numbounce)/vz_final))/(1500000000.0*1.2)); //Calculate the n-nbar sensitivity
		    }
		  if(x_final*x_final + y_final*y_final > 0.0625)
		    {
		      //cout << "AND I DIDN'T HIT THE DETECTOR!" << endl;
		    }
		  //Print out all pertinent transported data to .dat file
		  reflcout << setprecision(3) << fixed << showpoint << j << " " << c << " " << a << " " << reflect  << 
		  " " << bounce << " " << X_numbounce << " " << Y_numbounce << " " << Z_numbounce << " " << weight <<
		   " " << get_x << " " << get_y << " " << get_z << " " << x_final << " " << y_final << " " << z_final <<
		    " " << get_vx << " " << get_vy << " " << get_vz << " " << v_initial << " " << l_initial << " " << 
			vx_final  << " " << vy_final  << " " << vz_final  << " " << v_final << " " << l_final << endl;
		}
	      //Reset all reflection/bounce values
	      reflect = 0; 
		  bounce = 0; 
		  X_numbounce = 0; 
		  Y_numbounce = 0; 
		  Z_numbounce = 0;
	    }
	}
      cout << "The ellipsoid configuration is (a,b,c,f) = (" << get_a << "," << get_b << "," << get_c << "," << get_f << ")" << endl;
      cout << "The FINAL TOTAL value of the INCIDENT flux is: " << inc_int << " n/s" << endl;
      cout << "The pre-value of the UNREFLECTED detector flux is: " << det_int_unreflected << " n/s" << endl;
      cout << "The FINAL TOTAL value of the UNREFLECTED+REFLECTED detector flux is: " << det_int_reflected + det_int_unreflected << " n/s" << endl;
      cout << "The UNREFLECTED n-nbar <Nt^2> sensitivity is :" << unreflected_sensitivity << endl;
      cout << "The UNREFLECTED+REFLECTED n-nbar <Nt^2> sensitivity is :" << reflected_sensitivity + unreflected_sensitivity<< endl;

      inc_int = 0; det_int_reflected = 0; det_int_unreflected = 0; unreflected_sensitivity = 0; reflected_sensitivity = 0;
    }
}
