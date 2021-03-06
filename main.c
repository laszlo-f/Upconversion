#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <time.h>		//random seed
#include <assert.h>

#include "gasdev.c"
#include "nrutil.c"
#include "splint.c"
#include "spline.c"
#include "ran1.c"
#include "prabsorbed.c"
int
main (int argc, char *argv[])
{

#include "declarations.c"
#include "read_input.c"
  if (min_emission<min_emitter_absorption || min_emission<min_absorption) {
	  printf("The minimum emission of %f is beyond the domain of the emissiom abosoption spectrum of %f and/or absorption spectrum of %f\n",min_emission, min_emitter_absorption, min_absorption);
  }
  if (max_emission>max_emitter_absorption || max_emission>max_absorption) {
	  printf("The maximum emission of %f is beyond the domain of the emissiom abosoption spectrum of %f and/or absorption spectrum of %f\n",max_emission, max_emitter_absorption, max_absorption);
  }
  esc = 0;
  total_energy = 0.0;
  N_phot = 1e4 * bin_N;		//check me for convergence - recommend at least 10^4 bin_N
  for (i = 1; i <= bin_N; i++)	//initialize array
    {
      depth_bin[i] = 0;
      reabs_bin[i] = 0;
      new_reabs_bin[i] = 0;
      singlet_depth_bin[i]=0;
      singlet_bin[i]=0;
      new_singlet_bin[i]=0;
    }

  for (i = 1; i <= N_phot; i++)	//number of incident photons
{
      splint (spec_x, spec_y, spec_y2, spec_N, ran1 (), &lambda);	//choose a wavelength
	  total_energy += 6.626e-34 * 2.9979e8 / lambda / 1e-9;
	  if (ran1()> pr_absorbed(lambda,min_wavelength,sctA))	//test if photon is absorbed by solar cell
	{
		if (lambda > max_emitter_absorption || lambda > max_absorption || lambda < min_emitter_absorption || lambda < min_absorption) {//test if photon is in suitable absorption range
		//if not assume absorbtion is 0 so photon is not absorbed by upconverter
		A1=0;
		emit_absorption=0;	
		}

	  else {
	  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	  splint (emit_abs_x, emit_abs_y, emit_abs_y2, emit_abs_N, lambda, &emit_absorption);	//get absorption coefficienct
	  A1 = fabs(A1); //eliminate negative values caused by noise
	  
	  if (emit_absorption < emit_abs_noise){ //eliminate values caused by noise
		  emit_absorption =0;
		  }
	  
	  A1 *= c;		//sensitizer concentration
	  emit_absorption *= emit_concentration;
	  }
	  
	  z = -1.0 * log (ran1 ()) / log (10.0) / (A1+emit_absorption);	//check to see if photon is absorbed
	  if (z < d && z > 0.0)	//incoming photon is absorbed
	    {
	      j = z / d * bin_N + 1;
	      if(ran1()>(emit_absorption/(A1+emit_absorption))){
		      //photon was absorbed by sensitizer
	        depth_bin[j] += 1;
			
	      } else {
		      //photon was absorbed by emitter, singlet is excited
		      singlet_depth_bin[j] += 1;
			  
	      }
	    }
	  else
	    {
	      //z = d;
	      if (LA == 1)
		dx = gasdev ();
	      else
		dx = 0.0;
	      if (LA == 1)
		dy = gasdev ();
	      else
		dy = 0.0;
	      
		  dz = gasdev ();
		  
		  
	      //normalize
	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dx = dx / dr;
	      dy = dy / dr;
	      dz = dz / dr;
	      dr = sqrt (dx * dx + dy * dy + dz * dz);

	      dst = -1.0 * log (ran1 ()) / log (10.0) / (A1+emit_absorption);
	      z = d - sqrt (dst * dz * dst * dz);	//increment z
		//distance above the mirror that the photon was absorbed

	      if (z > 0.0 && z < d && BR == 1)	//reabsorb off the back reflector
		{
		  j = z / d * bin_N + 1;
	      if(ran1()>(emit_absorption/(A1+emit_absorption))){
		      //photon was absorbed by sensitizer
		    depth_bin[j] += 1;
	      } else {
		      //photon was absorbed by emitter, singlet is excited
		      singlet_depth_bin[j] += 1;	  
	      }
		}
	    }			//if absorbed
	}
    } 				//i
  //one sun is 100 mW/cm^2
  run_time = total_energy / 100e-3 / C;

//NOW GET A SELF CONSISTENT SOLUTION
  for (iterate = 1; iterate <= reabsorptioncycles; iterate++)
    {
      esc = 0;

      for (i = 1; i <= bin_N; i++)
	{
	  reabs_bin[i] = new_reabs_bin[i];//copy data
	  new_reabs_bin[i] = 0.0;//reset

	  singlet_bin[i] = new_singlet_bin[i];//copy data
          new_singlet_bin[i] = 0.0;//reset
	}

      for (i = 1; i <= bin_N; i++)
	{
	  k_phiS = (depth_bin[i] + reabs_bin[i]) / run_time / (d / bin_N);	//add in reabsorption term
	  
	  NT = (emit_concentration/(emit_concentration+c*exp(-delta_E/(kB*temperature)))) * (-k1 + sqrt (k1 * k1 + 4 * k_phiS * k2)) / 2 / k2; 
	  //triplet number based on rate equation at steady state
	  //boltzmann factor of (1/(1+exp(-delta_E/(kB*temperature)))) is included
	  //because we assume triplets in the sensitizer cannot be used
	  //if you don't want this assumption just make delta_E big
	  
	  f2 = k2 * NT / (k2 * NT + k1);
	
	  //propagate photons produced by upconversion out
	  for (j = 1; j <= depth_bin[i] + reabs_bin[i]; j++)
	    {
	      splint (emi_x, emi_y, emi_y2, emi_N, ran1 (), &lambda);	//choose a wavelength

		  if (lambda > max_emitter_absorption || lambda > max_absorption || lambda < min_emitter_absorption || lambda < min_absorption) {//test if photon is in suitable absorption range
		  //if not assume absorbtion is 0 so photon is not absorbed by upconverter
		  A1=0;
		  emit_absorption=0;	
		  }

	      else {
		  
		  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	      splint (emit_abs_x, emit_abs_y, emit_abs_y2, emit_abs_N, lambda, &emit_absorption); //get absorption coefficient
	      A1 = fabs(A1); //eliminate negative values caused by noise
	      
		  if (emit_absorption < emit_abs_noise){ //eliminate values caused by noise
		  emit_absorption =0;
		  }
	      
		  A1 *= c;		//sensitizer concentration
	      emit_absorption *= emit_concentration;
		  }

	      dx = gasdev ();
	      dy = gasdev ();
	      dz = gasdev ();

	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dx = dx / dr;
	      dy = dy / dr;
	      dz = dz / dr;
	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dst = -1.0 * log (ran1 ()) / log (10.0) / (A1+emit_absorption); //beer-lambert law

	      z = d * (1.0 * i - 0.5) / (1.0 * bin_N) + dz * dst;//how far the emitted photon travels

	      if (z < 0.0)
		{		//if we escaped
		  if (ran1()< pr_absorbed(lambda,min_wavelength,sctA))//use tauc theory to test if photon is absorbed
		    {		//photon energy must be absorbed by solar cell
		      esc += eta_c * f2 * 0.5;		
		    }
		}
	      else if (z < d)
		{
		  h = z / d * bin_N + 1;//which bin did absorption occur in?  cast float to int
			if(ran1()>(emit_absorption/(A1+emit_absorption))){
		  		new_reabs_bin[h] += eta_c * f2 * 0.5;
			 } else {
				new_singlet_bin[h]+= eta_c * f2 * 0.5;
			}
		}

	      if (z > d && BR == 1)
		{
		  if (LA == 1)
		    {
		      dx = gasdev ();
		      dy = gasdev ();
		      dz = gasdev ();

		      dr = sqrt (dx * dx + dy * dy + dz * dz);
		      dx = dx / dr;
		      dy = dy / dr;
		      dz = dz / dr;
		      dr = sqrt (dx * dx + dy * dy + dz * dz);
			}

		  dst = -1.0 * log (ran1 ()) / log (10.0) /(A1+emit_absorption);

		  z = d - sqrt (dst * dz * dst * dz);	//increment z
		// distance the photon travels from the mirror.  It starts at the mirror because we know it wasn't absorbed on the back side of the mirror
		  if (z < 0.0)
		    {		//if we escaped
		      if (ran1()< pr_absorbed(lambda,min_wavelength,sctA))//use tauc theory to test if photon is absorbed
			{	//photon energy must be absorbed by solar cell
			  esc += eta_c * f2 * 0.5;
			}
		    }
		  else
		    {
		      h = z / d * bin_N + 1;
			if(ran1()>(emit_absorption/(A1+emit_absorption))){
                    		  //photon was absorbed by sensitizer
                      		new_reabs_bin[h] += eta_c * f2 * 0.5;
						
                  	} else {
                      		//photon was absorbed by emitter, singlet is excited
                  		new_singlet_bin[h]+= eta_c * f2 * 0.5;

			}
                    }           //reabsorption and escapage is fractional photons, more efficient that way for calcualtion
		}
	    }

//propagate photons from emitter absorption-emission out
	  for (j = 1; j <= singlet_depth_bin[i]+singlet_bin[i]; j++){
	      splint (emi_x, emi_y, emi_y2, emi_N, ran1 (), &lambda);	//choose a wavelength
	      
		  if (lambda > max_emitter_absorption || lambda > max_absorption || lambda < min_emitter_absorption || lambda < min_absorption) {//test if photon is in suitable absorption range
		  //if not assume absorbtion is 0 so photon is not absorbed by upconverter
		  A1=0;
		  emit_absorption=0;	
		  }

	      else {
		  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	      splint (emit_abs_x, emit_abs_y, emit_abs_y2, emit_abs_N, lambda, &emit_absorption); //get absorption coefficient
	      A1 = fabs(A1); //eliminate negative values caused by noise
	      
		  if (emit_absorption < emit_abs_noise){ //eliminate values caused by noise
		  emit_absorption =0;
		  }
		  
	      A1 *= c;		//sensitizer concentration
	      emit_absorption *= emit_concentration;
		  }

	      dx = gasdev ();
	      dy = gasdev ();
	      dz = gasdev ();

	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dx = dx / dr;
	      dy = dy / dr;
	      dz = dz / dr;
	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dst = -1.0 * log (ran1 ()) / log (10.0) / (A1+emit_absorption); //beer-lambert law

	      z = d * (1.0 * i - 0.5) / (1.0 * bin_N) + dz * dst;

	      if (z < 0.0)
		{		//if we escaped
		  if (ran1()< pr_absorbed(lambda,min_wavelength,sctA))//use tauc theory to test if photon is absorbed
		    {		//photon energy must be absorbed by solar cell
		      esc += 1;
			  
		    }
		}
	      else if (z < d)
		{
		  h = z / d * bin_N + 1;//which bin did absorption occur in?  cast float to int
			if(ran1()>(emit_absorption/(A1+emit_absorption))){
		  		new_reabs_bin[h] += 1;
				
			 } else {
				new_singlet_bin[h]+= 1;
				
			}
		}

	      if (z > d && BR == 1)//reflected emission 
		{
		  if (LA == 1)
		    {
		      dx = gasdev ();
		      dy = gasdev ();
		      dz = gasdev ();

		      dr = sqrt (dx * dx + dy * dy + dz * dz);
		      dx = dx / dr;
		      dy = dy / dr;
		      dz = dz / dr;
		      dr = sqrt (dx * dx + dy * dy + dz * dz);
		    }

		  dst = -1.0 * log (ran1 ()) / log (10.0) /(A1+emit_absorption);

		  z = d - sqrt (dst * dz * dst * dz);	//increment z
		//distance above the mirror that the photon was absorbed

		  if (z < 0.0)
		    {		//if we escaped
		      if (ran1()< pr_absorbed(lambda,min_wavelength,sctA))//use tauc theory to test if photon is absorbed
			{	//photon energy must be absorbed by solar cell
			  esc += 1;
			  
			}
		    }
		  else
		    {
		      h = z / d * bin_N + 1;
			if(ran1()>(emit_absorption/(A1+emit_absorption))){
                    		  //photon was absorbed by sensitizer
                      		new_reabs_bin[h] += 1;
							
                  	} else {
                      		//photon was absorbed by emitter, singlet is excited
                  		new_singlet_bin[h]+= 1;
						
			}
                    }           //reabsorption and escapage is fractional photons, more efficient that way for calcualtion
		}
	    }
	}
	
	}				//end iterate
  //current mAcm-2: bins: depth (cm): k1: k2: eta_c: sensitizer concentration: emitter concentration: solar concentration factor: bandgap: solar cell absorbsion constant: delta_E: temperature
 printf
    ("%1.6e\t%d\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\n",
     esc / run_time * 1.602e-16, bin_N, d, k1, k2, eta_c, c, emit_concentration,  C,
     min_wavelength, sctA, delta_E, temperature);
  
  return 0;
}



	