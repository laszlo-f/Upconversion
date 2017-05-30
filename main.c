#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <time.h>		//random seed

#include "gasdev.c"
#include "nrutil.c"
#include "splint.c"
#include "spline.c"
#include "ran1.c"

int
main (int argc, char *argv[])
{

#include "declarations.c"
#include "read_input.c"

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
      splint (spec_x, spec_y, spec_y2, spec_N, ran1 (&k), &lambda);	//choose a wavelength
      total_energy += 6.626e-34 * 2.9979e8 / lambda / 1e-9;
      if (lambda > min_wavelength && lambda < max_wavelength)	//test if photon is in the user selected window of the spectrum
	{
	  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	  splint (emit_abs_x, emit_abs_y, emit_abs_y2, emit_abs_N, lambda, &emit_absorption);	//get absorption coefficienct
	  A1 *= c;		//sensitizer concentration
	  emit_absorption *= emit_concentration;

	  z = -1.0 * log (ran1 (&k)) / log (10.0) / (A1+emit_absorption);	//check to see if photon is absorbed
	  if (z < d && z > 0.0)	//incoming photon is absorbed
	    {
	      j = z / d * bin_N + 1;
	      if(ran1(&k)>(emit_absorption/(A1+emit_absorption))){
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
		dx = gasdev (&k);
	      else
		dx = 0.0;
	      if (LA == 1)
		dy = gasdev (&k);
	      else
		dy = 0.0;
	      dz = gasdev (&k);

	      //normalize
	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dx = dx / dr;
	      dy = dy / dr;
	      dz = dz / dr;
	      dr = sqrt (dx * dx + dy * dy + dz * dz);

	      dst = -1.0 * log (ran1 (&k)) / log (10.0) / (A1+emit_absorption);
	      z = d - sqrt (dst * dz * dst * dz);	//increment z

	      if (z > 0.0 && z < d && BR == 1)	//reabsorb off the back reflector
		{
		  j = z / d * bin_N + 1;
	      if(ran1(&k)>(emit_absorption/(A1+emit_absorption))){
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
	  NT = (-k1 + sqrt (k1 * k1 + 4 * k_phiS * k2)) / 2 / k2;
	  f2 = k2 * NT / (k2 * NT + k1);

	  //propagate photons out
	  for (j = 1; j <= depth_bin[i] + reabs_bin[i]; j++)
	    {
	      splint (emi_x, emi_y, emi_y2, emi_N, ran1 (&k), &lambda);	//choose a wavelength
	      splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	      splint (emit_abs_x, emit_abs_y, emit_abs_y2, emit_abs_N, lambda, &emit_absorption); //get absorption coefficient
	      A1 *= c;		//sensitizer concentration
	      emit_absorption *= emit_concentration;

	      dx = gasdev (&k);
	      dy = gasdev (&k);
	      dz = gasdev (&k);

	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dx = dx / dr;
	      dy = dy / dr;
	      dz = dz / dr;
	      dr = sqrt (dx * dx + dy * dy + dz * dz);
	      dst = -1.0 * log (ran1 (&k)) / log (10.0) / (A1+emit_absorption); //beer-lambert law

	      z = d * (1.0 * i - 0.5) / (1.0 * bin_N) + dz * dst;

	      if (z < 0.0)
		{		//if we escaped
		  if (lambda < min_wavelength)
		    {		//photon energy must be absorbed by solar cell
		      esc += eta_c * f2 * 0.5;
		    }
		}
	      else if (z < d)
		{
		  h = z / d * bin_N + 1;//which bin did absorption occur in?  cast float to int
			if(ran1(&k)>(emit_absorption/(A1+emit_absorption))){
		  		new_reabs_bin[h] += eta_c * f2 * 0.5;
			 } else {
				new_singlet_bin[h]+= eta_c * f2 * 0.5;
			}
		}

	      if (z > d && BR == 1)
		{
		  if (LA == 1)
		    {
		      dx = gasdev (&k);
		      dy = gasdev (&k);
		      dz = gasdev (&k);

		      dr = sqrt (dx * dx + dy * dy + dz * dz);
		      dx = dx / dr;
		      dy = dy / dr;
		      dz = dz / dr;
		      dr = sqrt (dx * dx + dy * dy + dz * dz);
		    }

		  dst = -1.0 * log (ran1 (&k)) / log (10.0) /(A1+emit_absorption);

		  z = d - sqrt (dst * dz * dst * dz);	//increment z
		  if (z < 0.0)
		    {		//if we escaped
		      if (lambda < min_wavelength)
			{	//photon energy must be absorbed by solar cell
			  esc += eta_c * f2 * 0.5;
			}
		    }
		  else
		    {
		      h = z / d * bin_N + 1;
			if(ran1(&k)>(emit_absorption/(A1+emit_absorption))){
                    		  //photon was absorbed by sensitizer
                      		new_reabs_bin[h] += eta_c * f2 * 0.5;
                  	} else {
                      		//photon was absorbed by emitter, singlet is excited
                  		new_singlet_bin[h]+= eta_c * f2 * 0.5;
			}
                    }           //reabsorption and escapage is fractional photons, more efficient that way for calcualtion
		}
	    }
	}
    }				//end iterate
  //current, mAcm-2: bins: depth (cm): k1: k2: eta_c: sensitizer concentration: solar concentration factor
  printf
    ("%1.6e\t%d\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\n",
     esc / run_time * 1.602e-16, bin_N, d, k1, k2, eta_c, c, C,
     min_wavelength, max_wavelength);
  return 0;
}
