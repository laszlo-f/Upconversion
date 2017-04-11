#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <time.h> //random seed


#include "gasdev.c"
#include "nrutil.c"
#include "splint.c"
#include "spline.c"
#include "ran1.c"

int
main (int argc, char *argv[])
{
  FILE *input, *output;
  FILE *out_incoming, *out_absorbed, *out_emitted, *out_detected, *out_depth;
    //*out_loop;
  //double x_dim, y_dim, z_dim;
  const unsigned int commandlineargs=7; //including name of program
  double *spec_x, *spec_y, *spec_y2;
  double *abs_x, *abs_y, *abs_y2;
  double *emi_x, *emi_y, *emi_y2;
  double lambda, c, A1, A2;// ratio;
  //double random, lambda, c, A1, A2;// ratio;
  double stokes;// shift;
  //double Phi_F, thresh, stokes;// shift;
  double d, z, dst;	//base calculations on a 1m2 concentrator
  //double G, d, x, y, z, n, dst;	//base calculations on a 1m2 concentrator
  double total_energy = 0, run_time = 0;
  double dx, dy, dz, dr;
  double k1, k2, f2;		//cm3/s
  double k_phiS, NT, em, esc, C, eta_c, *new_reabs_bin, T;

  int i, j, h, abs;// trp, wfg, 
  int N_phot;// N_curr, W_curr;
  long k = (short)time(NULL);//seed random number generator using the time
	//printf("Random seed %ld\n",k);

  //int reabs;
  int spec_N, abs_N, emi_N;// alive, out, cheat;
  int bin_N, *depth_bin, *reabs_bin, BR = 1, LA = 1, iterate;
  int reabsorptioncycles = 5; //number of loops over the reabsorption calculator: recommend 5


//AM1.5G data file has 1813
  input = fopen ("spectrum.txt", "r");
  fscanf (input, "%d\n", &spec_N);
  spec_x = dvector (1, spec_N);
  spec_y = dvector (1, spec_N);
  spec_y2 = dvector (1, spec_N);


  for (i = 1; i <= spec_N; i++)
    {
      fscanf (input, "%lf", &spec_x[i]);
      fscanf (input, "%lf", &spec_y[i]);
      //printf("%d\t%1.4f\t%1.3f\n",i,spec_x[i],spec_y[i]);
    }
  fclose (input);

  spline (spec_x, spec_y, spec_N, 1e30, 1e30, spec_y2);
//end set up solar spectrum

//set up absorption spectrum
  input = fopen ("absorption.txt", "r");
  fscanf (input, "%d\n", &abs_N);
  abs_x = dvector (1, abs_N);
  abs_y = dvector (1, abs_N);
  abs_y2 = dvector (1, abs_N);

  for (i = 1; i <= abs_N; i++)
    {
      fscanf (input, "%lf", &abs_x[i]);
      fscanf (input, "%lf", &abs_y[i]);
      //printf("%d\t%1.4f\t%1.3f\n",i,abs_x[i],abs_y[i]);
    }
  fclose (input);

  spline (abs_x, abs_y, abs_N, 1e30, 1e30, abs_y2);
//end set up absorption spectrum


//set up emission spectrum
  input = fopen ("emission.txt", "r");
  fscanf (input, "%d\n", &emi_N);
  emi_x = dvector (1, emi_N);
  emi_y = dvector (1, emi_N);
  emi_y2 = dvector (1, emi_N);

  for (i = 1; i <= emi_N; i++)
    {
      fscanf (input, "%lf", &emi_x[i]);
      fscanf (input, "%lf", &emi_y[i]);
      //printf("%d\t%1.4f\t%1.3f\n",i,emi_x[i],emi_y[i]);
    }
  fclose (input);

  spline (emi_x, emi_y, emi_N, 1e30, 1e30, emi_y2);
//end set up emission spectrum

//Get command line arguments
  if(argc!=commandlineargs)
  {//test number of arguments
	printf("Wrong number of command line arguments\n");
	exit(-1);
  }
  printf("%s\n",argv[1]);
  bin_N = atoi(argv[1]);			//number of bins
	//laszlo recomends 1e5, but make sure bin size is << absorption length
  d = atof(argv[2]); //depth of the cell
/*
  k1 = ;			//triplet nonradiative decay
  k2 = ;			//annihilation rate cm3/s 1.7e-13 for rubrene
  eta_c = ;			//proportion of annihilation events which lead to the singlet state
  d = ;			//thickness in cm //only used if the for loop is disabled
  c = ;			//concentration of sensitizer
  C = ;			//solar concentration factor
*/
  for(i=1;i<=commandlineargs;i++)
  {
	printf("%d\n",i);
  }
  
  //bin_N = 1e1;
  depth_bin = ivector (1, bin_N);
  reabs_bin = ivector (1, bin_N);
  new_reabs_bin = dvector (1, bin_N);

  abs = 0;
  esc = 0;
  //trp = 0;
  //wfg = 0;
  //reabs = 0;
  em = 0;
  //N_curr = 0;
  //W_curr = 0;

  //cheat = 0;


  output = fopen ("parameter.dat", "w");

  out_incoming = fopen ("incoming.dat", "w");
  out_absorbed = fopen ("absorbed.dat", "w");
  out_emitted = fopen ("emitted.dat", "w");
  out_detected = fopen ("detected.dat", "w");
  out_depth = fopen ("depth.dat", "w");



  N_phot = 1e4*bin_N;			//check me for convergence - recommend at least 10^4
  c = 1e-3;			//concentration of sensitizer
//  d = 0.0008;			//thickness in cm //only used if the for loop is disabled
  //do not define bin_N after you use it in ivector() bin_N = 5000;
  BR = 1;
  LA = 1;
  stokes = 0.0;			//don't use this
  C = 1.0;			//solar concentration factor
  eta_c = 1.0;			//proportion of annihilation events which lead to the singlet state
  T = 1.0;			//solar cell/front surface transparency
  k1 = 1e4;			//triplet nonradiative decay
  k2 = 1.7e-12;			//annihilation rate cm3/s 1.7e-13 for rubrene


  //out_loop = fopen ("looper.dat", "w");
  //
  //loop over depth
  //for (d = 0.0001; d <= .01; d *= 10)
    {
      total_energy = 0.0;
      for (i = 1; i <= bin_N; i++)
      {
	  depth_bin[i] = 0;
	}
      for (i = 1; i <= N_phot; i++)	//number of incident photons
	{
	  splint (spec_x, spec_y, spec_y2, spec_N, ran1 (&k), &lambda);	//choose a wavelength
	  //fprintf (out_incoming, "%1.1f\n", lambda);

	  //printf("%1.1f\n",lambda);
	  total_energy += 6.626e-34 * 2.9979e8 / lambda / 1e-9;

	  if (lambda > 600.0 && lambda < 1300.0 && ran1 (&k) < T)	//transmitted to sample and absorbable
	    {
	      splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	      splint (abs_x, abs_y, abs_y2, abs_N, lambda - stokes, &A2);	//get absorption coefficienct
	      A1 = (A1 + A2) / 2;
	      A1 *= c;		//concentration

	      z = -1.0 * log (ran1 (&k)) / log (10.0) / A1;	//check to see if photon is absorbed

	      if (z < d && z > 0.0)	//incoming photon is absorbed
		{
		  abs++;
		  //fprintf (out_absorbed, "%1.1f\n", lambda);
		  //printf("%1.1f is absorbed\n",lambda);

		  j = z / d * bin_N + 1;
		  depth_bin[j] += 1;

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

		  dst = -1.0 * log (ran1 (&k)) / log (10.0) / A1;

		  z = d - sqrt (dst * dz * dst * dz);	//increment z

		  if (z > 0.0 && z < d && BR == 1)	//reabsorb off the back reflector
		    {
		      abs++;
		      //fprintf (out_absorbed, "%1.1f\n", lambda);
		      //printf("%1.1f is absorbed\n",lambda);
		      //fprintf(out_depth,"%1.4e\n",z);
		      j = z / d * bin_N + 1;
		      depth_bin[j] += 1;
		    }

		}		//if absorbed
	    }
	}			//i

      for (i = 1; i <= bin_N; i++)
	{
	  //fprintf (out_depth, "%d\t%d\n", i, depth_bin[i]);
	}

      //printf ("Was that fun?\n");
      //printf ("total energy was %1.4e J\n", total_energy);
      //printf ("run time was therefore %1.4e s\n", run_time = total_energy / 100e-3 / C);
      //
      //one sun is 100 mW/cm^2
      run_time = total_energy / 100e-3 / C;

      //don't do assignments in print statements
      //printf ("... so at front of cuvette, we absorbed photons at %1.4e cm-3s-1\n",
	 //k_phiS = depth_bin[1] / run_time / (d / bin_N));
	 k_phiS = depth_bin[1] / run_time / (d / bin_N);
      //printf ("... triplet concentration is about %1.4e cm-3, which is about %1.4e M\n",
	// NT = k_phiS / k1, k_phiS / k1 * 1000 / 6.022e23);

      NT = (-k1 + sqrt (k1 * k1 + 4 * k_phiS * k2)) / 2 / k2;//not used?

      //printf ("solving a quadratic,  we get that the triplet concentration at the front of the cuvette is %1.4e cm-3\n", NT);
      //printf ("f2 at cuvette front is approximately %1.3e\n", k2 * NT / (k2 * NT + k1));

//NOW GET A SELF CONSISTENT SOLUTION



      for (iterate = 1; iterate <= reabsorptioncycles; iterate++)
	{
	  esc = 0;
	  em = 0;
	  //printf ("%d\t\t%1.3f\t%d\n", depth_bin[1], new_reabs_bin[1], reabs_bin[1]);

	  for (i = 1; i <= bin_N; i++)
	    {
	      reabs_bin[i] = new_reabs_bin[i];
	      new_reabs_bin[i] = 0.0;
	    }

	  for (i = 1; i <= bin_N; i++)
	    {
	      k_phiS = (depth_bin[i] + reabs_bin[i]) / run_time / (d / bin_N);	//add in reabsorption term
	      NT = (-k1 + sqrt (k1 * k1 + 4 * k_phiS * k2)) / 2 / k2;
	      f2 = k2 * NT / (k2 * NT + k1);
	      em += eta_c * depth_bin[i] * f2 * 0.5;

	      //propagate photons out
	      for (j = 1; j <= depth_bin[i] + reabs_bin[i]; j++)
		{
		  splint (emi_x, emi_y, emi_y2, emi_N, ran1 (&k), &lambda);	//choose a wavelength
		  //fprintf (out_emitted, "%1.1f\n", lambda);
		  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
		  splint (abs_x, abs_y, abs_y2, abs_N, lambda - stokes, &A2);	//get absorption coefficienct
		  A1 = (A1 + A2) / 2;
		  A1 *= c;	//concentration

		  dx = gasdev (&k);
		  dy = gasdev (&k);
		  dz = gasdev (&k);

		  dr = sqrt (dx * dx + dy * dy + dz * dz);
		  dx = dx / dr;
		  dy = dy / dr;
		  dz = dz / dr;
		  dr = sqrt (dx * dx + dy * dy + dz * dz);

		  dst = -1.0 * log (ran1 (&k)) / log (10.0) / A1;

		  z = d * (1.0 * i - 0.5) / (1.0 * bin_N) + dz * dst;

		  if (z < 0.0)
		    esc += eta_c * f2 * 0.5;
		  else if (z < d)
		    {
		      h = z / d * bin_N + 1;
		      new_reabs_bin[h] += eta_c * f2 * 0.5;
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

		      dst = -1.0 * log (ran1 (&k)) / log (10.0) / A1;

		      z = d - sqrt (dst * dz * dst * dz);	//increment z
		      if (z < 0.0)
			esc += eta_c * f2 * 0.5;
		      else
			{
			  h = z / d * bin_N + 1;
			  new_reabs_bin[h] += eta_c * f2 * 0.5;
			}	//reabsorption and escapage is fractional photons, more efficient that way for calcualtion
		    }

		}
	    }

//printf("*****\nInternal Quantum Efficiency is %1.4f out of 0.5000\n",em/abs);
	  //printf ("External Quantum Efficiency is %1.4f out of 0.5000\n", esc / abs);

	}			//end iterate
      printf ("d is %1.4e, Emitted Current is %1.4e mAcm-2\n", d, esc / run_time * 1.602e-16);
      //fprintf (out_loop, "%1.4e\t%1.4e\n", d, esc / run_time * 1.602e-16);
    }				// end d loop

  fclose (output);
  fclose (out_incoming);
  fclose (out_absorbed);
  fclose (out_emitted);
  fclose (out_detected);
  return 0;
}
