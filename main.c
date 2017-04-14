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

//initialize variables
  FILE *input;
  const unsigned int commandlineargs=10; //including name of program

//solar spectru
  double *spec_x, *spec_y, *spec_y2;

//sensitizer absorption spectrum
  double *abs_x, *abs_y, *abs_y2;

//emission spectrum
  double *emi_x, *emi_y, *emi_y2;

  double lambda;//random sample wavelength 
  double c;//sensitizer concentration
  double  A1; //interpolated absorption
  double d, z, dst;	//base calculations on a 1m2 concentrator
  double total_energy = 0, run_time = 0;
  double dx, dy, dz, dr;
  double k1, k2, f2;		//cm3/s
  double k_phiS, NT, em, esc, C, eta_c, *new_reabs_bin, T;

  double min_wavelength, max_wavelength;//define spectral window

  int i, j, h, abs;
  int N_phot;
  long k = (short)time(NULL);//seed random number generator using the time

  int spec_N, abs_N, emi_N;//number of items in data file
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
    }
  fclose (input);

  spline (emi_x, emi_y, emi_N, 1e30, 1e30, emi_y2);
//end set up emission spectrum

//Get command line arguments
  if(argc!=commandlineargs)
  {//test number of arguments
	printf("Wrong number of command line arguments: %d\n",argc);
	exit(-1);
  }
  bin_N = atoi(argv[1]);			//number of bins
	//laszlo recomends 1e5, but make sure bin size is << absorption length
  d = atof(argv[2]); //depth of the cell
  k1 = atof(argv[3]);			//triplet nonradiative decay
  k2 = atof(argv[4]);			//annihilation rate cm3/s 1.7e-13 for rubrene
  eta_c = atof(argv[5]);			//proportion of annihilation events which lead to the singlet state
  c = atof(argv[6]);			//concentration of sensitizer
  C = atof(argv[7]);			//solar concentration factor
  min_wavelength = atof(argv[8]);			//short wavelength cutoff
  max_wavelength = atof(argv[9]);			//long wavelength cutoff
  
  depth_bin = ivector (1, bin_N);
  //do not define bin_N after you use it in ivector() bin_N = 5000;
  reabs_bin = ivector (1, bin_N);
  new_reabs_bin = dvector (1, bin_N);

  abs = 0;
  esc = 0;
  em = 0;


  N_phot = 1e4*bin_N;			//check me for convergence - recommend at least 10^4 bin_N
  BR = 1;
  LA = 1;
  T = 1.0;			//solar cell/front surface transparency


      total_energy = 0.0;
      for (i = 1; i <= bin_N; i++)
      {
	  depth_bin[i] = 0;
	}
      for (i = 1; i <= N_phot; i++)	//number of incident photons
	{
	  splint (spec_x, spec_y, spec_y2, spec_N, ran1 (&k), &lambda);	//choose a wavelength

	  total_energy += 6.626e-34 * 2.9979e8 / lambda / 1e-9;

	  if (lambda > min_wavelength && lambda < max_wavelength && ran1 (&k) < T)	//transmitted to sample and absorbable
	    {
	      splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
	      A1 *= c;		//concentration

	      z = -1.0 * log (ran1 (&k)) / log (10.0) / A1;	//check to see if photon is absorbed

	      if (z < d && z > 0.0)	//incoming photon is absorbed
		{
		  abs++;

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
		      j = z / d * bin_N + 1;
		      depth_bin[j] += 1;
		    }

		}		//if absorbed
	    }
	}			//i

      

      //printf ("total energy was %1.4e J\n", total_energy);
      //printf ("run time was therefore %1.4e s\n", run_time = total_energy / 100e-3 / C);

      //one sun is 100 mW/cm^2
      run_time = total_energy / 100e-3 / C;


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
		  splint (abs_x, abs_y, abs_y2, abs_N, lambda, &A1);	//get absorption coefficienct
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


	}			//end iterate
	//current, mAcm-2: bins: depth (cm): k1: k2: eta_c: sensitizer concentration: solar concentration factor
      printf("%1.6e\t%d\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\n",  esc / run_time * 1.602e-16, bin_N,d,k1,k2,eta_c,c,C,min_wavelength,max_wavelength);
  


  return 0;
}
