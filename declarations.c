
//initialize variables
  FILE *input;
  const unsigned int commandlineargs=11; //including name of program

//solar spectru
  double *spec_x, *spec_y, *spec_y2;

//sensitizer absorption spectrum
  double *abs_x, *abs_y, *abs_y2;
//emitter_absorption spectrum
  double *emit_abs_x, *emit_abs_y, *emit_abs_y2;

//emission spectrum
  double *emi_x, *emi_y, *emi_y2;

  double lambda;//random sample wavelength 
  double c;//sensitizer concentration
  double  A1; //interpolated absorption of sensitizer
  double emit_concentration=0;
  double emit_absorption=0;
  double d, z, dst;	//base calculations on a 1m2 concentrator
  double total_energy = 0, run_time = 0;
  double dx, dy, dz, dr;
  double k1, k2, f2;		//cm3/s
  double k_phiS, NT, esc;//arithmetic helpers 
  double  C;//solar concentration 
  double eta_c;//proportion of triplet triplet annihilation events that produce a singlet 

  double min_wavelength, max_wavelength;//define spectral window

  int i, j, h;
  int N_phot;

  srand(time(NULL));//seed random number generator using the time

  int spec_N, abs_N, emi_N, emit_abs_N;//number of items in data file
  int bin_N, *depth_bin, *singlet_depth_bin; 
  double *reabs_bin, *new_reabs_bin; //store reabsoprtion by sensitizer
  double *singlet_bin, *new_singlet_bin; //store reabsorption by emitter singlet
  const int BR = 1;//use back reflector 
  const int  LA = 1; //if true use lambertian reflector, otherwise use specular
  int iterate;
  int reabsorptioncycles = 5; //number of loops over the reabsorption calculator: recommend 5
