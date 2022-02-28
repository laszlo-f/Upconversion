
//initialize variables
  FILE *input;
  const unsigned int commandlineargs=13; //including name of program

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
  double k2;		//  cm3/s
  double k1, f2;		//  /s
  double k_phiS, NT, esc;//arithmetic helpers 
  double  C;//solar concentration 
  double eta_c;//proportion of triplet triplet annihilation events that produce a singlet 

  double min_wavelength;//define spectral window

  int i, j, h;
  long int N_phot;

  srand(time(NULL));//seed random number generator using the time

  int spec_N, abs_N, emi_N, emit_abs_N;//number of items in data file
  int bin_N, *depth_bin, *singlet_depth_bin; 
  double *reabs_bin, *new_reabs_bin; //store reabsoprtion by sensitizer
  double *singlet_bin, *new_singlet_bin; //store reabsorption by emitter singlet
  const int BR = 1;//use back reflector 
  const int  LA = 1; //if true use lambertian reflector, otherwise use specular
  int iterate;
  int reabsorptioncycles = 5; //number of loops over the reabsorption calculator: recommend 5
  
  double max_emitter_absorption, max_absorption, max_emission ,min_emitter_absorption, min_absorption, min_emission;//largest/smallest wavelenghth in a file

  double temperature; //usually 300 K.
  double delta_E; //energy difference between the sensitizer triplet and the emitter triplet. Typically between 0.3 and -0.1 eV.
  const double kB = 8.61733e-5; //Boltzmann constant.  (eV/K)
 
  double sctA;//(eV^1/2) product of thickness of solar cell(m) and square root of proportionality constant A from Tauc model
  //determines strength of solar cell absorption
  
  double pr_absorbed(double emission_wavelength, double short_wavelength,double sctA);//calculates the probability of photon being absorbed by solar cell
  
   const double emit_abs_noise = 1000; //amount of noise for emitter absorption 
  //int l; // index for loops not always used
  //double sum_emit_ab, sum_ab,sum_emit_abn, sum_abn; // for calcultating  emitter absorbtion and sensitizer absorbtion n is for after reflection fpr absorption
  //double sum_emit_reab, sum_reab;// for calcultating  emitter absorbtion and sensitizer absorbtion r is for after reflection for reabsorption by emitter
  //double sum_emit_reabr, sum_reabr; // for calcultating  emitter absorbtion and sensitizer absorbtion r is for after reflection for reabsorption by emitter
  //double sum_emit_reabs, sum_reabs,sum_emit_reabrs, sum_reabrs; // for calcultating  emitter absorbtion and sensitizer absorbtion r is for after reflection for reabsorption by sensitizer
  //double sum_reab_bin,sum_singlet_bin;//photons in reabsoprtion bin
  //double rando;
  
  