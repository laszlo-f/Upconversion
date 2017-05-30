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

//set up absorption spectrum of sensitizer
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
  
 //read emitter's absorption spectrum 
  input = fopen ("emitter_absorption.txt", "r");
  fscanf (input, "%d\n", &emit_abs_N);
  emit_abs_x = dvector (1, emit_abs_N);
  emit_abs_y = dvector (1, emit_abs_N);
  emit_abs_y2 = dvector (1, emit_abs_N);
  for (i = 1; i <= abs_N; i++)
    {
      fscanf (input, "%lf", &emit_abs_x[i]);
      fscanf (input, "%lf", &emit_abs_y[i]);
    }
  fclose (input);
  spline (emit_abs_x, emit_abs_y, emit_abs_N, 1e30, 1e30, emit_abs_y2);
  //end read emitter's absorption spectrum



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
  emit_concentration = atof(argv[7]);
  C = atof(argv[8]);			//solar concentration factor
  min_wavelength = atof(argv[9]);			//short wavelength cutoff
  max_wavelength = atof(argv[10]);			//long wavelength cutoff
  
  depth_bin = ivector (1, bin_N);
  singlet_depth_bin = ivector (1, bin_N);
  //do not define bin_N after you use it in ivector() bin_N = 5000;
  reabs_bin = dvector (1, bin_N);
  new_reabs_bin = dvector (1, bin_N);
  singlet_bin = dvector (1, bin_N);
  new_singlet_bin = dvector (1, bin_N);

