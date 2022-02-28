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
//first line is the number of wavelengths
//must be sorted by increasing wavelength
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
  min_absorption=abs_x[1];
  max_absorption=abs_x[abs_N];
  fclose (input);
  spline (abs_x, abs_y, abs_N, 1e30, 1e30, abs_y2);
//end set up absorption spectrum
  
//first line is the number of wavelengths
//must be sorted by increasing wavelength
 //read emitter's absorption spectrum 
  input = fopen ("emitter_absorption.txt", "r");
  fscanf (input, "%d\n", &emit_abs_N);
  emit_abs_x = dvector (1, emit_abs_N);
  emit_abs_y = dvector (1, emit_abs_N);
  emit_abs_y2 = dvector (1, emit_abs_N);
  for (i = 1; i <= emit_abs_N; i++)
    {
      fscanf (input, "%lf", &emit_abs_x[i]);
      fscanf (input, "%lf", &emit_abs_y[i]);
    }
  min_emitter_absorption=emit_abs_x[1];
  max_emitter_absorption=emit_abs_x[emit_abs_N];
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
  min_emission=emi_y[1];
  max_emission=emi_y[emi_N];
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

  sctA = atof(argv[10]);//(eV^1/2) product of thickness of solar cell(m) and square root of proportionality constant A from Tauc model
  //determines strength of solar cell absorption
  delta_E = atof(argv[11]); //energy difference between the sensitizer triplet and the emitter triplet. Typically between 0.3 and -0.1 eV.
  temperature = atof(argv[12]); //usually 300 K.
  
  depth_bin = ivector (1, bin_N);
  singlet_depth_bin = ivector (1, bin_N);
  //do not define bin_N after you use it in ivector() bin_N = 5000;
  reabs_bin = dvector (1, bin_N);
  new_reabs_bin = dvector (1, bin_N);
  singlet_bin = dvector (1, bin_N);
  new_singlet_bin = dvector (1, bin_N);

