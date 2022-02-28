//function that returns the probability of a photon being absorbed by the solar cell using Tauc model
//emission_wavelength: the wavlength of the photon passing through the solar cell
//short_wavelength: the band gap of the solar cell 
//sctA:constant that is the product of the thickness of the solar cell and the square root of A (eV)^1/2
double pr_absorbed(double emission_wavelength, double short_wavelength, double sctA) {
	
	const double plancks_constant=4.135667662e-15;//eVs
    const double speed_light=2.99792458e8;//(m/s)
    
	
	double alpha;//absorption coefficient divided by the square root of proportionality constant "square root of A"
    double pr; //probability of photon being absorbed by solar cell, equals 1 - T_lambda
    double T_lambda; //probability of transmittance
	
	if (emission_wavelength <= short_wavelength) {//cant take the square root of a negative number
	alpha=sqrt(emission_wavelength*(1e-9)/(plancks_constant*speed_light)*(1-(emission_wavelength/short_wavelength))); //Uses Tauc model, emission_wavelength is converted from nanometres to metres
	T_lambda=pow(10,-(sctA*alpha));
	pr=1-T_lambda;
  }
  else {
	pr=0;//assuming square root af negative number is zero so T_lambda=10^0=1
  }
  return pr;
}