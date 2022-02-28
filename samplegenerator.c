#include <stdio.h>
#include <math.h>
// this generator file is an example.
//this sample uses a range of different sensitizer concentrations from 1 micromolar to 1mM evenly spacing datapoints for a log scale
//
int main(){
  double i;//used for iteration 

  //individual values can be changed here

  int bins=100;
  double thickness=0.1;//thickness of anabathmophore 
  double min_wavelength=480;// or bandgap(nm)
  double emit_concentration=1e-2; // (M)
  double sens_concentration=1e-3; // (M)
  double kB = 8.61733e-5; //Boltzmann constant.  (eV/K)
  double temperature =300; //usually 300 K.
  double kq=4.8e7;//quenching constant. per second per M.
  double kE1=2e3; //energy loss rate of the emitter in the absence of sensitizer. per second.
  double kS=8550; //energy loss rate of the sensitizer in the absence of emitter. per second.
  double delta_E=0.3; //energy difference between the sensitizer triplet and the emitter triplet. eV.
  double SCTA;//product of thickness (m) and square root of proportionality constant for solarcell
  double k2=4.7e-12; // triplet annihalation.  cm^3/s
  
  double Boltzmann;
  double k1_new; //It is an energy loss rate.(cm3/s)
  
  double plancks_constant=4.135667662e-15;// (eVs)
  double speed_light=2.99792458e8;//(m/s)
    
  emit_concentration= 1e-2;
  min_wavelength =480;
  for (i = 1e-6; i <=1e-3 ; i*=1.29154966501){//concentration iteration     (10^(1/9)=1.29)
	  sens_concentration=i;//reassign values when iterating of that particular value
	  SCTA=2*sqrt((0.1*min_wavelength*(1e-9)+plancks_constant*speed_light)/(min_wavelength*(1e-9)*(1-((plancks_constant*speed_light)/(0.1*min_wavelength*(1e-9)+plancks_constant*speed_light)))));
	  //to not use tauc model add line "SCTA=1e9;" here (eg any large number)
	  Boltzmann=exp(-1*(delta_E)/(kB*temperature)); 
	  k1_new=((sens_concentration*kS*Boltzmann)+(emit_concentration*(kE1+kq*sens_concentration)))/((sens_concentration*Boltzmann)+(emit_concentration));
	  printf( "%d %1.6e %1.6e %1.6e 1. %1.6e %1.6e 1. %1.6e %1.6e %1.6e %1.6e\n",bins, thickness, k1_new, k2, sens_concentration, emit_concentration, min_wavelength, SCTA, delta_E,temperature);
	  //bins	depth (cm)	k1	k2	eta_c	sensitizer concentration (M)	emitter concentration (M)	solar concentration factor	minwavelength (nm)	solar cell absorption parameter thickness*sqrt(A)	delta_E	temperature (K)

  }

  return (0);
}  
