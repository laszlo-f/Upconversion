echo "Current (mA/cm2)	bins	depth (cm)	k1	k2	eta_c	sensitizer concentration (M)	solar concentration factor" > $1.out.dat

while read p; do
  exec ./a.out $p >>$1.out.dat &
done <$1

echo "Started Calculation"
