#runs program to generate a bunch of configurations, then runs simulation

starttime=$(date +%d%m%H%M%S);
echo "Current (mA/cm2), bins, depth(cm), k1, k2, eta_c, sensitizer concentration(M), emitter concentration(M), solar concentration factor, minwavelength(nm), solar cell constant, delta_E, temperature" > $starttime.out.dat

while read p; do
  exec ./a.out $p >>$starttime.out.dat &
done < <(./generator.out)

echo "Started Calculation"

FAIL=0


for job in `jobs -p`
do
    wait $job || let "FAIL+=1"
done

echo $FAIL

if [ "$FAIL" == "0" ];
then
echo "Success"
else
echo "FAIL! ($FAIL)"
fi
