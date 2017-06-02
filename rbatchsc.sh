#runs R script to generate a bunch of configurations, then runs simulation
#supercomputer version
#only runs one configuration based on the ARRAYID

starttime=$(date +%d%m%H%M%S);
echo "Current (mA/cm2)	bins	depth (cm)	k1	k2	eta_c	sensitizer concentration (M)	emitter concentration (M)	solar concentration factor	minwavelength	maxwavelength" > $starttime.out.dat

#unnecessary loop goes over the row selected by PBS_ARRAYID
while read p; do
  exec ./a.out $p >>$starttime.out.dat &
done < <(Rscript generator.R | sed -n "$PBS_ARRAYID"'p')

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
