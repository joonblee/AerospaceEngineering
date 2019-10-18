# source respective.sh ${MODE_NUMBER} ${phi1} ${MASS_SIGN1} ${phi2} ${MASS_SIGN2}

for ((i=2; i<6; i++)); do
  MODE_NUMBER=$i
  source respective.sh $MODE_NUMBER + +
  source respective.sh $MODE_NUMBER + -
done
