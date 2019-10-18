# source respective.sh ${MODE_NUMBER} ${phi1} ${MASS_SIGN1} ${phi2} ${MASS_SIGN2}

#for ((i=2; i<4; i++)); do
#  MODE_NUMBER=$i
#  source respective.sh $MODE_NUMBER pi/3 +
#  source respective.sh $MODE_NUMBER pi/3 -
#  source respective.sh $MODE_NUMBER pi/3 + 4/3*pi +
#  source respective.sh $MODE_NUMBER pi/3 + 4/3*pi -
#done

sum() {
python - <<END
result=$1+$2;
print result;
END
}

divide() {
python - <<END
result=float($1)/$2;
print result;
END
}

for (( i = 2; i < 4; i++ )); do
  MODE_NUMBER=$i
  theta1=$(divide 0 18)
  #source respective.sh $MODE_NUMBER ${theta1}*pi +
  source respective.sh $MODE_NUMBER ${theta1}*pi -

  #for (( j = 0; j < 13; j++ )); do
  #  Dtheta=$(divide $j 6)
  #	 theta2=$(sum ${theta1} ${Dtheta})

  #  source respective.sh $MODE_NUMBER ${theta1}*pi + ${theta2}*pi +
  #  source respective.sh $MODE_NUMBER ${theta1}*pi + ${theta2}*pi -
  #done
done


