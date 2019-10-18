#!/bin/bash
# source respective.sh ${mode_number} ${theta1} ${mass_sign1} ${theta2} ${mass_sign2}

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


DIR=figures

if [ -d "$DIR" ]; then # -f: DIRs
  echo "Directory '$DIR' already exists."
else
  echo "There is no '$DIR' directory."
  mkdir $DIR
  echo "'$DIR' is created."
fi

nTheta=2;

theta1=$(divide 0 $nTheta);
 #source respective.sh ${theta1}*pi +
 #source respective.sh ${theta1}*pi -

for ((j=0; j<(nTheta+1); j++)); do
  Dtheta=$(divide $j $nTheta);
  theta2=$(sum ${theta1} ${Dtheta});

  source respective.sh ${theta1}*pi + ${theta2}*pi +
  source respective.sh ${theta1}*pi + ${theta2}*pi -
done

echo "Finish, run.sh"

