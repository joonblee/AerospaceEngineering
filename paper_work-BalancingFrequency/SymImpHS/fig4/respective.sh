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

if [ -d "$DIR" ]; then # -f: regular file, -e: file exists
  echo "Directory '$DIR' already exists."
else
  echo "There is no '$DIR' directory."
  mkdir $DIR
  echo "'$DIR' is created."
fi

if [ -z "$2" ]; then
  echo 'If you want to run this code, try below:'
  echo 'source respective.sh $1 $2 ($3) ($4 $5)'
  echo '$1: mode number'
  echo '$2: theta(1), $4: theta(2) (ex. 0, pi)'
  echo '$3: mass sign(1), $5: mass sign(2) [+: attahced mass, -: hole]'
  echo 'You can omit $4 and $5, your model has sigle attached mass'
else
  rm run.m

  touch temp

  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '%% options                      %%' >> temp
  echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' >> temp
  echo '' >> temp
  echo "n=$1 %% mode number" >> temp
  echo '' >> temp

  # imperfection #
  MASS=0.01
  PHI=pi/2
  # ------------ #

  if ! [ -z "$4" ]; then
    echo "phi_i   = [${PHI}    ${PHI}   ] %% imperfections" >> temp
    echo "theta_i = [$2        $4       ]" >> temp
    echo "m_i     = [$3${MASS} $5${MASS}]" >> temp
    echo "SIGN_i  = ['$3'      '$5'     ];"     >> temp
  else
    echo "phi_i   = [${PHI}   ]" >> temp
	echo "theta_i = [$2       ]" >> temp
	echo "m_i     = [$3${MASS}]" >> temp
	echo "SIGN_i  = ['$3'     ]" >> temp
  fi

  echo '' >> temp
  echo 'phi   = pi/2;' >> temp
  echo 'theta = [0:359]/180*pi;' >> temp
  echo '' >> temp
  echo 'a=10; % radius' >> temp
  echo 'D=.1; % mode shape amplitude' >> temp
  echo '' >> temp

  cat temp main.m >> run.m
  rm temp

  echo "Run code"
  echo "matlab run ($1 $2 $3 $4 $5)"


  for ((i=0; i<9999; i++)); do 
    if [ ! -e "log$i" ]; then # -f: DIRs
      matlab run > log$i < /dev/null
      break
    fi
  done

  rm run.m
fi

echo "Finished"

