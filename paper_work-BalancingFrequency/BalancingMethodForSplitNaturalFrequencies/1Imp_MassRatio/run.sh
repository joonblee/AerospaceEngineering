#!/bin/bash
# source run.sh 

DIR=figures

if [ -d "$DIR" ]; then # -f: DIRs
  echo "'$DIR' already exists."
else
  echo "There is no $DIR directory."
  mkdir $DIR
  echo "'$DIR' is created."
fi

matlab main > log < /dev/null

echo "Finished"

