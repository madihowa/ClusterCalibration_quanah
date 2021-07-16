#!/usr/bin/env bash

# setting XDG dir and permissions -> MIH
export XDG_RUNTIME_DIR=/lustre/scratch/madihowa/xdg
chown -R madihowa $XDG_RUNTIME_DIR
chmod -R 0700 $XDG_RUNTIME_DIR
chmod 000700 $XDG_RUNTIME_DIR 

# setting QT var -> MIH
export QT_QPA_PLATFORM='offscreen'


# load runType as a variable
RCONTENT=$1

interactive=
DATE=$(date '+%Y-%m-%d')
RUN=-1
counter=1
RUNPYTHON="1"
echo $date

while [ "$1" != "" ]; do
  case $1 in
    -d | --date )   shift
    DATE=$1
    ;;
    -r | --run )    shift
    RUN=$1
    ;;
  esac
  shift
done

echo "$date"

for FILE in `ls -l`
do
  if test -d $FILE
  then
    echo $FILE
    if [[ "$FILE" == *"$DATE"* ]]
    then
      if [[ "$FILE" == *run_"$RUN"* ]]
      then
          for FILE2 in `ls $FILE -l`
          do
            if [[ "$FILE2" == *root* ]]
            then
              RUNPYTHON=0
              echo "Found Root FILE not rerunning network"
              rm $FILE/$FILE2 #Remove root file then regenerate it later.
            fi
          done
      fi
      ((counter++))
    fi

  fi
done
#printf "\033c"
clear
# echo $RUNPYTHON

if (($RUNPYTHON == "1"));
then
  # echo $RUN
  python Master.py -r $RUN
  #python Master.py -r $RUN -d $DATE
fi

if (("$RUN" == "-1")); then
  #echo "New Run"
  RUN=$counter
fi

StringRes=Results_
StringRun=_run_


# Comment Start -> MIH
#if (($RUNPYTHON == "1"));
#then
#  mv callback_history.csv $DIRECTORY
#fi
# Comment End -> MIH
# Reason Comment -> implemented this feature in the .py code itself

#echo $DIRECTORY
#echo $DATE

DIRECTORY="$StringRes$DATE$StringRun$RUN"
cd $DIRECTORY

echo $RCONTENT >> runType.txt

root -l ../runT_new.C

cd ..
