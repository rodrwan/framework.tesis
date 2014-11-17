#/bin/bash

MODEL=$1
PARTS=$2
nohup matlab -nodesktop -nosplash -nodisplay -r "training('$MODEL', $PARTS);" -logfile "logfile_$MODEL.out" </dev/null &
