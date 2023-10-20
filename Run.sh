#!/bin/sh

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
PROJECTNAME=$(basename "$SCRIPTPATH")
#BUILDTYPE=Debug
#BUILDTYPE=RelWithDebInfo
BUILDTYPE=build

cd ${SCRIPTPATH}
./${BUILDTYPE}/Data_processing | tee ${BUILDTYPE}/Log.txt

