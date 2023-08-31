#!/bin/sh

# Necessary input for CMake
BOOSTPATH=${HOME}/Software/boost_1_67_0/

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
PROJECTNAME=$(basename "$SCRIPTPATH")
#BUILDTYPE=Debug
#BUILDTYPE=RelWithDebInfo
BUILDTYPE=Release

# Clear build directory and cd to it
rm -rf ${SCRIPTPATH}/${BUILDTYPE}
mkdir ${SCRIPTPATH}/${BUILDTYPE}
cd ${SCRIPTPATH}/${BUILDTYPE}
# set -x displays cmake command. Brackets create subshell so that there is no need to call set +x
(set -x; cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=${BUILDTYPE} -DBOOST_ROOT=${BOOSTPATH} -DCMAKE_ECLIPSE_GENERATE_SOURCE_PROJECT=TRUE -DCMAKE_ECLIPSE_MAKE_ARGUMENTS=-j6 ../)


