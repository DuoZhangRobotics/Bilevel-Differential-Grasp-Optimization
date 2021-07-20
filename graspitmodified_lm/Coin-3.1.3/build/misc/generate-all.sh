#!/bin/sh

SolutionName=coin3

cd $(dirname $0)/../..
export COINDIR=$(cygpath -w $(pwd))
cd build/misc
for i in $(seq 6 9)
do
  ./generate.bat msvc$i ${SolutionName}
done

#msvc6 is only our starting point, but cannot be built.
#Uncomment to keep the files.
#rm -rf msvc6
