#/bin/bash
path=$(pwd)
if [ ! -d "~/graspitmodified-build" ]; then
  mkdir ~/graspitmodified-build
fi
cd ~/graspitmodified-build

if [ ! -d "graspit" ]; then
  mkdir graspit
fi
cd graspit

export COINDIR=$(pwd)/../coin/install-custom
cmake $path/graspit -DQHULL_ROOT=$(pwd)/../qhull/install-custom/usr/local
make
