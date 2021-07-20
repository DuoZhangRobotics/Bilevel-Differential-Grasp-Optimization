#!/bin/bash
path=$(pwd)
if [ ! -d "~/graspitmodified-build" ]; then
  mkdir ~/graspitmodified-build
fi
cd ~/graspitmodified-build

if [ ! -d "coin" ]; then
  mkdir coin
fi
cd coin

$path/Coin-3.1.3/configure --prefix=$(pwd)/install-custom
make install
