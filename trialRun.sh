#/bin/sh

set -eu
while getopts m OPT
  do 
    ./configure
    make 
    clear
  done

./build/cadical -q -t 60 --compact=false ../instances/sc2019/ex051_9-sc2018.cnf.xz
