#/bin/sh

set -eu
while getopts m OPT
  do 
    ./configure
    make 
    clear
  done

./build/cadical -q -t 60 --compact=false ../instances/sc2021/b04_s_unknown_pre.cnf.xz
