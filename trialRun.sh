#/bin/sh

set -eu
while getopts m OPT
  do 
    ./configure
    make 
    clear
  done

./build/cadical -q -t 60 --compact=false ../SAT_instances/sc2021/ktf_TF-3.tf_3_0.02_24.cnf.xz
