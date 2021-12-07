#/bin/sh

set -eu
while getopts m OPT
  do 
    ./configure
    make 
    clear
  done

#./build/cadical -q -t 60 --compact=false ../instances/b04_s_unknown_pre.cnf
#./build/cadical -q -t 60 --compact=false ../instances/sc2019/ex051_9-sc2018.cnf.xz
#./build/cadical -q -t 60 --compact=false ../instances/sc2019/aes_24_4_keyfind_2-sc2013.cnf.xz
#./build/cadical -q -t 60 --compact=false ../instances/sc2019/Pb-chnl15-16_c18.cnf.xz
#./build/cadical -q -t 60 --compact=false ../instances/sc2020/170225515.cnf.xz
./build/cadical -q -t 60 --compact=false ../instances/sc2021/Kakuro-easy-051-ext.xml.hg_4.cnf.xz
