#/bin/sh

year=$1 #第一引数で"sc20xx"を指定
timeout=1200

for file in $(ls ../instances/${year}); do
	
	echo "${file}" >> ./output/${year}/progress.txt
	./build/cadical -q --compact=false -t ${timeout} ../instances/${year}/${file} >> ./output/${year}/${file}.output

done
