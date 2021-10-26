# -*- coding: utf-8 -*-

#How To Use:
# CL > python checkResult.py "year(or something directory name"

import csv
import sys
import os

args = sys.argv
year = args[1]

dirpath = "./output/"+year #outputファイル保存先Dir

#dirPathにあるファイルのファイル名をすべて取得
files = os.listdir(dirpath)
outputFiles = [f for f in files if os.path.isfile(os.path.join(dirpath, f))]

for outputFileName in outputFiles:
  path = dirpath+"/"+outputFileName
  base, ext = os.path.splitext(path);
  if ext != ".output": continue
  print(outputFileName)
  result = "UNKNOWN" #規定時間内に解けていない場合、target1がresultに入っていないためUNKNOWNを入れておく

	#各ファイルから必要なラインのみ取得	
  with open(path) as f:
    lines = f.readlines()
    for str in lines:
      target1 = "s "
      if str.startswith(target1):
        result = str.replace(target1,'').rstrip('\n') 

      target2 = "c total real time since initialization: "
      if str.startswith(target2):
        time = str.replace(target2,'').rstrip('seconds\n').replace(' ','')
  f.close()

  #取得した内容をCSVに書き出し
  csvPath = "./output/result.csv" #集計結果出力先
  with open(csvPath, "a") as f:
    writer = csv.writer(f)
    name = outputFileName.rstrip('.output')
    writer.writerow([year,name,result,time])
    #print(time+"\t"+result+"\t"+name)

  f.close()
