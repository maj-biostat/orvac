#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01.log           -n 500 -s 9676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null02" -l null_02.log           -n 500 -s 9676 -a 30 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null03" -l null_03.log           -n 500 -s 9676 -a 50 -d 1.0 -b 0.3 -p 0.3 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null04" -l null_04.log           -n 500 -s 9676 -a 30 -d 1.0 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin01" -l clin_01.log          -n 1000 -s 8676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin02" -l clin_02.log          -n 1000 -s 7676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin03" -l clin_03.log          -n 1000 -s 6676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 45 
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin04" -l clin_04.log          -n 1000 -s 8676 -a 30 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin05" -l clin_05.log          -n 1000 -s 7676 -a 30 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Clin06" -l clin_06.log          -n 1000 -s 6676 -a 30 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 45 
