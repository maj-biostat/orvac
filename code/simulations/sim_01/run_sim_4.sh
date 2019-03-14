#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N04" -l i_01.log       -n 1000 -s 9441  -a 20 -d 1 -b 0.8 -p 0.80 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N05" -l i_01.log       -n 1000 -s 9422  -a 40 -d 1 -b 0.8 -p 0.80 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N06" -l i_01.log       -n 1000 -s 9482  -a 60 -d 1 -b 0.8 -p 0.80 -m 30 -t 30

# baseline 0.2, 30 months, no info delay, accrual from 20 to 60
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I17" -l i_01.log       -n 1000 -s 11917  -a 40 -d 1 -b 0.8 -p 0.90 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I18" -l i_01.log       -n 1000 -s 11918  -a 40 -d 1 -b 0.8 -p 0.90 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I19" -l i_01.log       -n 1000 -s 11919  -a 40 -d 1 -b 0.8 -p 0.90 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I20" -l i_01.log       -n 1000 -s 11920  -a 40 -d 1 -b 0.8 -p 0.90 -m 30 -t 50


