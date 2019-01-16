#!/bin/bash

#/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroSmallTrt" -l log_99.log -n 1 -s 9676 -a 50 -b 0.3 -p 0.3 -m 30 -t 35


# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_99.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.7 -m 30 -t 30


/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.0 -m 30 -t 30
