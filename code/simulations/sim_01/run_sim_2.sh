#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_01.log       -n 500 -s 9676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_01.log -n 500 -s 8676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_01.log -n 500 -s 7676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 39
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_01.log -n 500 -s 6676 -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 45 

