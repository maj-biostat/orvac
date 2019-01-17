#!/bin/bash

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Test_null" -l test_null_01.log -n 100 -s 9676 -a 50 -b 0.4 -p 0.4 -m 30 -t 30
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Test_comb" -l test_comb_01.log -n 100 -s 4561 -a 50 -b 0.4 -p 0.6 -m 30 -t 35
