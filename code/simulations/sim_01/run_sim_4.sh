#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N04" -l i_01.log       -n 1000 -s 9401  -a 20 -d 1 -b 0.2 -p 0.20 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N05" -l i_01.log       -n 1000 -s 9402  -a 40 -d 1 -b 0.2 -p 0.20 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N06" -l i_01.log       -n 1000 -s 9402  -a 60 -d 1 -b 0.2 -p 0.20 -m 30 -t 30

# baseline 0.2, 30 months, no info delay, accrual from 20 to 60
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I13" -l i_01.log       -n 1000 -s 1913  -a 40 -d 0 -b 0.2 -p 0.25 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I14" -l i_01.log       -n 1000 -s 1914  -a 40 -d 0 -b 0.2 -p 0.25 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I15" -l i_01.log       -n 1000 -s 1915  -a 40 -d 0 -b 0.2 -p 0.25 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I16" -l i_01.log       -n 1000 -s 1916  -a 40 -d 0 -b 0.2 -p 0.25 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I17" -l i_01.log       -n 1000 -s 1917  -a 40 -d 0 -b 0.2 -p 0.30 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I18" -l i_01.log       -n 1000 -s 1918  -a 40 -d 0 -b 0.2 -p 0.30 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I19" -l i_01.log       -n 1000 -s 1919  -a 40 -d 0 -b 0.2 -p 0.30 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I20" -l i_01.log       -n 1000 -s 1920  -a 40 -d 0 -b 0.2 -p 0.30 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I21" -l i_01.log       -n 1000 -s 1881  -a 40 -d 0 -b 0.2 -p 0.40 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I22" -l i_01.log       -n 1000 -s 1842  -a 40 -d 0 -b 0.2 -p 0.40 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I23" -l i_01.log       -n 1000 -s 1883  -a 40 -d 0 -b 0.2 -p 0.40 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I24" -l i_01.log       -n 1000 -s 1824  -a 40 -d 0 -b 0.2 -p 0.40 -m 30 -t 50

# baseline 0.2, 30 months, 1 month info delay, accrual from 20 to 60
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D13" -l i_01.log       -n 1000 -s 6213  -a 40 -d 1.0 -b 0.2 -p 0.25 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D14" -l i_01.log       -n 1000 -s 6214  -a 40 -d 1.0 -b 0.2 -p 0.25 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D15" -l i_01.log       -n 1000 -s 6215  -a 40 -d 1.0 -b 0.2 -p 0.25 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D16" -l i_01.log       -n 1000 -s 6216  -a 40 -d 1.0 -b 0.2 -p 0.25 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D17" -l i_01.log       -n 1000 -s 6217  -a 40 -d 1.0 -b 0.2 -p 0.30 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D18" -l i_01.log       -n 1000 -s 6218  -a 40 -d 1.0 -b 0.2 -p 0.30 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D19" -l i_01.log       -n 1000 -s 6219  -a 40 -d 1.0 -b 0.2 -p 0.30 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D20" -l i_01.log       -n 1000 -s 6210  -a 40 -d 1.0 -b 0.2 -p 0.30 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D21" -l i_01.log       -n 1000 -s 6221  -a 40 -d 1.0 -b 0.2 -p 0.40 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D22" -l i_01.log       -n 1000 -s 6222  -a 40 -d 1.0 -b 0.2 -p 0.40 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D23" -l i_01.log       -n 1000 -s 6223  -a 40 -d 1.0 -b 0.2 -p 0.40 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D24" -l i_01.log       -n 1000 -s 6224  -a 40 -d 1.0 -b 0.2 -p 0.40 -m 30 -t 50

