#!/bin/bash

# effects of accrual
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "A01" -l a_01.log       -n 500 -s 1  -a 20 -d 0.80 -b 0.3 -p 0.4 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "A02" -l a_02.log       -n 500 -s 2  -a 30 -d 0.80 -b 0.3 -p 0.4 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "A03" -l a_03.log       -n 500 -s 3  -a 40 -d 0.80 -b 0.3 -p 0.4 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "A04" -l a_04.log       -n 500 -s 4  -a 50 -d 0.80 -b 0.3 -p 0.4 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "A05" -l a_05.log       -n 500 -s 5  -a 60 -d 0.80 -b 0.3 -p 0.4 -m 30 -t 35



# null case
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01.log -n 500 -s 1  -a 50 -d 0.0 -b 0.2 -p 0.2 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null02" -l null_02.log -n 500 -s 2  -a 30 -d 0.0 -b 0.2 -p 0.2 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null03" -l null_04.log -n 500 -s 3  -a 30 -d 1.0 -b 0.2 -p 0.2 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null04" -l null_05.log -n 500 -s 4  -a 50 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null05" -l null_06.log -n 500 -s 5  -a 30 -d 0.0 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null06" -l null_08.log -n 500 -s 6  -a 30 -d 1.0 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null07" -l null_05.log -n 500 -s 7  -a 50 -d 0.0 -b 0.3 -p 0.3 -m 50 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null08" -l null_06.log -n 500 -s 8  -a 30 -d 0.0 -b 0.3 -p 0.3 -m 50 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null09" -l null_08.log -n 500 -s 9  -a 30 -d 1.0 -b 0.3 -p 0.3 -m 50 -t 50
#
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I01" -l i_01.log       -n 500 -s 10  -a 30 -d 1.0 -b 0.2 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I02" -l i_03.log       -n 500 -s 11  -a 30 -d 1.0 -b 0.2 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I03" -l i_01.log       -n 500 -s 12  -a 30 -d 1.0 -b 0.3 -p 0.4 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I04" -l i_03.log       -n 500 -s 13  -a 30 -d 1.0 -b 0.3 -p 0.6 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I05" -l i_01.log       -n 500 -s 14  -a 30 -d 1.0 -b 0.4 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I06" -l i_03.log       -n 500 -s 15  -a 30 -d 1.0 -b 0.4 -p 0.7 -m 30 -t 30
#
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC01" -l i_01.log       -n 500 -s 16  -a 30 -d 1.0 -b 0.2 -p 0.3 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC02" -l i_03.log       -n 500 -s 17  -a 30 -d 1.0 -b 0.2 -p 0.5 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC03" -l i_01.log       -n 500 -s 18  -a 30 -d 1.0 -b 0.3 -p 0.4 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC04" -l i_03.log       -n 500 -s 19  -a 30 -d 1.0 -b 0.3 -p 0.6 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC05" -l i_01.log       -n 500 -s 20  -a 30 -d 1.0 -b 0.4 -p 0.5 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC06" -l i_03.log       -n 500 -s 21  -a 30 -d 1.0 -b 0.4 -p 0.7 -m 30 -t 35
#
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC07" -l i_01.log       -n 500 -s 22  -a 30 -d 1.0 -b 0.2 -p 0.3 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC08" -l i_03.log       -n 500 -s 23  -a 30 -d 1.0 -b 0.2 -p 0.5 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC09" -l i_01.log       -n 500 -s 24  -a 30 -d 1.0 -b 0.3 -p 0.4 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC10" -l i_03.log       -n 500 -s 25  -a 30 -d 1.0 -b 0.3 -p 0.6 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC11" -l i_01.log       -n 500 -s 26  -a 30 -d 1.0 -b 0.4 -p 0.5 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC12" -l i_03.log       -n 500 -s 27  -a 30 -d 1.0 -b 0.4 -p 0.7 -m 30 -t 40
#
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC13" -l i_01.log       -n 500 -s 28  -a 30 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 32
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC14" -l i_03.log       -n 500 -s 29  -a 30 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC15" -l i_01.log       -n 500 -s 30  -a 30 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 37
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC16" -l i_03.log       -n 500 -s 31  -a 30 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "IC17" -l i_01.log       -n 500 -s 32  -a 30 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 41


