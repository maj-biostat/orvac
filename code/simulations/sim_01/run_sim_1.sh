#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_02.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_03.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_04.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_05.log -n 1000 -s 9676 -a 50 -d 0.75 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_06.log -n 1000 -s 1170 -a 50 -d 0.75 -b 0.5 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_07.log -n 1000 -s 3625 -a 50 -d 0.75 -b 0.3 -p 0.3 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_08.log -n 1000 -s 6147 -a 50 -d 0.75 -b 0.5 -p 0.5 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_09.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_10.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_11.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_12.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_13.log -n 1000 -s 7850 -a 30 -d 0.75 -b 0.3 -p 0.3 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_14.log -n 1000 -s 2414 -a 30 -d 0.75 -b 0.5 -p 0.5 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_15.log -n 1000 -s 6809 -a 30 -d 0.75 -b 0.3 -p 0.3 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "Null01" -l null_01_16.log -n 1000 -s 9951 -a 30 -d 0.75 -b 0.5 -p 0.5 -m 40 -t 40

## 25% increase in prob of conversio, 0 month increase in tte
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_02.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_03.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_04.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_05.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_06.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_07.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_08.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 30
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_09.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_10.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_11.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "SeroTrt01" -l sero_trt01_12.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 40

## 0% increase in prob of conversio, 5 month increase in tte
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_02.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_03.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_04.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_05.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_06.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_07.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_08.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_09.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_10.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_11.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt01" -l clin_trt01_12.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 45
#
## 0% increase in prob of conversio, 9 month increase in tte
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_02.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_03.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_04.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_05.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_06.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_07.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_08.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_09.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_10.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_11.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.3 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinTrt02" -l clin_trt02_12.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.5 -m 40 -t 49
#
## 25% increase in prob of conversio, 5 month increase in tte
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_02.log -n 1000 -s 1426 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_03.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_04.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_05.log -n 1000 -s 7560 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_06.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_07.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_08.log -n 1000 -s 2290 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_09.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_10.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_11.log -n 1000 -s 4565 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_12.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_13.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_14.log -n 1000 -s 6523 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_15.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_16.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_17.log -n 1000 -s 5906 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l clin_sero_trt01_18.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 45

# 45% increase in prob of conversio, 5 month increase in tte
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_02.log -n 1000 -s 1426 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_03.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_04.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_05.log -n 1000 -s 7560 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_06.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_07.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_08.log -n 1000 -s 2290 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_09.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_10.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_11.log -n 1000 -s 4565 -a 30 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_12.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_13.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.435 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_14.log -n 1000 -s 6523 -a 30 -d 0.50 -b 0.4 -p 0.580 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_15.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 35 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_16.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.435 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_17.log -n 1000 -s 5906 -a 30 -d 0.50 -b 0.4 -p 0.580 -m 40 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l clin_sero_trt02_18.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 40 -t 45

# 25% increase in prob of conversio, 9 month increase in tte
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_01.log -n 1000 -s 9676 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 39
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_02.log -n 1000 -s 1426 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_03.log -n 1000 -s 1170 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_04.log -n 1000 -s 2061 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_05.log -n 1000 -s 7560 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_06.log -n 1000 -s 1652 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_07.log -n 1000 -s 3625 -a 50 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_08.log -n 1000 -s 2290 -a 50 -d 0.50 -b 0.4 -p 0.500 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_09.log -n 1000 -s 6147 -a 50 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 49
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_10.log -n 1000 -s 7850 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_11.log -n 1000 -s 4565 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_12.log -n 1000 -s 2414 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_13.log -n 1000 -s 6037 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_14.log -n 1000 -s 6523 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_15.log -n 1000 -s 8620 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 35 -t 44
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_16.log -n 1000 -s 6809 -a 30 -d 0.50 -b 0.3 -p 0.375 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_17.log -n 1000 -s 5906 -a 30 -d 0.50 -b 0.4 -p 0.500 -m 40 -t 49
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l clin_sero_trt03_18.log -n 1000 -s 9951 -a 30 -d 0.50 -b 0.5 -p 0.625 -m 40 -t 49

# 45% increase in prob of conversio,  9 month increase in tte
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_01.log -n 1000 -s 6887 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_02.log -n 1000 -s 7513 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_03.log -n 1000 -s 8823 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_04.log -n 1000 -s 1348 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_05.log -n 1000 -s 9020 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_06.log -n 1000 -s 3435 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_07.log -n 1000 -s 5032 -a 50 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_08.log -n 1000 -s 5709 -a 50 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_09.log -n 1000 -s 2376 -a 50 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_10.log -n 1000 -s 2556 -a 30 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_11.log -n 1000 -s 7026 -a 30 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_12.log -n 1000 -s 2122 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 39
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_13.log -n 1000 -s 1030 -a 30 -d 0.75 -b 0.3 -p 0.435 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_14.log -n 1000 -s 4221 -a 30 -d 0.75 -b 0.4 -p 0.580 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_15.log -n 1000 -s 8250 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 46
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_16.log -n 1000 -s 2762 -a 30 -d 0.50 -b 0.3 -p 0.435 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_17.log -n 1000 -s 9919 -a 30 -d 0.50 -b 0.4 -p 0.580 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_18.log -n 1000 -s 8119 -a 30 -d 0.50 -b 0.5 -p 0.725 -m 30 -t 52
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_19.log -n 1000 -s 4221 -a 30 -d 0.75 -b 0.4 -p 0.580 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "ClinSeroTrt04" -l clin_sero_trt04_20.log -n 1000 -s 4221 -a 30 -d 0.75 -b 0.4 -p 0.580 -m 15 -t 50

