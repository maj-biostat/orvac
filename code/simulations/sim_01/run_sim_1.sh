#!/bin/bash

# null case
#/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.3 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.4 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.5 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.6 -m 30 -t 30
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.3 -m 35 -t 35
#/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.4 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.5 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.6 -m 35 -t 35
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.3 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.4 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.5 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.6 -m 40 -t 40
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.3 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.4 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.5 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.6 -m 30 -t 30
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.3 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.4 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.5 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.6 -m 35 -t 35
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.3 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.4 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.5 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "Null01" -l log_null01_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.6 -m 40 -t 40


# sero effect
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.7 -m 30 -t 30
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.7 -m 35 -t 35
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.7 -m 40 -t 40
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.7 -m 30 -t 30
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.7 -m 30 -t 30
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.7 -m 35 -t 35
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.7 -m 35 -t 35
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.7 -m 40 -t 40
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "SeroTrt01" -l log_serotrt01_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.7 -m 40 -t 40


# clin effect 20% increase in median time to event
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.3 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.4 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.5 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.6 -m 30 -t 36
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.3 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.4 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.5 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.6 -m 35 -t 42
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.3 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.4 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.5 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.6 -m 40 -t 48
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.3 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.4 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.5 -m 30 -t 36
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.6 -m 30 -t 36
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.3 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.4 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.5 -m 35 -t 42
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.6 -m 35 -t 42
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.3 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.4 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.5 -m 40 -t 48
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt03" -l log_clintrt01_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.6 -m 40 -t 48

# clin effect 30% increase in median time to event
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.3 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.4 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.5 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.6 -m 30 -t 39
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.3 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.4 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.5 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.6 -m 35 -t 46
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.3 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.4 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.5 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.6 -m 40 -t 52
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.3 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.4 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.5 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.6 -m 30 -t 39
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.3 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.4 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.5 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.6 -m 35 -t 46
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.3 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.4 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.5 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinTrt02" -l log_clintrt01_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.6 -m 40 -t 52
#
# # clin effect 30% increase in median time to event trt sero blanket 0.7
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.7 -m 30 -t 39
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.7 -m 35 -t 46
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.7 -m 40 -t 52
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.7 -m 30 -t 39
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.7 -m 30 -t 39
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.7 -m 35 -t 46
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.7 -m 35 -t 46
#
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.7 -m 40 -t 52
# /usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt01" -l log_clintrt01_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.7 -m 40 -t 52



# clin effect 30% increase in median time to event trt sero blanket 0.8
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_01.log -n 1000 -s 9676 -a 50 -b 0.3 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_02.log -n 1000 -s 1426 -a 50 -b 0.4 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_03.log -n 1000 -s 1170 -a 50 -b 0.5 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_04.log -n 1000 -s 7012 -a 50 -b 0.6 -p 0.8 -m 30 -t 39

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_05.log -n 1000 -s 2061 -a 50 -b 0.3 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_06.log -n 1000 -s 7560 -a 50 -b 0.4 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_07.log -n 1000 -s 1652 -a 50 -b 0.5 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_08.log -n 1000 -s 8488 -a 50 -b 0.6 -p 0.8 -m 35 -t 46

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_09.log -n 1000 -s 3625 -a 50 -b 0.3 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_10.log -n 1000 -s 2290 -a 50 -b 0.4 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_11.log -n 1000 -s 6147 -a 50 -b 0.5 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_12.log -n 1000 -s 9272 -a 50 -b 0.6 -p 0.8 -m 40 -t 52

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_13.log -n 1000 -s 7850 -a 30 -b 0.3 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_14.log -n 1000 -s 4565 -a 30 -b 0.4 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_15.log -n 1000 -s 2414 -a 30 -b 0.5 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_16.log -n 1000 -s 6923 -a 30 -b 0.6 -p 0.8 -m 30 -t 39

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_17.log -n 1000 -s 6037 -a 30 -b 0.3 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_18.log -n 1000 -s 6523 -a 30 -b 0.4 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_19.log -n 1000 -s 8620 -a 30 -b 0.5 -p 0.8 -m 35 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_20.log -n 1000 -s 1820 -a 30 -b 0.6 -p 0.8 -m 35 -t 46

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_21.log -n 1000 -s 6809 -a 30 -b 0.3 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_22.log -n 1000 -s 5906 -a 30 -b 0.4 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_23.log -n 1000 -s 9951 -a 30 -b 0.5 -p 0.8 -m 40 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt02" -l log_clintrt02_24.log -n 1000 -s 7658 -a 30 -b 0.6 -p 0.8 -m 40 -t 52

# clin effect large increases in median time to event trt sero blanket 0.8
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_01.log -n 1000 -s 6887 -a 50 -b 0.3 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_02.log -n 1000 -s 7513 -a 50 -b 0.4 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_03.log -n 1000 -s 8823 -a 50 -b 0.5 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_04.log -n 1000 -s 5255 -a 50 -b 0.6 -p 0.8 -m 30 -t 39

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_05.log -n 1000 -s 1348 -a 50 -b 0.3 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_06.log -n 1000 -s 9020 -a 50 -b 0.4 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_07.log -n 1000 -s 3435 -a 50 -b 0.5 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_08.log -n 1000 -s 2802 -a 50 -b 0.6 -p 0.8 -m 30 -t 46

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_09.log -n 1000 -s 5032 -a 50 -b 0.3 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_10.log -n 1000 -s 5709 -a 50 -b 0.4 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_11.log -n 1000 -s 2376 -a 50 -b 0.5 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_12.log -n 1000 -s 3678 -a 50 -b 0.6 -p 0.8 -m 30 -t 52

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_13.log -n 1000 -s 2556 -a 30 -b 0.3 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_14.log -n 1000 -s 7026 -a 30 -b 0.4 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_15.log -n 1000 -s 2122 -a 30 -b 0.5 -p 0.8 -m 30 -t 39
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_16.log -n 1000 -s 5549 -a 30 -b 0.6 -p 0.8 -m 30 -t 39

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_17.log -n 1000 -s 1030 -a 30 -b 0.3 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_18.log -n 1000 -s 4221 -a 30 -b 0.4 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_19.log -n 1000 -s 8250 -a 30 -b 0.5 -p 0.8 -m 30 -t 46
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_20.log -n 1000 -s 2028 -a 30 -b 0.6 -p 0.8 -m 30 -t 46

/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_21.log -n 1000 -s 2762 -a 30 -b 0.3 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_22.log -n 1000 -s 9919 -a 30 -b 0.4 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_23.log -n 1000 -s 8119 -a 30 -b 0.5 -p 0.8 -m 30 -t 52
/usr/bin/Rscript main_2.R -f cfg1.yaml -o T -i "ClinSeroTrt03" -l log_clintrt03_24.log -n 1000 -s 5010 -a 30 -b 0.6 -p 0.8 -m 30 -t 52
