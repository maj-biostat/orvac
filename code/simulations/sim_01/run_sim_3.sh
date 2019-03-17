#!/bin/bash

# null case
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N01" -l i_01.log       -n 500 -s 7404  -a 20 -d 0 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N02" -l i_01.log       -n 500 -s 7405  -a 40 -d 0 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N03" -l i_01.log       -n 500 -s 7406  -a 60 -d 0 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N04" -l i_01.log       -n 500 -s 7407  -a 20 -d 1 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N05" -l i_01.log       -n 500 -s 7708  -a 40 -d 1 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N06" -l i_01.log       -n 500 -s 7409  -a 60 -d 1 -b 0.3 -p 0.30 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N07" -l i_01.log       -n 500 -s 7410  -a 40 -d 0 -b 0.2 -p 0.20 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N08" -l i_01.log       -n 500 -s 7411  -a 40 -d 0 -b 0.4 -p 0.40 -m 30 -t 30
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N09" -l i_01.log       -n 500 -s 7412  -a 40 -d 0 -b 0.2 -p 0.20 -m 20 -t 20
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N10" -l i_01.log       -n 500 -s 7413  -a 40 -d 0 -b 0.2 -p 0.20 -m 20 -t 20
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N11" -l i_01.log       -n 500 -s 7414  -a 40 -d 0 -b 0.4 -p 0.40 -m 40 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "N12" -l i_01.log       -n 500 -s 7415  -a 40 -d 0 -b 0.4 -p 0.40 -m 40 -t 40

# baseline 0.3, 30 months, no info delay, accrual from 20 to 60
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I01" -l i_01.log       -n 500 -s 661001  -a 20 -d 0 -b 0.3 -p 0.35 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I02" -l i_01.log       -n 500 -s 661002  -a 20 -d 0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I03" -l i_01.log       -n 500 -s 661003  -a 20 -d 0 -b 0.3 -p 0.35 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I04" -l i_01.log       -n 500 -s 661004  -a 20 -d 0 -b 0.3 -p 0.35 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I05" -l i_01.log       -n 500 -s 661005  -a 20 -d 0 -b 0.3 -p 0.40 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I06" -l i_01.log       -n 500 -s 661006  -a 20 -d 0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I07" -l i_01.log       -n 500 -s 661007  -a 20 -d 0 -b 0.3 -p 0.40 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I08" -l i_01.log       -n 500 -s 661008  -a 20 -d 0 -b 0.3 -p 0.40 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I09" -l i_01.log       -n 500 -s 661009  -a 20 -d 0 -b 0.3 -p 0.50 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I10" -l i_01.log       -n 500 -s 661010  -a 20 -d 0 -b 0.3 -p 0.50 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I11" -l i_01.log       -n 500 -s 661011  -a 20 -d 0 -b 0.3 -p 0.50 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I12" -l i_01.log       -n 500 -s 661012  -a 20 -d 0 -b 0.3 -p 0.50 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I13" -l i_01.log       -n 500 -s 661013  -a 40 -d 0 -b 0.3 -p 0.35 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I14" -l i_01.log       -n 500 -s 661014  -a 40 -d 0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I15" -l i_01.log       -n 500 -s 661015  -a 40 -d 0 -b 0.3 -p 0.35 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I16" -l i_01.log       -n 500 -s 661016  -a 40 -d 0 -b 0.3 -p 0.35 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I17" -l i_01.log       -n 500 -s 661017  -a 40 -d 0 -b 0.3 -p 0.40 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I18" -l i_01.log       -n 500 -s 661018  -a 40 -d 0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I19" -l i_01.log       -n 500 -s 661019  -a 40 -d 0 -b 0.3 -p 0.40 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I20" -l i_01.log       -n 500 -s 661020  -a 40 -d 0 -b 0.3 -p 0.40 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I21" -l i_01.log       -n 500 -s 661021  -a 40 -d 0 -b 0.3 -p 0.50 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I22" -l i_01.log       -n 500 -s 661022  -a 40 -d 0 -b 0.3 -p 0.50 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I23" -l i_01.log       -n 500 -s 661023  -a 40 -d 0 -b 0.3 -p 0.50 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I24" -l i_01.log       -n 500 -s 661024  -a 40 -d 0 -b 0.3 -p 0.50 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I25" -l i_01.log       -n 500 -s 661025  -a 60 -d 0 -b 0.3 -p 0.35 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I26" -l i_01.log       -n 500 -s 661026  -a 60 -d 0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I27" -l i_01.log       -n 500 -s 661027  -a 60 -d 0 -b 0.3 -p 0.35 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I28" -l i_01.log       -n 500 -s 661028  -a 60 -d 0 -b 0.3 -p 0.35 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I29" -l i_01.log       -n 500 -s 661029  -a 60 -d 0 -b 0.3 -p 0.40 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I30" -l i_01.log       -n 500 -s 661030  -a 60 -d 0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I31" -l i_01.log       -n 500 -s 661031  -a 60 -d 0 -b 0.3 -p 0.40 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I32" -l i_01.log       -n 500 -s 661032  -a 60 -d 0 -b 0.3 -p 0.40 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I33" -l i_01.log       -n 500 -s 661033  -a 60 -d 0 -b 0.3 -p 0.50 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I34" -l i_01.log       -n 500 -s 661034  -a 60 -d 0 -b 0.3 -p 0.50 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I35" -l i_01.log       -n 500 -s 661035  -a 60 -d 0 -b 0.3 -p 0.50 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "I36" -l i_01.log       -n 500 -s 661036  -a 60 -d 0 -b 0.3 -p 0.50 -m 30 -t 50


# baseline 0.3, 30 months, 1 month info delay, accrual from 20 to 60
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D01" -l i_01.log       -n 500 -s 2001  -a 20 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D02" -l i_01.log       -n 500 -s 2002  -a 20 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D03" -l i_01.log       -n 500 -s 2003  -a 20 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D04" -l i_01.log       -n 500 -s 2004  -a 20 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D05" -l i_01.log       -n 500 -s 2005  -a 20 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D06" -l i_01.log       -n 500 -s 2006  -a 20 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D07" -l i_01.log       -n 500 -s 2007  -a 20 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D08" -l i_01.log       -n 500 -s 2008  -a 20 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D09" -l i_01.log       -n 500 -s 2009  -a 20 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D10" -l i_01.log       -n 500 -s 2000  -a 20 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D11" -l i_01.log       -n 500 -s 2011  -a 20 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D12" -l i_01.log       -n 500 -s 2012  -a 20 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D13" -l i_01.log       -n 500 -s 2013  -a 40 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D14" -l i_01.log       -n 500 -s 2014  -a 40 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D15" -l i_01.log       -n 500 -s 2015  -a 40 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D16" -l i_01.log       -n 500 -s 2016  -a 40 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 50
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D17" -l i_01.log       -n 500 -s 2017  -a 40 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D18" -l i_01.log       -n 500 -s 2018  -a 40 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D19" -l i_01.log       -n 500 -s 2019  -a 40 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 45
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D20" -l i_01.log       -n 500 -s 2010  -a 40 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D21" -l i_01.log       -n 500 -s 2021  -a 40 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 35
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D22" -l i_01.log       -n 500 -s 2022  -a 40 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 40
/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D23" -l i_01.log       -n 500 -s 2023  -a 40 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D24" -l i_01.log       -n 500 -s 2024  -a 40 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D25" -l i_01.log       -n 500 -s 2025  -a 60 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D26" -l i_01.log       -n 500 -s 2026  -a 60 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D27" -l i_01.log       -n 500 -s 2027  -a 60 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D28" -l i_01.log       -n 500 -s 2028  -a 60 -d 1.0 -b 0.3 -p 0.35 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D29" -l i_01.log       -n 500 -s 2029  -a 60 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D30" -l i_01.log       -n 500 -s 2020  -a 60 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D31" -l i_01.log       -n 500 -s 2031  -a 60 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D32" -l i_01.log       -n 500 -s 2032  -a 60 -d 1.0 -b 0.3 -p 0.40 -m 30 -t 50
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D33" -l i_01.log       -n 500 -s 2033  -a 60 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 35
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D34" -l i_01.log       -n 500 -s 2034  -a 60 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 40
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D35" -l i_01.log       -n 500 -s 2035  -a 60 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 45
#/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i "D36" -l i_01.log       -n 500 -s 2036  -a 60 -d 1.0 -b 0.3 -p 0.50 -m 30 -t 50

