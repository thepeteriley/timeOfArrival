# timeOfArrival
Code and data for CCMC Scoreboard Time of arrival forecasts analysis

The following is taken from the header information of the .R script....

# Driver routine to read in the results of the CCMC shock arrival times forecasts, perform some statistics on them, 
# and make the figures found in the paper: 
#
# "Forecasting the Arrival Time of Coronal Mass Ejections: Analysis of the CCMC CME Scoreboard"
#
# to be, or already published in "Space Weather". 
#
# by:
#
# Pete Riley\affil{1}, Leila Mays\affil{2}, Jesse Andries\affil{3}, Tanja Amerstorfer\affil{4},  
# Douglas Biesecker\affil{5}, V\'eronique Delouille\affil{3}, Mateja Dumbovi\'{c}\affil{6,7}, Xueshang Feng\affil{8}, Edmund Henley\affil{9}, Jon A. Linker\affil{1}, Christian M\"{o}stl\affil{4}, Marlon Nu\~{n}ez\affil{10},Vic Pizzo\affil{5}, Manuela Temmer\affil{4}, W.K. Tobiska\affil{11}, C. Verbeke\affil{12}, Matthew J West\affil{3}, and Xinhua Zhao\affil{6}
}
# \affiliation{1}{Predictive Science Inc., San Diego, USA}
# \affiliation{2}{NASA/GSFC, Greenbelt, MD 20771, USA}
# \affiliation{3}{Solar-Terrestrial Center of Excellence, Royal Observatory of Belgium, Ringlaan 3, B-1180 Brussels, Belgium}
# \affiliation{4}{Space Research Institute, Austrian Academy of Sciences, 8042 Graz, Austria}
# \affiliation{5}{Space Weather Prediction Center, NOAA, Boulder, Colorado, USA}
# \affiliation{6}{Institute of Physics, University of Graz, Graz, Austria}
# \affiliation{7}{Hvar Observatory, Faculty of Geodesy, University of Zagreb, Zagreb, Croatia}
# \affiliation{8}{SIGMA Weather Group, State Key Laboratory of Space Weather, National Space Science Center, Chinese Academy of Sciences, Beijing 100190, China}
# \affiliation{9}{Met Office, FitzRoy Road, Exeter, Devon, UK}
# \affiliation{10}{Department of Languages and Computer Sciences, Universidad de M\'{a}laga, Málaga, Spain}
# \affiliation{11}{Space Environment Technologies, Pacific Palisades, CA 90272, USA}
# \affiliation{12}{Centre for Mathematical Plasma-Astrophysics, KU Leuven, Leuven, Belgium}
#
# Written 09/21/17 by PR. Continuously updated through July 18, 2018, at which point it was frozen for publication. 
#
# This is intended to provide support for the analysis presented in the aforementioned paper; however,
# you are free to take/modify this code for their own purpose(s). 
# 
# Disclaimer: 
#
# It includes all the gory details and transient changes so that I, or someone else 
# can recover all/earlier results. It's not supposed to be "production" code. 
#
# Please contact Pete Riley (pete@predsci.com) with any issues or suggested changes - 
# I'm more than happy to help set this up if you want to test YOUR forecasts against 
# those in the CCMC Scoreboard. 
#
# At some point, I'll update this to be a general purpose tool for testing/assessing other 
# forecasts to try to improve on the ones that are currently in the CCMC. 
#
# USAGE: 
#
# Point the "file" variable to the location of the csv file (included with the GitHub repository).
#
# choose one of the plot/analysis options (myPlot)
#
# Dependencies: Uses some R packages that are not installed by default. 
# You'll be prompted to install (via error messages) if you don't have them. 
