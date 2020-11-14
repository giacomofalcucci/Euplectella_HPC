#!/opt/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 4    last modified 2018-06-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
set terminal postscript landscape enhanced defaultplex \
   leveldefault color colortext \
   dashlength 1.0 linewidth 1.0 pointscale 1.0 butt noclip \
   nobackground \
   palfuncparam 2000,0.003 \
   "Times-Roman" 23  fontscale 1.0 
set output 'C_d.eps'
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5 unsorted
set cntrparam firstlinetype 0
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit logscale autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit logscale autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "Re" offset 2
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 10.0000 : 10000.0 ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "C_D" offset 2
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 0.800000 : 14.000 ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
set logscale y 10
set logscale x 10
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "it_IT.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
GNUTERM = "aqua"
#set key top right maxrows 7 width 0.1
#set key box height 1
#set boxwidth -2
#set key box lt -1 lw 2
set key box vertical width -6 height 0.2 maxcols 1 spacing 1.05
set pointsize 1.5
set style line 1 lw 1.5
set ytics (0.8,1.0,2.0,4.0,6.0,8.0,10.0)
set grid
## Last datafile plotted: "C_d_Data_for_all_geos.dat"
#  rgb '#00ff5a'
#plot 'Wieselsberger_data_from_Cheng.csv' u 1:2 w p pt 4 lc 6 lw 2 title 'Wieselberger, 1922 - Cylinder', 'Tritton_data_from_Cheng.csv' u 1:2 w p pt 8 lc 4 lw 2.5 title 'Tritton, 1959 - Cylinder', 'Clift_data_from_Cheng.csv'  u 1:2 w p pt 12 lw 2 lc 1 title 'Clift, 1972 - Cylinder','Henderson_data_from_Sahin.csv' u 1:2 w p pt 14 lc 3 lw 2.5 title 'Henderson, 1995 - Cylinder', 'Posdziech_Grundmann_data_from_Shain.csv' u 1:2 w p pt 16 lw 2 lc 2 title 'Posdziech and Grundmann, 2001 - Cylinder', 'Shain_data.csv' u 1:2 w p pt 10 lc 7 lw 2.5 title 'Sahin and Owens, 2004 - Cylinder', 'Eq_1_from_Cheng.csv' u 1:2 w p pt 1 lc 16 ps 1.6 lw 2 title 'Cheng, 2013 - Cylinder', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 7 lc 7 title 'LBM - Cylinder', 'C_d_Data_for_all_geos.dat' u 1:3 w p pt 5 lc 4 title 'LBM - Cylinder with Ridges', 'C_d_Data_for_all_geos.dat' u 1:4 w p pt 13 lc rgb '#f800a7' title 'LBM - Reticular Sponge', 'C_d_Data_for_all_geos.dat' u 1:5 w p pt 15 lc 2 title 'LBM - {/Times-Italic Euplectella Aspergillum}', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 7 lc 7 title '', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 6 lc 8 lw 1 title ''
plot 'Kawamura_1986.csv' u 1:2 w p pt 8 lc rgb '#f800a7' lw 3 title 'Kawamura, 1986', 'Henderson_data_from_Sahin.csv' u 1:2 w p pt 14 lc 3 lw 3 title 'Henderson, 1995', 'Posdziech_Grundmann_data_from_Shain.csv' u 1:2 w p pt 16 lw 2 lc rgb '#404040' title 'Posdziech and Grundmann, 2001', 'Shain_data.csv' u 1:2 w p pt 10 lc 7 lw 3 title 'Sahin and Owens, 2004', 'Fujisawa_2005.csv' u 1:2 w p pt 4 lc rgb '#ff8400' lw 3 title 'Fujisawa, 2005', 'Hanchi_1999.csv' u 1:2 w p pt 12 lw 3 lc 1 title 'Hanchi, 1999', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 7 lc 7 title 'LBM - Cylinder', 'C_d_Data_for_all_geos.dat' u 1:3 w p pt 5 lc rgb '#1e00ff' title 'LBM - Cylinder with Ridges', 'C_d_Data_for_all_geos.dat' u 1:4 w p pt 13 lc 4 ps 1.6 title 'LBM - Reticular Sponge', 'C_d_Data_for_all_geos.dat' u 1:5 w p pt 15 lc rgb '#009019' ps 1.6 title 'LBM - {/Times-Italic Euplectella Aspergillum}', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 7 lc 7 title '', 'C_d_Data_for_all_geos.dat' u 1:2 w p pt 6 lc 8 lw 1 title ''
#    EOF


