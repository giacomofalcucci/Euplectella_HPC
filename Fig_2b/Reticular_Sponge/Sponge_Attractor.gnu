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
   "Times-Roman" 24  fontscale 1.0 
set output 'Sponge_Attractor.eps'
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
set xtics  norangelimit autofreq 
set xtics 0.04
set xtics ("0" 0, "3.4" 0.04, "6.8" 0.08, "10.2" 0.12)
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq
set ytics 0.075
set ytics ("-12.8" -0.15, "-6.4" -0.075, "0" 0, "6.4" 0.075, "12.8" 0.15)
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit 500
set ztics ("100" 100, "500" 500, "1000" 1000, "1500" 1500, "2000" 2000)
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
set xlabel "u_x" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 0.00000 : 0.140000 ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "u_y" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [-0.15:0.15] noreverse writeback
set y2range [ * : * ] noreverse writeback
#set zlabel "Re" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [100:2000] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
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
## Last datafile plotted: "RUN_SINGLE-2000-T4000000/probe.dat"
set view 60,8
set ticslevel 0
set size 0.8,1.1
unset key
unset ztics
unset border
set border 15 # E' la SOMMA!!! guarda http://gnuplot.sourceforge.net/docs_4.2/node162.html
#splot 'RUN_SINGLE-100-T1000000/probe.dat' u 3:4:1 title 'Re =   100', 'RUN_SINGLE-500-T1500000/probe.dat' u 3:4:1 title 'Re =   500', 'RUN_SINGLE-1000-T2500000/probe.dat' u 3:4:1 title 'Re = 1000' lc 7, 'RUN_SINGLE-1500-T3000000/probe.dat' u 3:4:1 title 'Re = 1500' ps 0.75, 'RUN_SINGLE-2000-T4000000/probe.dat' u 3:4:1 title 'Re = 2000' w p pt 14 lc 3 ps 0.75
splot 'probe_g.0052_Re_100.dat' u 2:3:6 pt 15 title 'Re = 100', 'probe_g.0052_Re_500.dat' u 3:4:1 pt 13 title 'Re = 500', 'probe_g.0052_Re_1000.dat' u 2:3:6 title 'Re = 1000' lc 7 pt 7, 'probe_g.0052_Re_1500.dat' u 2:3:6 title 'Re = 1500' lc 6 pt 9, 'probe_g.0052_Re_2000.dat' u 2:3:6 lc 4 title 'Re = 2000' 
#    EOF
