reset
#Check list of files
fgeos=system("ls geo*.dat")

print "file should be geo0.dat"
print fgeos


#from input file in mready
tscale=5e-5

do for [file in fgeos]{

reset
unset terminal 

name=substr(file,4,4)
geoplot=sprintf("NUKT%s.eps",name)
print geoplot

ftch=sprintf("time_changes%s.cvs",name)
print ftch

plot file every ::1 using ($3*tscale):6 w points t 'U'
hrange=(GPVAL_DATA_Y_MAX-GPVAL_DATA_Y_MIN)/2
xmax=GPVAL_DATA_X_MAX
xmin=GPVAL_DATA_X_MIN

plot file every ::1 using ($3*tscale):10 w points t 'T'
ave=GPVAL_DATA_Y_MIN+(GPVAL_DATA_Y_MAX-GPVAL_DATA_Y_MIN)/2

set terminal postscript eps size 4.0,4.0 enhanced color font 'Helvetica,15' linewidth 1 
set output geoplot 

#set xrange [0:xmax]
set xrange [xmin:xmax]
#set xrange [0:200000]

#set xtics 0.01



## set number of plots
NX=1; NY=4
DX=0.13; DY=0.08

SX=0.83/NX; SY=0.87/NY



set bmargin 0.6; set tmargin 0.6 ; set lmargin 0.; set rmargin 0.
## set the margin of each side of a plot as small as possible
## to connect each plot without space

set pointsize 0.2

set size 1,1
set format y "%2.2f"
set multiplot
set size SX,SY
#set xrange [:6000]
#set xrange [930*5:980*5]
##First Figure bottom
#set yrange[-0.3000:-0.2970]
set yrange [ave-hrange:ave+hrange]
#set autoscale y

#set yrange [-36.63947:-36.63941]


set origin DX,DY;

#in time
set xlabel "time / ps" 
#in steps 
#set xlabel "time / steps" 


set label 1 "T" at graph 0.05,0.8 center


#set ytics 0.1 

#in time
#plot fgeo using ($4/1000000):10 w points pt 5 t 'T'



#in steps
  
plot file using ($3*tscale):10 w l lc rgb "red" t ''



###Second Figure middle

set label 1 "K"

#set yrange [6:16]
#set yrange[6.5:10]
#set format y "%2.3f"

set autoscale y
set notitle
set origin DX,DY+SY;
set xtics format " " 
set xlabel " "
set ylabel 'Energy / E_h'

#reading from time
#plot fgeo using ($4/1000000):8 w points pt 5 t 'K'
#reading from steps
plot file using ($3*tscale):8 w points  lc rgb "red" t ''


## Third Figure-top


set label 1 "U" at graph 0.05,0.2 center


#set yrange [-36:-26]

#set title "100 O2; 100 H2; 20 H; 3000 K "

#set yrange[-35.5:-32]


set origin DX,DY+SY*2
set xlabel " "
set ylabel " "
#in time
#plot fgeo using ($4/1000000):6 w points pt 5 t 'U'
#in steps


plot file using ($3*tscale):6 w points lc rgb "red" t ''




unset label 1

set datafile separator ";"




set origin DX,DY+SY*3

set key left center 

set title"Temperature K / pressure atm"

#set xlabel "time / ps" 

#set format y ""
#set ytics(0,1,2)
set format y "%2.f"


#set key autotitle columnhead

set ylabel "{N}_{i}" rotate by 0
set logscale y
 

#set yrange [0.9:300]

#set label 'ox'   at graph 0.1, 0.6 textcolor rgb "red"     
set label 'hy'   at graph 0.1, 0.7 textcolor rgb "red"   
#set label 'H_2'  at graph 0.2, 0.6 textcolor rgb "blue"    
set label 'O_2'  at graph 0.2, 0.7 textcolor rgb "green" 
#set label 'HO'   at graph 0.3, 0.6 textcolor rgb "cyan" 
set label 'HO_2' at graph 0.3, 0.7 textcolor rgb "blue"
#set label 'H_2O' at graph 0.3, 0.4 textcolor rgb "black"

plot ftch u ($1*tscale):4  w l ls 1 lw 3 lc rgb "red"    t '' ,\
     ftch u ($1*tscale):7  w l ls 1 lw 3 lc rgb "green"  t '' ,\
     ftch u ($1*tscale):9  w l ls 1 lw 3 lc rgb "blue"   t '' # ,\
#     ftch u ($1*tscale):7  w l ls 1 lw 3 lc rgb "orange" t '' ,\
#     ftch u ($1*tscale):8  w l ls 1 lw 3 lc rgb "cyan"   t '' ,\
#     ftch u ($1*tscale):9  w l ls 1 lw 3 lc rgb "black"  t '' #,\


unset multiplot


}
