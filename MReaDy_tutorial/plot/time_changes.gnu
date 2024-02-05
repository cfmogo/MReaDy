reset

fgeos=system("ls time_changes*.cvs")

print "file should be time_changes0.cvs"
print fgeos


#from input file in mready
tscale=5e-5

do for [file in fgeos]{

reset
unset terminal

name=substr(file,13,13)
geoplot=sprintf("time_changes%s.eps",name)
print geoplot

ftch=sprintf("time_changes%s.cvs",name)
print ftch


set terminal postscript eps size 4.0,4.0 enhanced color font 'Helvetica,15' linewidth 1 
set output geoplot 


set datafile separator ";"


set key left center 

set title"Reading time\\\_changes.cvs"

set xlabel "time / ps" 

set format y "%2.f"

set key autotitle columnhead

set ylabel "{N}_{i}" rotate by 0
set logscale y
 

set yrange [0.9:150]
 
#set label 'ox'       at graph 0.1, 0.8 textcolor rgb "red"     
#set label 'hy'       at graph 0.1, 0.6 textcolor rgb "green"   
#set label 'H2'       at graph 0.2, 0.8 textcolor rgb "blue"    
#set label 'O2'       at graph 0.2, 0.6 textcolor rgb "magenta" 
#set label 'HO'       at graph 0.3, 0.8 textcolor rgb "orange"  
#set label 'HO2'      at graph 0.3, 0.6 textcolor rgb "black"   
#set label 'H2O2(1A)' at graph 0.4, 0.8 textcolor rgb "grey"  


plot for [i=3:22] ftch u ($1*tscale):i w l# ls 1 lw 3 lc rgb "red"    
#     ftch u ($1*tscale):6 w l ls 1 lw 3 lc rgb "green"     ,\
#     ftch u ($1*tscale):7 w l ls 1 lw 3 lc rgb "blue"      ,\
#     ftch u ($1*tscale):9 w l ls 1 lw 3 lc rgb "magenta"   ,\
#     ftch u ($1*tscale):8 w l ls 1 lw 3 lc rgb "orange"    ,\
#     ftch u ($1*tscale):10 w l ls 1 lw 3 lc rgb "black"    ,\
#     ftch u ($1*tscale):15 w l ls 2 lw 3 lc rgb "grey"  
}
