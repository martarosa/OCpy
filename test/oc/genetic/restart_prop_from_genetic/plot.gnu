set terminal pngcairo enhanced dashed size 640,480  color lw 2 font 'Arial'

set xlabel 'time (fs)' font "Arial,18"
set ylabel 'Local Field (au)' font "Arial,18"
#set palette defined
set xrange[500*0.024188:550*0.024188]
#set yrange[-2.:2.]
#set ytics 1.0  
#set xtics 1.   
#set style line 1 lt 1 lc rgb "blue"  
#set style line 2 lt 1 lc rgb "red" 
#set style line 3 lt 1 lc rgb "black" 

set ylabel "Polarization field (au)"
set output "IEF-PCM.png"
p   "tmp/medium_t_1.dat"     u ($2*0.024188):5 w l ls 1 title "Test",\
    "out/medium_t_1.dat" u ($2*0.024188):5 w l dt 2 title "Reference"

