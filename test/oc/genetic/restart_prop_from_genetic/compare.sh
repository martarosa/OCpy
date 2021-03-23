#!/bin/bash

#gnuplot plot.gnu

threshold0=0.000000001
threshold1=0.003

diff0_medium=$(awk 'BEGIN{d2=0}
            FNR==NR&&FNR>1&&$2>120{a[FNR-1]=$5;next} 
	    FNR!=NR&&FNR>1&&$2>120{d2+=($5-a[FNR-1])^2}
	    END{print d2/(FNR-1)}' out/medium_t_1.dat medium_t_1.dat)
#echo $diff0_medium

diff1_medium=$(awk 'BEGIN{d2=0}
            $2>120{d2+=($4-$6)^2} 
	    END{print sqrt(d2/(NR-1))}' medium_t_1.dat)


