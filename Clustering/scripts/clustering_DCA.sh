#!/bin/bash
data_name=("sc_CELseq2" "sc_10x" "sc_Droseq" "sc_10x_5cl" "sc_Celseq2_5cl_p1" "sc_Celseq2_5cl_p2" "sc_Celseq2_5cl_p3")
for((j=6;j<11;j++));do
    for((i=2;i<6;i++));do
        dca --nocheckcounts --nonorminput --nosizefactors /home/suyanchi/project/scWMC/Clustering/data/${data_name[i]}.csv /home/suyanchi/project/scWMC/Clustering/DCA/${j}/${data_name[i]}/
    done
done