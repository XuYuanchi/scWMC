import magic
import pandas as pd


data_name = ["sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3"]
for j in range(1,11):
    for i in range(2,6):
        filename = "/home/suyanchi/project/scWMC/Clustering/data/" + data_name[i] + ".csv"
        X = pd.read_csv(filename, sep=',', header=0, index_col = 0)
        magic_operator = magic.MAGIC()
        X_magic = magic_operator.fit_transform(X.T)
        filename = "/home/suyanchi/project/scWMC/Clustering/MAGIC/" + str(j) + "/" + data_name[i] + ".csv"
        X_magic.to_csv(filename, sep = ",")

