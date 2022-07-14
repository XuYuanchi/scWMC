import magic
import pandas as pd

filename = "/home/grads/fuzhowang2/suyanchi/scWMC/DEG/data/DEG.raw.tsv"
X = pd.read_csv(filename, sep="\t", header=0, index_col = 0)
magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(X.T)
filename = "/home/grads/fuzhowang2/suyanchi/scWMC/DEG/Results/MAGIC.csv"
X_magic.to_csv(filename, sep = ",")