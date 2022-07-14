import magic
import pandas as pd

filename = "/home/suyanchi/project/scWMC/ICT/data/Deng_raw.csv"
X = pd.read_csv(filename, sep=",", header=0, index_col = 0)
magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(X.T)
filename = "/home/suyanchi/project/scWMC/ICT/Results/Deng/MAGIC.csv"
X_magic.to_csv(filename, sep = ",")