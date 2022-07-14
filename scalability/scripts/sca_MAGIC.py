import magic
import pandas as pd

filename = "/home/grads/fuzhowang2/suyanchi/scWMC/scalability/data/100k.csv"
X = pd.read_csv(filename, index_col = 0)
magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(X.T)
filename = "/home/grads/fuzhowang2/suyanchi/scWMC/scalability/Results/MAGIC.csv"
X_magic.to_csv(filename, sep = ",")