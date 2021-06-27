import uproot3

base_dir = "pi_zero/pi0"
csv_name = "pi_zero.csv"

events= uproot3.open("{}/topo-cluster.pool.root".format(base_dir))

columns = events["ClusterTree"].keys()

fcolumns = []

for column in columns:
    column = str(column)
    fcolumn = column.split('\'')[1]
    fcolumns.append(fcolumn)


df = events["ClusterTree"].pandas.df(fcolumns)
df.to_csv("{}/{}".format(base_dir,csv_name))
