import pandas as pd

df = pd.read_csv('tvtst1.csv')
df1 = df.iloc[:, 1:3]
df1.to_csv('tvws_data.csv', index=False, header=False)
