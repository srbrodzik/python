import pandas as pd
from datetime import datetime

times = [datetime(2019,7,31,12),datetime(2019,7,31,13)]
test_values = [1,2]
df = pd.DataFrame(index=times, columns=['test'])
df['test'] = test_values
print(df)

# new dataframe (in this example one of the rows repeats
new_times = [datetime(2019,7,31,13),datetime(2019,7,31,14),datetime(2019,7,31,15)]
new_values = [2,3,4]
df2 = pd.DataFrame(index=new_times, columns=['test'])
df2['test'] = new_values
print(df2)

# add new df to original df
df = pd.concat([df,df2])
print(df)

# drop duplicates
print('dropping duplicates')
df = df.drop_duplicates()

print(df)
print('yay')


