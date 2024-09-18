import csv
import pandas as pd
from collections import Counter
file1 = pd.read_csv(input())
x = file1["Iteration_Sum_ChiSquared_O_E_Difference"]
y = round(x,1)
d = Counter(y)
print(sorted(d.items()))
a = sorted(d.items())
with open('230426_CAP257_Decay_41week_Halflife_Rounded_Counts.csv', 'w') as f:
        writer = csv.writer(f)
        header = ['Sum', 'Count']
        writer.writerow(header)
        for value, count in a:
            output = [value, count]
            writer.writerow(output)

