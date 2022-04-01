import sys
import pandas as pd

x = {'size':[],'range':[]}

for line in sys.stdin:
	row = line.strip().split()
	if 'WES100' not in line:
		continue
	s = int(row[2]) - int(row[1])
	x['size'].append(s)
	x['range'].append(row[0] + ':'+row[1]+'-'+row[2])
df = pd.DataFrame(x)

print(df.sort_values('size'))
