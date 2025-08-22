import pandas as pd
from io import StringIO

# Data from the user's table, reconstructed manually
raw_data = '''ID T X Y Z R C
1 1 -1.5 -8.2 1.3 3.8 -1
2 3 0.8 4.2 1.4 2.3 1
3 3 21.1 41.1 1.4 1.3 2
4 3 89.2 115.1 -5.3 1.1 3
5 3 153.6 163.5 3.0 0.6 4
6 3 184.9 218.1 -2.3 0.5 5
7 3 220.9 238.7 -13.7 0.3 6
8 3 159.9 183.9 -23.6 1.1 4
9 3 211.9 270.9 -35.2 0.7 8
10 3 0.6 58.1 -2.6 1.1 2
11 3 -0.5 75.9 -18.6 1.1 10
12 3 -40.5 180.0 -55.3 1.0 11
13 3 -156.4 281.8 -65.9 0.8 12
14 3 0.3 186.8 -7.0 0.8 12
15 3 -51.5 235.7 -27.0 0.8 14
16 3 -85.2 297.6 -41.5 0.5 15
17 3 75.3 198.3 -54.6 0.8 10
18 3 130.6 228.6 -81.3 0.5 17'''

# Read data into DataFrame
sio = StringIO(raw_data)
df = pd.read_csv(sio, sep=' ')

# Output as space-separated CSV
df.to_csv('simpleSWCFile.swc', sep=' ', index=False)