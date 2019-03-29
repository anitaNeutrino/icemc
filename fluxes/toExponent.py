## WebPlotDigitizer -> explicit icemc input format
## How to use: python toExponent.py yourFileNameWithoutExtension
import csv
import math
import os
import sys

new_rows = []
count=0
fileName=str(sys.argv[1])

with open(fileName+".dat", 'rb') as f, open(fileName+"C.dat", 'wb') as g:
    reader = csv.reader(f, delimiter=',')
    writer = csv.writer(g, delimiter=' ', skipinitialspace=True)
    g.write(str(sum(1 for line in open(fileName+".dat")))+"\n")
    for row in reader:
        count += 1
        row[0] = math.log10(float(row[0])*1e9) #Convert to eV and get log format
        row[1] = math.log10(float(row[1]))
        writer.writerow(row)
 
