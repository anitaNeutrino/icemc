######
# Code to access and download the 3FGL catalog and process
# the associated FITS file. The catalog is then converted to CSV.
# Optional: The light curves are then grabbed for each source
######
## fetch
import urllib2
## process fit
from astropy.io import fits
from astropy.table import Table
import csv
import json
import datetime
import sys
import os
import string
###### Get file and show download status
url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4yr_catalog/gll_psc_v16.fit"
file_name = url.split('/')[-1]
u = urllib2.urlopen(url)
f = open(file_name, 'wb')
meta = u.info()
file_size = int(meta.getheaders("Content-Length")[0])
cat_name = '3FGL catalog'
print "Downloading: %s. Bytes: %s" % (cat_name, file_size)

file_size_dl = 0
block_sz = 8192
while True:
    buffer = u.read(block_sz)
    if not buffer:
        break

    file_size_dl += len(buffer)
    f.write(buffer)
    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
    status = status + chr(8)*(len(status)+1)
    print status,

f.close()

###### Get data and print
data = fits.getdata('gll_psc_v16.fit',1)
table = Table(data)
print('Printing summary info about two sources. This should look reasonable.')
print(table['Source_Name','RAJ2000','DEJ2000','CLASS1','ASSOC1'][3:5])
print('')
###### Convert to csv
print('Writing catalogue to csv...')
with open('dataPre_3fgl.csv','w') as csvfile:       
    spamwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_ALL)
    for i in range(0,len(table)):
        # Select which variables we want for csv
        desVars = table['Source_Name','ASSOC1','RAJ2000','DEJ2000','CLASS1','TEVCAT_FLAG','ASSOC_TEV'][i]
        spamwriter.writerow(desVars)

## Removing white space
with open("dataPre_3fgl.csv") as f:
    reader = csv.reader(f, delimiter=",")
    with open("data_3fgl.csv", "w") as fo:
        writer = csv.writer(fo,quoting=csv.QUOTE_ALL)
        # add header
        writer.writerow(["Source Name","Association","RA","Dec","Classification (code)","TEV Cat", "TEV Association"])
        for rec in reader:
            writer.writerow(map(string.strip, rec))
os.remove("gll_psc_v16.fit")
print 'Done! Saved catalog to data_3fgl.csv'

##### Grabbing light curves for each object we just got the csv of
checker = True # Set to true if you want to grab light curves
fermi_start_date = datetime.datetime(2001,1,1)
time_offset = int(fermi_start_date.strftime("%s"))
count = 0

file = open('dataPre_3fgl.csv')
numline = len(file.readlines())
directory = './lightCurves'
if checker:
    if not os.path.exists(directory):
        os.makedirs(directory)
        print('Processing light curves for objects. This will take some time...')
    with open('dataPre_3fgl.csv') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        # For each source
        for row in reader:
            count = count + 1
            sys.stdout.write('\rOn object: %s of %s' % (count, numline))
            sys.stdout.flush()
            if count == 10: # remove if you want to get curves for all objects, this limits it to 10
                break # ^
            # One light curve file per source
            lcFileName = "./lightCurves/lightcurve_%s_%s.csv" % (row[2], row[3]) # Make sure you create a lightCurves dir
            with open(lcFileName,'w') as csvfile:
                writer = csv.writer(csvfile,delimiter=',')
                #print(row[2])
                #print(row[3])
                webpage = urllib2.urlopen("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/queryDB_Lightcurve.php?ra=%s&dec=%s" % (row[2],row[3]))
                #print("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/queryDB_Lightcurve.php?ra=%s&dec=%s" % (row[1],row[2]))
                js = json.loads(webpage.read()) 
                N = len(js)
                # For each data point 
                for i in range(N):
                    met = float(js[i]['time']) 
                    unixTime = met + time_offset ;
                    le_count = float(js[i]['nev'])
                    le_bg = float(js[i]['avnev'])
                    le_sig = float(js[i]['sigma'])
                    he_count = float(js[i]['he_nev'])
                    he_bg = float(js[i]['he_avnev'])
                    he_sig = float(js[i]['he_sigma'])
                    
                    writer.writerow([met,unixTime,le_count,le_bg,le_sig,he_count,he_bg,he_sig])
                    

os.remove("dataPre_3fgl.csv")
print('     All done!')
