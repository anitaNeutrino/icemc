######
# Code to access and download the 3FGL catalog and process
# the associated FITS file. The catalog is then converted to CSV.
######
## fetch
import urllib2
## process fit
from astropy.io import fits
from astropy.table import Table
import csv
## root
###### Get file and show download status
url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4yr_catalog/gll_psc_v16.fit"
file_name = url.split('/')[-1]
u = urllib2.urlopen(url)
f = open(file_name, 'wb')
meta = u.info()
file_size = int(meta.getheaders("Content-Length")[0])
print "Downloading: %s Bytes: %s" % (file_name, file_size)

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
print(table['Source_Name','DEJ2000','RAJ2000','CLASS1','ASSOC1'][3:5])
print('')
###### Convert to csv
print('Writing catalogue to csv...')
with open('data_3fgl.csv','w') as csvfile:       
    spamwriter = csv.writer(csvfile, delimiter=',')
    for i in range(0,len(table)):
        # Select which variables we want
        spamwriter.writerow(table['Source_Name','DEJ2000','RAJ2000','GLON','GLAT','CLASS1','ASSOC1'][i])
print 'Done! Saved catalog to data_3fgl.csv'
