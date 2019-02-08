#----
# Program to pull supernovae info from the Open Supernova Catalog.
# Uses variable-length, nested tables, so some default values are used.
# This is rather annoying to pass as many dates / co-ord systems have completely different formats!
# Some data is also input incorrectly, so there are adjustments for that...
#---
# Default string: Unknown -> now just ""
# Default number: -999
#---

import urllib2
import csv
import json
import time
import calendar
import sys
import os
import string
import re

print('Scanning catalog...')
snFileName = "supernovae.csv"
with open(snFileName,'w') as csvfile:
    writer = csv.writer(csvfile,delimiter='%',quoting=csv.QUOTE_MINIMAL)
    webpage = urllib2.urlopen("https://raw.githubusercontent.com/astrocatalogs/supernovae/master/output/catalog.json")
    js = json.loads(webpage.read()) 
    N = len(js)
    print('Processing entries...')
    # For each data point 
    for i in range(N):
        # need the encode part or it won't work:
        # i.e. 'ascii' codec can't encode character u'\u2011'
        name = js[i]['name'].encode('ascii', 'ignore').decode('ascii')
        # Nested -> pick first entry for alias
        ### We have a nested variable length table, need to check if key exists before assigning it a default value
        if("alias" in js[i]):
            alias = js[i]['alias'][0]['value'].encode('ascii', 'ignore').decode('ascii') # Only choose the first field for now
        else:
            alias = ""
            
        if("discoverer" in js[i]):
            discoverer = js[i]['discoverer'][0]['value'].encode('ascii', 'ignore').decode('ascii')
        else:
            discoverer = ""

        ## Checking for proper date format
        if("discoverdate" in js[i]):
            discoverdate = js[i]['discoverdate'][0]['value']
            # 01/01/1000 format
            if re.compile('.*/.*/.*').match(discoverdate) is not None:
                # One entry is using MM/DD/YYYY format... account for this
                discoveryearMismatch = discoverdate.split('/')[2]
                discoverdayMismatch = discoverdate.split('/')[1]
                discovermonthMismatch = discoverdate.split('/')[0]
                if(int(discoveryearMismatch) > 31):
                    discoverdate = discoveryearMismatch + "/" + discovermonthMismatch + "/" + discoverdayMismatch
                    # 01/2001 format
            elif re.compile('.*/.*').match(discoverdate) is not None:
                # set to middle of month
                discoverdate = discoverdate + '/15'
                # 2001 format
            elif re.compile('.*').match(discoverdate) is not None:
                # set to middle of year
                discoverdate = discoverdate + '/07/02'
                # if none of the above formats, just set to the default
            else:
                discoverdate = "-999"
        else:
            discoverdate = "-999"

        # Some dates are super old and predate unix time. Reject those, or those with erroneous dates.
        if(discoverdate == "-999"):
            discoverUnixTime = -999

        else:
            discoveryear = discoverdate.split('/')[0]
            if(int(discoveryear) < 1971 or int(discoveryear) > 2020):
                discoverUnixTime = -999
            else:
                discoverUnixTime = calendar.timegm(time.strptime(discoverdate,'%Y/%m/%d'))
                True
                # end
        
        ## Checking for proper date format
        if("maxdate" in js[i]):
            maxdate = js[i]['maxdate'][0]['value']
            # 01/01/1000 format
            if re.compile('.*/.*/.*').match(maxdate) is not None:
                # One entry is using MM/DD/YYYY format... account for this
                maxyearMismatch = maxdate.split('/')[2]
                maxdayMismatch = maxdate.split('/')[1]
                maxmonthMismatch = maxdate.split('/')[0]
                if(int(maxyearMismatch) > 31):
                    maxdate = maxyearMismatch + "/" + maxmonthMismatch + "/" + maxdayMismatch
                    # 01/2001 format
            elif re.compile('.*/.*').match(maxdate) is not None:
                # set to middle of month
                maxdate = maxdate + '/15'
                # 2001 format
            elif re.compile('.*').match(maxdate) is not None:
                # set to middle of year
                maxdate = maxdate + '/07/02'
                # if none of the above formats, just set to the default
            else:
                maxdate = "-999"
        else:
            maxdate = "-999"

        # Some dates are super old and predate unix time. Reject those.
        if(maxdate == "-999"):
            maxUnixTime = -999

        else:
            maxyear = maxdate.split('/')[0]
            if(int(maxyear) < 1971):
                maxUnixTime = -999
            else:
                maxUnixTime = calendar.timegm(time.strptime(maxdate,'%Y/%m/%d'))
                True
                # end
                
        if("maxappmag" in js[i]):
            maxappmag = float(js[i]['maxappmag'][0]['value'])
        else:
            maxdappmag = -999
            
        if("host" in js[i]):
            host = js[i]['host'][0]['value']
        else:
            host = ""

        if("ra" in js[i]):
            ra = js[i]['ra'][0]['value']
            # get ra in better format
            if re.compile('.*:.*:.*').match(ra) is not None:
                True
            elif re.compile('.*:.').match(ra) is not None:
                ra = ra + ':00'
            else:
                ra = ra + ':00' + ':00'
                
            raH = float(ra.split(':')[0])
            raM = float(ra.split(':')[1])
            raS = float(ra.split(':')[2])
            ra = 15*(raH + raM/60 + raS/3600)
            ra = "{0:.4f}".format(ra) # limit precision to 4 dp, it gets ugly otherwise
        else:
            ra = -999

        if("dec" in js[i]):
            dec = js[i]['dec'][0]['value']
            # get dec in better format
            if re.compile('.*:.*:.*').match(dec) is not None:
                True
            elif re.compile('.*:.').match(dec) is not None:
                dec = dec + ':00'
            else:
                dec = dec + ':00' + ':00'

            # We must be very careful about the minus sign for declination
            decH = dec.split(':')[0]
            decM = float(dec.split(':')[1])
            decS = float(dec.split(':')[2])
            if '-' in decH:
                dec = -1*(-float(decH) + decM/60 + decS/3600)
            else:
                dec = float(decH) + decM/60 + decS/3600
            dec = "{0:.4f}".format(dec)
        else:
            dec = -999

        if("redshift" in js[i]):
            redshift = float(js[i]['redshift'][0]['value'])
        else:
            redshift = -999

        if("claimedtype" in js[i]):
            claimedtype = js[i]['claimedtype'][0]['value']
        else:
            claimedtype = ""
            
        writer.writerow([name,alias,discoverer,discoverUnixTime,maxUnixTime,maxappmag,host,ra,dec,redshift,claimedtype])

print('Written to catalog %s') % snFileName
