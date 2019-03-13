#----
#1.
# Program to pull supernovae info from the Open Supernova Catalog.
# Uses variable-length, nested tables, so some default values are used.
# This is rather annoying to pass as many dates / co-ord systems have completely different formats!
# Some data is also input incorrectly, so there are adjustments for that...
#---
# Default string: Unknown -> now just ""
# Default number: -999
#---
#2. Automatically pull TNS tsv

import urllib2
import csv
import json
import time
import calendar
import sys
import os
import string
import re

#Initials
dataDir = './data'
if not os.path.exists(dataDir):
    os.makedirs(dataDir)

#1.
print('Scanning Open Source Catalog...')
snFileName = "./data/supernovaeOSC.tsv"
with open(snFileName,'w') as csvfile:
    writer = csv.writer(csvfile,delimiter='\t',quoting=csv.QUOTE_MINIMAL)
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

print('Catalog written to: %s') % snFileName
print('------------------------------------')
print('')

#2.
### Same for TNS
print('Scanning Transient Name Server...')
tnsFileName = "./data/supernovaeTNS.tsv"
### CSV already made via URL manipulation # this one is for classified SNe between 01/06/2016 - 01/06/2017
### Can change variables wanted, and their range by editing the below URL (limited to 1000 lines...):
# VERY IMPORTANT NOTE: Some info may be missing for the older data (the TNS only became the official SNe reporting method in 2016). If you want other flights (their associated times will be 00:00:00), check more thoroughly...
webpage0 = urllib2.urlopen("https://wis-tns.weizmann.ac.il/search?&name=&name_like=0&isTNS_AT=all&public=all&unclassified_at=0&classified_sne=1&ra=&decl=&radius=&coords_unit=arcsec&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=2016-06-01&date_end%5Bdate%5D=2017-06-01&discovery_mag_min=&discovery_mag_max=&internal_name=&redshift_min=&redshift_max=&spectra_count=&discoverer=&classifier=&discovery_instrument%5B%5D=&classification_instrument%5B%5D=&hostname=&associated_groups%5B%5D=null&ext_catid=&num_page=1000&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&format=tsv")
webContent0 = webpage0.read();
f = open(tnsFileName,'w')
f.write(webContent0)
#this one is for classified SNe between 01/06/2015 - 31/05/2016 
webpage1 = urllib2.urlopen("https://wis-tns.weizmann.ac.il/search?&name=&name_like=0&isTNS_AT=all&public=all&unclassified_at=0&classified_sne=1&ra=&decl=&radius=&coords_unit=arcsec&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=2015-06-01&date_end%5Bdate%5D=2016-05-31&discovery_mag_min=&discovery_mag_max=&internal_name=&redshift_min=&redshift_max=&spectra_count=&discoverer=&classifier=&discovery_instrument%5B%5D=&classification_instrument%5B%5D=&hostname=&associated_groups%5B%5D=null&ext_catid=&num_page=1000&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&format=tsv")
webContent1 = webpage1.read()
f.write(webContent1)
# (this one is for classified SNe between 01/06/2014 - 31/05/2015)
webpage2 = urllib2.urlopen("https://wis-tns.weizmann.ac.il/search?&name=&name_like=0&isTNS_AT=all&public=all&unclassified_at=0&classified_sne=1&ra=&decl=&radius=&coords_unit=arcsec&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=2014-06-01&date_end%5Bdate%5D=2015-05-31&discovery_mag_min=&discovery_mag_max=&internal_name=&redshift_min=&redshift_max=&spectra_count=&discoverer=&classifier=&discovery_instrument%5B%5D=&classification_instrument%5B%5D=&hostname=&associated_groups%5B%5D=null&ext_catid=&num_page=1000&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&format=tsv")
webContent2 = webpage2.read()
f.write(webContent2)
f.close

print('Catalog written to: %s') % tnsFileName
print('------------------------------------')
print('Finished processing all SNe catalogs!')
