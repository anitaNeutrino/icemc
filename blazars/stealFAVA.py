#! /usr/bin/env python

# Script to grab fava events
# Since I'm downloading from a web site and need to parse json, it was easier to just do it in python. 
# You need PyROOT obviously to run this... 



from ROOT import TFile,TTree, gROOT, TObjString
import urllib 
import json 
import datetime
import csv
import sys 

## Limit this if you want a lower run time / more time constrained data
week_start = 330; 
week_end = 340; 


gROOT.ProcessLine("#include \"fava.h\"") 
from ROOT import FAVAEntry; 

fava = FAVAEntry() 


f = TFile("fava.root","recreate"); 

t = TTree("fava","FAVA tree"); 

t.Branch('fava', fava); 

# put some in by hand, 
# for RGB J2243+203, VERITAS sets an upper limit of ~1 and optical sets a lower limit of 0.4, so let's call it 0.7
# TXS 0506+056 is measured but not in the NED catalog 
# PKS 1441+25 has a crazy preferred value in NED... 

zs = {"RGB J2243+203":0.7, "TXS 0506+056":0.3365, "PKS 1441+25":0.939} 

fermi_start_date = datetime.datetime(2001,1,1)
time_offset = int(fermi_start_date.strftime("%s"))



#catalog = open("3fgl.csv")
catalog = open("data_3fgl.csv")
reader = csv.DictReader(catalog)
## Make reading whitespaced data easier
#reader = (
#    dict((k, v.strip()) for k, v in row.items() if v) for row in reader)

source_dict = {}

for row in reader: 
 source_dict[row["Source Name"]] = (row["Classification (code)"].lower(),row["TEV Cat"])




for week in range(week_start,week_end+1): 
  url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/queryDB_2FAV.php?typeOfRequest=SourceList&week=%d&threshold=6Sigma" % week
  fava.week = week 
  page = urllib.urlopen(url); 
  js = json.loads(page.read()) 

  N = len(js) 
  for i in range(N): 

    fava.num = int(float(js[i]['num']))
    fava.ra = float(js[i]['best_ra'])
    fava.dec = float(js[i]['best_dec'])
    fava.r95 = float(js[i]['best_r95'])
    fava.met_tmin = int(float(js[i]['tmin']))
    fava.met_tmax = int(float(js[i]['tmax']))
    fava.unix_tmin = fava.met_tmin + time_offset
    fava.unix_tmax = fava.met_tmax + time_offset
    fava.le_nev =  int(float(js[i]['nev']))
    fava.he_nev =  int(float(js[i]['he_nev']))
    fava.le_avnev =  float(js[i]['avnev'])
    fava.he_avnev =  float(js[i]['he_avnev'])
    fava.le_sigma =  float(js[i]['sigma'])
    fava.he_sigma =  float(js[i]['he_sigma'])
    fava.le_ts =  float(js[i]['le_ts'])
    fava.he_ts =  float(js[i]['he_ts'])
    fava.le_ts_sigma =  float(js[i]['le_tssigma'])
    fava.he_ts_sigma =  float(js[i]['he_tssigma'])
    fava.le_flux =  float(js[i]['le_flux'])
    fava.he_flux =  float(js[i]['he_flux'])
    fava.le_flux_err =  float(js[i]['le_fuxerr'])
    fava.he_flux_err =  float(js[i]['he_fuxerr'])
    fava.le_index =  float(js[i]['le_index'])
    fava.he_index =  float(js[i]['he_index'])
    assoc =  js[i]['assoc']
    fava.association = TObjString(assoc)
    glassoc = js[i]['fglassoc'] 
    fava.association_3fgl = TObjString(glassoc)
    fava.tevcat = False 
    if glassoc in source_dict: 
      fava.source_class = TObjString(source_dict[glassoc][0])
      if source_dict[glassoc][1] != "N" :
        fava.tevcat = True
    else:
      fava.association_source_class=TObjString("")
    

    ## try to to find the Z! 

    if assoc in zs: 
      fava.z = zs[assoc]
    else: 
      encodeme = {'name': js[i]['assoc']}; 
      ned_url="http://ned.ipac.caltech.edu/srs/ObjectLookup?"+urllib.urlencode(encodeme); 
      ned_page = urllib.urlopen(ned_url); 
      ned_js = json.loads(ned_page.read()) 

      fava.z = -1; 

      if 'Preferred' in ned_js: 
        if 'Redshift' in ned_js['Preferred']:
          if ned_js['Preferred']['Redshift']['Value']!=None:
            fava.z = float(ned_js['Preferred']['Redshift']['Value']) 

      sys.stdout.write('\rOn object: %s, z = %s' % (assoc, repr(fava.z)))
      sys.stdout.flush()
      zs[assoc] = fava.z 

    t.Fill() 




t.Write()


