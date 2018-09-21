#! /usr/bin/env python

# Script to grab fava events
# Since I'm downloading from a web site and need to parse json, it was easier to just do it in python. 
# You need PyROOT obviously to run this... 



from ROOT import TFile,TTree, gROOT, TObjString
import urllib 
import json 
import datetime
import csv 


week_start = 300; 
week_end = 500; 


gROOT.ProcessLine("#include \"fava.h\"") 
from ROOT import FAVAEntry; 

fava = FAVAEntry() 


f = TFile("fava.root","recreate"); 

t = TTree("fava","FAVA tree"); 

t.Branch('fava', fava); 

fermi_start_date = datetime.datetime(2001,1,1)
time_offset = int(fermi_start_date.strftime("%s"))



catalog = open("3fgl.csv") 
reader = csv.DictReader(catalog) 


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
    fava.le_index =  float(js[i]['le_index'])
    fava.he_index =  float(js[i]['he_index'])
    fava.association = TObjString(js[i]['assoc'])
    glassoc = js[i]['fglassoc'] 
    fava.association_3fgl = TObjString(glassoc)
    fava.tevcat = False 
    if glassoc in source_dict: 
      fava.source_class = TObjString(source_dict[glassoc][0])
      if source_dict[glassoc][1] != "N" :
        fava.tevcat = True
    else:
      fava.association_source_class=TObjString("")
    

    t.Fill() 








t.Write()


