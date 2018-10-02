#! /usr/bin/env python

# Script to grab light curves for TeV blazars 
# Since I'm downloading from a web site and need to parse json, it was easier to just do it in python. 
# You need PyROOT obviously to run this... 


from ROOT import TFile,TGraph,TTree,TNamed
from array import array 
import urllib 
import json 
import datetime


# Dictionary of TeV blazars within +/-25 degrees dec (From TeVCat... someone should check that I didn't typo anyting) 
blazars = {} 


blazars["1ES 1440+122"]     = (15*(14 + 43./60 + 15./3600), 12 + 11./3600) 
blazars["1ES 1101-232"]     = (15*(11 + 3./60 + 36.5/3600), -(23 + 29./60 + 45./3600))
blazars["1ES 0414+009"]     = (15*(4 + 16./60 + 52.96/3600), (1 + 5./60 + 20.4/3600))
blazars["1ES 1741+196"]     = (15*(17 + 44./60 + 1.2/3600), (19 + 32./60 + 47./3600))
blazars["1ES 0229+200"]     = (15*(2 + 32./60 + 53.2/3600), (20 + 16./60 + 21./3600))
blazars["1ES 0347-121"]     = (15*(3 + 49./60 + 23.0/3600), -(11 + 58./60 + 38./3600))
blazars["3C 279"]           = (15*(12 + 56./60 + 11.1/3600), -(5 + 47./60 + 22./3600))
blazars["4C +21.35"]        = (15*(12 + 24./60 + 54.4/3600), (21 + 22./60 + 46./3600))
blazars["AP Lib"]           = (15*(15 + 17./60 + 41.8/3600), -(24 + 22./60 + 19./3600))
blazars["H 1722+119"]       = (15*(17 + 24./60 + 4.3/3600), (11 + 52./60 + 15./3600))
blazars["HESS J1943+213"]   = (15*(19 + 43./60 + 55./3600), (21 + 18./60 + 8./3600))
blazars["KUV 00311-1938"]   = (15*(0 + 33./60 + 34.2/3600), -(19 + 21./60 + 33./3600))
blazars["MS 1221.8+2452"]   = (15*(12 + 24./60 + 24.2/3600), (24 + 36./60 + 24./3600))
blazars["OJ 287"]           = (15*(8 + 54./60 + 49.1/3600), (20 + 5./60 + 58.89/3600))
blazars["OT 081"]           = (15*(17 + 51./60 + 32.82/3600), (9 + 39./60 + 0.73/3600))
blazars["PG 1553+113"]      = (15*(15 + 55./60 + 44.7/3600), (11 + 11./60 + 41./3600))
blazars["PKS 0301-243"]     = (15*(3 + 3./60 + 23.49/3600), -(24 + 7./60 + 35.86/3600))
blazars["PKS 0736+017"]     = (15*(7 + 39./60 + 18.0/3600), (1 + 37./60 + 5./3600))
blazars["PKS 1424+240"]     = (15*(14 + 27./60 + 0./3600), (23 + 47./60 + 40./3600))
blazars["PKS 1510-089"]     = (15*(15 + 12./60 + 52.2/3600), -(9 + 6./60 + 21.6/3600))
blazars["RBS 0413"]         = (15*(3 + 19./60 + 47/3600), (18 + 45./60 + 42./3600))
blazars["RBS 0723"]         = (15*(8 + 47./60 + 12.9/3600), (11 + 33./60 + 50./3600))
blazars["RGB J0152+017"]    = (15*(1 + 52./60 + 33.5/3600), (1 + 46./60 + 40.3/3600))
blazars["RGB J2243+203"]    = (15*(22 + 43./60 + 52/3600), (20 + 19./60 + 12./3600))
blazars["RX J0638.7+1516"]  = (15*(6 + 48./60 + 45.6/3600), (15 + 16./60 + 12./3600))
blazars["S2 0109+22"]       = (15*(1 + 12./60 + 5.8/3600), (22 + 44./60 + 39./3600))
blazars["SHBL J001355.9-185406"] = (15*(0 + 13./60 + 52./3600), -(18 + 53./60 + 29./3600))
blazars["TXS 0506+056"] =    (15*(5 + 9./60 + 25.96/3600), (5 + 41./60 + 35.33/3600))
blazars["VER J0421+211"] =    (15*(5 + 21./60 + 45./3600), (21 + 12./60 + 51.4/3600))


f = TFile("blazars.root","RECREATE"); 
t = TTree("blazars","blazars") 

g_he_counts = TGraph() 
g_he_significance = TGraph() 
g_he_bg = TGraph() 
g_le_counts = TGraph() 
g_le_significance = TGraph() 
g_le_bg = TGraph() 

g_he_counts.SetName("g_he_counts") 
g_he_significance.SetName("g_he_significance") 
g_he_bg.SetName("g_he_bg") 
g_le_counts.SetName("g_le_counts") 
g_le_significance.SetName("g_le_significance") 
g_le_bg.SetName("g_le_bg") 

g_he_counts.GetXaxis().SetName("time") 
g_he_significance.GetXaxis().SetName("time") 
g_he_bg.GetXaxis().SetName("time") 
g_le_counts.GetXaxis().SetName("time") 
g_le_significance.GetXaxis().SetName("time") 
g_le_bg.GetXaxis().SetName("time") 

g_he_counts.GetYaxis().SetName("counts") 
g_he_significance.GetYaxis().SetName("#sigma") 
g_he_bg.GetYaxis().SetName("counts") 
g_le_counts.GetYaxis().SetName("counts") 
g_le_significance.GetYaxis().SetName("#sigma") 
g_le_bg.GetYaxis().SetName("counts") 


name = TNamed("name","blazar") 
#yuck 
ra = array('d',[0.])
dec = array('d',[0.])

t.Branch("he_counts", g_he_counts); 
t.Branch("he_bg", g_he_bg); 
t.Branch("he_significance", g_he_significance); 
t.Branch("le_counts", g_le_counts); 
t.Branch("le_bg", g_le_bg); 
t.Branch("le_significance", g_le_significance); 
t.Branch("name", name) 
t.Branch("RA", ra, "RA/D") 
t.Branch("dec", dec, "dec/D") 


fermi_start_date = datetime.datetime(2001,1,1)
time_offset = int(fermi_start_date.strftime("%s"))

for blazar in blazars.keys(): 

  name.SetTitle(blazar) 
  ra[0] = blazars[blazar][0] 
  dec[0] = blazars[blazar][1] 
  webpage = urllib.urlopen("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/FAVA/queryDB_Lightcurve.php?ra=%f&dec=%f" % (ra[0],dec[0]))
  js = json.loads(webpage.read()) 
  N = len(js) 
  g_he_counts.Set(N) 
  g_le_counts.Set(N) 
  g_he_bg.Set(N) 
  g_le_bg.Set(N) 
  g_he_significance.Set(N) 
  g_le_significance.Set(N) 

  g_he_counts.SetTitle("HE Counts for %s" % (blazar))
  g_le_counts.SetTitle("LE Counts for %s" % (blazar))
  g_he_bg.SetTitle("HE BG Counts for %s" % (blazar))
  g_le_bg.SetTitle("LE BG Counts for %s" % (blazar))
  g_he_significance.SetTitle("HE Significance for %s" % (blazar))
  g_le_significance.SetTitle("LE Significance for %s" % (blazar))



  for i in range(N): 
    stupid_time = float(js[i]['time']) 
    tm = stupid_time + time_offset ;
    he_count = float(js[i]['he_nev'])
    le_count = float(js[i]['nev'])
    he_bg = float(js[i]['he_avnev'])
    le_bg = float(js[i]['avnev'])
    he_sig = float(js[i]['he_sigma'])
    le_sig = float(js[i]['sigma'])

    g_he_counts.SetPoint(i,tm,he_count) 
    g_le_counts.SetPoint(i,tm,le_count) 
    g_he_bg.SetPoint(i,tm,he_bg) 
    g_le_bg.SetPoint(i,tm,le_bg) 
    g_he_significance.SetPoint(i,tm,he_sig) 
    g_le_significance.SetPoint(i,tm,le_sig) 


  t.Fill() 

t.Write() 
