const int n_ANITA = 7;
Double_t ANITA_x[n_ANITA]       = {18, 18.5, 19, 19.5, 20, 20.5, 21};

Double_t ANITA_1_effArea[n_ANITA] = {3.13E-4,
				     1.60E-2,
				     4.58E-1,
				     4.60,
				     2.45E1,
				     7.31E1,
				     1.85E2};

/* 1021.5 eV 4.30 × 102 km2 sr */
/* 1022 eV 8.23 × 102 km2 sr */
/* 1022.5 eV 1.41 × 103 km2 sr */
/* 1023 eV 2.54 × 103 km2 sr */

//double ANITA_4_eff[n_ANITA] = { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
double ANITA_4_eff[n_ANITA] = { 0.35, 0.78, 0.834, 0.836, 0.838, 0.84, 0.84 };
double ANITA_3_eff[n_ANITA] = { 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84 };
double ANITA_2_eff[n_ANITA] = { 0.67664 , 0.60748 , 0.60000 , 0.55701 , 0.47477, 0.40935 , 0.41495 };
double ANITA_1_eff[n_ANITA] = { 0.62041 , 0.71429 , 0.75306 , 0.68571 , 0.64694, 0.50000 , 0.43878 };

double ANITA_4_livetime = 24.5416*24*3600.; // 24.516 days
double ANITA_3_livetime = 17.4*24*3600.; // 17.4 days
double ANITA_2_livetime = 28.5*24*3600.; 
double ANITA_1_livetime = 17.4*24*3600.;


// From ANITA-2 Erratum
Double_t ANITA_2_effArea_published[n_ANITA] = {0.00043,
					       0.05000,
					       0.92000,
					       6.60,
					       36.00,
					       108.00,
					       259.00};
// From SuperMongo code
Double_t ANITA_2_effArea_icemc2010[n_ANITA]={0.00021,
					     0.017,
					     0.4013,
					     3.97748,
					     21.7782,
					     81.2713,
					     234.776};

Double_t ANITA_2_effArea_Peter[n_ANITA]={0.00063,
					 0.05400,
					 0.97,
					 6.70,
					 44.50,
					 122.00,
					 279.00};


// Icemc January 2018
Double_t ANITA_2_effVol[n_ANITA] = { 0.612838,
				     19.2151,
				     184.361,
				     1146.8,
				     3918.96,
				     9684.73,
				     18955.6};

Double_t ANITA_3_effVol_published[n_ANITA] = { 0.28516,
					       5.79936,
					       101.24 ,
					       784.763,
					       3244.84,
					       9398.71,
					       19776.2}; // in km^3 srt

// Icemc December 2018
Double_t ANITA_3_effVol[] = {0.544934,
			     3.8882,
			     131.667,
			     837.749,
			     3249.79,
			     9031.45,
			     18516.6};

Double_t ANITA_4_effVol[n_ANITA] = { 2.43441,
				     12.0405,
				     139.711,
				     989.925,
				     4744.19,
				     14142.7,
				     30423.9};

// Icemc December 2018 (number at E18 TBC)
Double_t ANITA_4_vol_tuffsoff[n_ANITA] = {0.8,
					  5.33191,
					  172.263,
					  1488.34,
					  6529.09,
					  17998.5,
					  36258.9};

Double_t ANITA_3_effVol_noMasking[n_ANITA] = { 0.910887,
					       22.,
					       248.993,
					       1581.04,
					       5768.12,
					       15041.8,
					       29690.}; // in km^3 srt

Double_t ANITA_3_effVol_justTau[n_ANITA] = {0.0183866,
					    0.255914,
					    30.5467,
					    263.894,
					    1010.1,
					    3290.64,
					    7159.01};


Double_t intLength_RENO[n_ANITA] = { 1131.34,
				     793.952,
				     557.179,
				     391.016,
				     274.407,
				     192.573 ,
				     135.143 }; // E21

Double_t intLength_CONNOLLY_nuNC[n_ANITA] = { 3826.5,
					      2617.44,
					      1817.98,
					      1279.9,
					      912.003,
					      656.936,
					      477.867}; // E21

  
Double_t intLength_CONNOLLY_nuCC[n_ANITA] = { 1544.84,
					      1065.44,
					      745.54,
					      528.431,
					      378.864,
					      274.445,
					      200.67}; // E21

  
Double_t intLength_CONNOLLY_nubarNC[n_ANITA] = { 3868.47,
						 2647.45,
						 1840.89,
						 1297.89,
						 926.252,
						 668.193 ,
						 486.701 }; // E21
  
Double_t intLength_CONNOLLY_nubarCC[n_ANITA] = { 1671.17,
						 1151.12,
						 805.022,
						 570.46,
						 408.962,
						 296.22,
						 216.548}; // E21


Double_t intLength_CONNOLLY_nuCCup[n_ANITA] = { 1854.212,
						1429.790,
						1143.231,
						945.714,
						807.683,
						710.780,
						643.375}; // E21

Double_t intLength_CONNOLLY_nuCClow[n_ANITA] = { 1276.548,
						 800.129,
						 500.581,
						 312.417,
						 194.457,
						 120.701,
						 74.718}; // E21


double Auger15_x[8] = { 16.70533, 17.21212, 17.71369, 18.22048, 18.72205, 19.22362, 19.73041, 20.22675 }; // values in log(eV)

double Auger15_y[8] = {-6.82741, -7.44332, -7.66667, -7.58545, -7.32826, -6.98985, -6.60406, -6.17090 }; // values in log( E^2F ) where E^2F is [GeV cm-2 s-1 sr-1]

 

double IcecubeOld_x[17] = { 14.7551   ,
			    14.932 ,
			    15.1769    ,
			    15.3628    ,
			    15.5941    ,
			    15.8707    ,
			    16.2426    ,
			    16.6327    ,
			    17.0227    ,
			    17.4127    ,
			    17.8481    ,
			    18.3288    ,
			    18.7732    ,
			    19.195 ,
			    19.5941    ,
			    19.8707    ,
			    20.0884 };

double IcecubeOld_y[17] = { -12.1491  ,
			    -12.5539   ,
			    -13.0652   ,
			    -13.3262   ,
			    -13.5819   ,
			    -13.8642   ,
			    -14.1624   ,
			    -14.4767   ,
			    -14.7856   ,
			    -15.0839   ,
			    -15.3768   ,
			    -15.6751   ,
			    -15.9414   ,
			    -16.1598   ,
			    -16.3089   ,
			    -16.3782   ,
			    -16.3515 };



double IceCube2017x[]={16.096,
		       16.4966      ,
		       16.8971      ,
		       17.2977 ,
		       17.6982      ,
		       18.0988      ,
		       18.4993      ,
		       18.8944      ,
		       19.2949      ,
		       19.6955      ,
		       19.9973   };

  
double IceCube2017y[]={-15.0108      ,
		       -15.2504      ,
		       -15.5284      ,
		       -15.883 ,
		       -16.2903      ,
		       -16.5443      ,
		       -16.8415      ,
		       -17.0571      ,
		       -17.2153      ,
		       -17.3542      ,
		       -17.3734   };
