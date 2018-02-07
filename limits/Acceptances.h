const int n_ANITA = 7;
Double_t ANITA_x[n_ANITA]       = {18, 18.5, 19, 19.5, 20, 20.5, 21};

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
// Icemc January 2018
Double_t ANITA_3_effVol[n_ANITA] = { 0.28516,
				     5.79936,
				     101.24 ,
				     784.763,
				     3244.84,
				     9398.71,
				     19776.2}; // in km^3 srt


Double_t ANITA_3_effVol_noMasking[n_ANITA] = { 0.910887,
					       22.,
					       248.993,
					       1581.04,
					       5768.12,
					       15041.8,
					       29690.}; // in km^3 srt

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
