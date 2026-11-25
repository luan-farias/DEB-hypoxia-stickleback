function [data, auxData, metaData, txtData, weights] = mydata_Gasterosteus_aculeatus

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Perciformes'; 
metaData.family     = 'Gasterosteidae';
metaData.species    = 'Gasterosteus_aculeatus'; 
metaData.species_en = 'Three-spined stickleback'; 

metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MA', 'MP', 'TH'};
metaData.ecoCode.habitat = {'piMc', '0iFe', '0iFl', '0iFp'};
metaData.ecoCode.embryo  = {'Fnm'};
metaData.ecoCode.migrate = {'Ma'};
metaData.ecoCode.food    = {'biCi'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(20); % K, body temp
metaData.data_0     = {'ah'; 'ab'; 'ap'; 'am'; 'Lh'; 'Lb'; 'Lp'; 'Lp_f'; 'Lp_m'; 'Li'; 'Li_f'; 'Li_m'; 'Wwb'; 'Wwp'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'t-L_f'; 'L-N'; 'Ww-JO_T'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};    
metaData.date_subm = [2014 04 24];              
metaData.email    = {'bas.kooijman@vu.nl'};            
metaData.address  = {'VU University, Amsterdam'};   

metaData.author_mod_1   = {'Bas Kooijman'};    
metaData.date_mod_1     = [2016 10 30];              
metaData.email_mod_1    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_1  = {'VU University, Amsterdam'};   

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'sta@akvaplan.niva.no'}; 
metaData.date_acc    = [2016 10 30]; 


%% set data
% zero-variate data

data.ah = 12.9;     units.ah = 'd';     label.ah = 'age at hatching';       % bibkey.ah = 'Leloutre2016';
  temp.ah = C2K(13.1); units.temp.ah = 'K'; label.temp.ah = 'temperature';
% data.thb = 48/24;   units.thb = 'd';    label.thb = 'interval between hatching and birth';  % bibkey.ah = 'Leloutre2016';
%   temp.thb = C2K(16); units.temp.thb = 'K'; label.temp.thb = 'temperature';
data.thb = 4;   units.thb = 'd';    label.thb = 'interval between hatching and birth';  % bibkey.ah = 'Swarup';
  temp.thb = C2K(19); units.temp.thb = 'K'; label.temp.thb = 'temperature';

data.aj = 22;     units.aj = 'd';  label.aj = 'age at metamorphosis';       % bibkey.aj = 'Swarup';   
  temp.aj = C2K(19);  units.temp.aj = 'K'; label.temp.aj = 'temperature'; 
data.ap = 4*30.5;     units.ap = 'd'; label.ap = 'age at puberty';            bibkey.ap = 'Paep1983';
  temp.ap = C2K(20);  units.temp.ap = 'K'; label.temp.ap = 'temperature';
data.am = 5*365;     units.am = 'd';  label.am = 'life span';                 % bibkey.am = 'Wootton1976';   
  temp.am = C2K(20);  units.temp.am = 'K'; label.temp.am = 'temperature'; 

% data.Lh = 5.72/10;     units.Lh = 'cm';    label.Lh = 'length at hatching';    % bibkey.Lh = 'deKermoysan2013';
data.Lh = 4.3/10;     units.Lh = 'cm';    label.Lh = 'length at hatching';    % bibkey.Lh = 'deKermoysan2013';
data.Lb = 0.58;     units.Lb = 'cm';    label.Lb = 'length at birth';    % bibkey.Lb = 'deKermoysan2013';
data.Lj = 1.1;     units.Lj = 'cm';    label.Lj = 'length at metamorphosis';    % bibkey.Lj = 'Swarup';
data.Lp = 30.3/10;   units.Lp = 'cm';  label.Lp = 'standard lenght at puberty for females';  % bibkey.Lp_f = 'Leloutre2016';
data.Lp_m = 33.1/10;   units.Lp_m = 'cm';  label.Lp_m = 'standard lenght at puberty for males';    % bibkey.Lp_m = 'Leloutre2016';
data.Li = 61.88/10;  units.Li = 'cm';  label.Li = 'ultimate standard length for females';    % bibkey.Li_f = 'Leloutre2016';
data.Li_m = 55.51/10;  units.Li_m = 'cm';  label.Li_m = 'ultimate standard length for males';      % bibkey.Li_m = 'Leloutre2016';

data.Wd0 = 3e-4;   units.Wd0 = 'g';   label.Wd0 = 'initial (egg) dry weight'; bibkey.Wwb = 'Woot1973';
data.Wwb = data.Lb^3 * 0.0057;   units.Wwb = 'g';  label.Wwb = 'wet weight at birth';      bibkey.Wwb = 'Leloutre2016';
  comment.Wwb = 'computed as Lb^3 * 0.0057';
data.Wwp = (data.Lp* 0.249)^3;  units.Wwp = 'g';   label.Wwp = 'wet weight at puberty';    bibkey.Wwp = 'Woot1973';
  comment.Wwp = 'computed as (Lp * 0.249)^3';
data.Wwp_m = (data.Lp_m* 0.249)^3;  units.Wwp_m = 'g';   label.Wwp_m = 'wet weight at puberty for males';    bibkey.Wwp_m = 'Woot1973';
  comment.Wwp_m = 'computed as (Lp_m * 0.249)^3';
data.Wwi = (data.Li* 0.249)^3;   units.Wwi = 'g';  label.Wwi = 'ultimate wet weight';  bibkey.Wwi = 'Leloutre2016';
  comment.Wwi = 'computed as (Li * 0.249)^3';

data.Wwp_f1 = 0.574;  units.Wwp_f1 = 'g';   label.Wwp_f1 = 'wet weight at puberty at f1';    bibkey.Wwp_f1 = 'Woot1973';
data.Wwp_f2 = 0.666;  units.Wwp_f2 = 'g';   label.Wwp_f2 = 'wet weight at puberty at f2';    bibkey.Wwp_f2 = 'Woot1973';
% data.Wwp_f3 = 0.831;  units.Wwp_f3 = 'g';   label.Wwp_f3 = 'wet weight at puberty at f3';    bibkey.Wwp_f3 = 'Woot1973';

data.Ri = 100 * 10 / 365; units.Ri  = '#/d';   label.Ri  = 'maximum reprod rate';   bibkey.Ri  = 'Paep1983';   
  temp.Ri = C2K(20);  units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

% uni-variate data
% t-L data from WrigWoot2004 at 19 C
data.tL_WrigWoot2004 = [... % time since hatch (d), standard length (mm)
0.556	4.320
0.556	4.737
0.556	5.155
0.556	5.504
2.222	5.784
2.222	6.202
2.222	6.411
1.944	6.690
4.167	5.996
3.889	6.413
3.889	6.692
3.889	7.040
5.556	6.554
5.833	6.903
5.833	7.321
6.389	7.670
6.389	7.949
9.722	8.161
9.722	8.510
10.000	8.928
10.000	9.416
13.611	10.325
13.611	10.534
13.611	10.743
13.611	11.161
17.500	11.932
18.056	12.351
18.056	12.560
18.056	12.978
18.056	13.326
21.389	13.190
21.389	13.539
21.667	13.887
21.667	14.584
21.667	15.141
25.833	14.032
25.833	14.519
25.556	15.076
25.278	15.563
25.000	15.981
30.556	15.918
30.556	16.405
30.000	16.753
29.444	17.101
29.444	18.006
34.167	16.828
34.167	17.106
33.611	18.011
32.778	19.334
34.167	20.450
37.778	18.225
37.778	18.852
37.500	20.175
36.944	21.289
46.111	20.464
46.111	22.066
45.278	23.180
45.278	23.876
53.889	23.189
53.889	23.468
53.889	24.652
53.611	25.627
61.944	24.383
61.667	24.801
61.389	26.403
61.389	28.492
69.444	25.924
69.444	26.273
69.444	27.596
69.167	29.546
97.778	29.440
97.778	29.788
97.500	32.366
97.500	32.505
96.944	34.524
125.556	34.488
125.278	36.717
125.278	37.831
125.556	40.966]; % cm, physical length at f and T
for i=2:size(data.tL_WrigWoot2004,1);if data.tL_WrigWoot2004(i,1)<=data.tL_WrigWoot2004(i-1,1);data.tL_WrigWoot2004(i,1)=data.tL_WrigWoot2004(i-1,1)+1e-3;end;end
data.tL_WrigWoot2004(:,2) = data.tL_WrigWoot2004(:,2)/ 10; % convert mm to cm
units.tL_WrigWoot2004   = {'d', 'cm'};  label.tL_WrigWoot2004 = {'time since hatching', 'standard length'};  
temp.tL_WrigWoot2004    = C2K(19);  units.temp.tL_WrigWoot2004 = 'K'; label.temp.tL_WrigWoot2004 = 'temperature';
bibkey.tL_WrigWoot2004 = 'WrigWoot2004';

% t-L data from JoneHyne1950 at 19 C?
data.tL_JoneHyne1950 = [ ... % time (d), length (mm)
0   	14.886
29.925	26.547
60.836	31.065
92.191	36.308
121.890	42.484
153.003	40.068
183.394	41.482
213.416	44.553
242.951	45.555
273.760	46.865
304.305	53.143
334.820	58.491
365.355	51.108
395.527	45.589];
data.tL_JoneHyne1950(:,2) = data.tL_JoneHyne1950(:,2)/ 10; % convert mm to cm
units.tL_JoneHyne1950   = {'d', 'cm'};  label.tL_JoneHyne1950 = {'time', 'standard length'};  
temp.tL_JoneHyne1950    = C2K(20);  units.temp.tL_JoneHyne1950 = 'K'; label.temp.tL_JoneHyne1950 = 'temperature';
bibkey.tL_JoneHyne1950 = 'JoneHyne1950';

% t-L data from Leloutre2016 at 293.17 K
data.tL_Leloutre2016 = [... % time since hatch (d), standard length (mm)
1	6.325233333
8	9.163663462
15	13.36921154
22	17.14696078];
data.tL_Leloutre2016(:,2) = data.tL_Leloutre2016(:,2)/ 10; % convert mm to cm
units.tL_Leloutre2016   = {'d', 'cm'};  label.tL_Leloutre2016 = {'time since hatching', 'standard length'};  
temp.tL_Leloutre2016    = 293;  units.temp.tL_Leloutre2016 = 'K'; label.temp.tL_Leloutre2016 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 7.61°C and f1 (f1<f2)
data.tL7_f1 = [... % time since birth (d), standard length (mm)
0	5.669105
4.2	6.16815];
data.tL7_f1(:,2) = data.tL7_f1(:,2)/ 10; % convert mm to cm
units.tL7_f1   = {'d', 'cm'};  label.tL7_f1 = {'time since birth', 'standard length'};  
temp.tL7_f1    = C2K(7.61);  units.temp.tL7_f1 = 'K'; label.temp.tL7_f1 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 7.61°C and f2 (f1<f2)
data.tL7_f2 = [... % time since birth (d), standard length (mm)
0	5.6728
4.2	6.470222];
data.tL7_f2(:,2) = data.tL7_f2(:,2)/ 10; % convert mm to cm
units.tL7_f2   = {'d', 'cm'};  label.tL7_f2 = {'time since birth', 'standard length'};  
temp.tL7_f2    = C2K(7.61);  units.temp.tL7_f2 = 'K'; label.temp.tL7_f2 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 19.61°C and f1 (f1<f2)
data.tL19_f1 = [... % time since birth (d), standard length (mm)
0	5.702619
4.2	6.878063];
data.tL19_f1(:,2) = data.tL19_f1(:,2)/ 10; % convert mm to cm
units.tL19_f1   = {'d', 'cm'};  label.tL19_f1 = {'time since birth', 'standard length'};  
temp.tL19_f1    = C2K(19.61);  units.temp.tL19_f1 = 'K'; label.temp.tL19_f1 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 19.61°C and f2 (f1<f2)
data.tL19_f2 = [... % time since birth (d), standard length (mm)
0	5.722842
4.2	7.086611];
data.tL19_f2(:,2) = data.tL19_f2(:,2)/ 10; % convert mm to cm
units.tL19_f2   = {'d', 'cm'};  label.tL19_f2 = {'time since birth', 'standard length'};  
temp.tL19_f2    = C2K(19.61);  units.temp.tL19_f2 = 'K'; label.temp.tL19_f2 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 15.5°C and f1 (f1<f2<f3<f4)
data.tL15_f1 = [... % time since birth (d), standard length (mm)
0	5.97
1.98	6.12
5.96	6.81];
data.tL15_f1(:,2) = data.tL15_f1(:,2)/ 10; % convert mm to cm
units.tL15_f1   = {'d', 'cm'};  label.tL15_f1 = {'time since birth', 'standard length'};  
temp.tL15_f1    = C2K(15.5);  units.temp.tL15_f1 = 'K'; label.temp.tL15_f1 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 15.5°C and f2 (f1<f2<f3<f4)
data.tL15_f2 = [... % time since birth (d), standard length (mm)
0	5.95
1.98	6.47
5.96	7.61];
data.tL15_f2(:,2) = data.tL15_f2(:,2)/ 10; % convert mm to cm
units.tL15_f2   = {'d', 'cm'};  label.tL15_f2 = {'time since birth', 'standard length'};  
temp.tL15_f2    = C2K(15.5);  units.temp.tL15_f2 = 'K'; label.temp.tL15_f2 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 15.5°C and f3 (f1<f2<f3<f4)
data.tL15_f3 = [... % time since birth (d), standard length (mm)
0	5.97
1.98	6.98
5.96	8.19];
data.tL15_f3(:,2) = data.tL15_f3(:,2)/ 10; % convert mm to cm
units.tL15_f3   = {'d', 'cm'};  label.tL15_f3 = {'time since birth', 'standard length'};  
temp.tL15_f3    = C2K(15.5);  units.temp.tL15_f3 = 'K'; label.temp.tL15_f3 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% t-L data from Leloutre2016 at 15.5°C and f4 (f1<f2<f3<f4)
data.tL15_f4 = [... % time since birth (d), standard length (mm)
0	6.02
1.98	7.1
5.96	8.3];
data.tL15_f4(:,2) = data.tL15_f4(:,2)/ 10; % convert mm to cm
units.tL15_f4   = {'d', 'cm'};  label.tL15_f4 = {'time since birth', 'standard length'};  
temp.tL15_f4    = C2K(15.5);  units.temp.tL15_f4 = 'K'; label.temp.tL15_f4 = 'temperature';
% bibkey.tL_Leloutre2016 = 'Leloutre2016';

% % Salinity dependent data
% % t-L data from DeFaveri2013 at 18 C and 5 PSU
% % control
% data.tL_DeFaveri2013_5 = [ ... % time since fertilization (d), length (cm)
% 48	1.532
% 112	3.417
% 172	3.938
% 232	4.33];
% units.tL_DeFaveri2013_5   = {'d', 'cm'};  label.tL_DeFaveri2013_5 = {'time since fertilization', 'standard length'};  
% temp.tL_DeFaveri2013_5    = C2K(18);  units.temp.tL_DeFaveri2013_5 = 'K'; label.temp.tL_DeFaveri2013_5 = 'temperature';
% % sal.tL_DeFaveri2013_5 = 5; units.sal.tL_DeFaveri2013_5 = 'PSU'; label.sal.tL_DeFaveri2013_5 = 'salinity';
% % bibkey.tL_DeFaveri_5 = 'DeFaveri2013';

% L-N data from Leloutre2016
data.LN_Leloutre2016 = [ ... % standard length (mm), clutch size (eggs)
40.5	50.51
41	    45.66
42.42	69.9
44.77	87.68
44.78	67.88
45.58	63.84
45.97	130.1
46.77	88.08
46.85	75.15
47.37	112.73
47.38	82.42
47.43	90.1
47.75	119.6
47.87	105.86
48.27	110.71
48.63	97.78
48.7	88.48
49	    102.63
49.1	98.59
49.85	139.8
49.95	108.69
49.97	116.77
50.37	115.96
50.38	66.67
50.95	88.89
50.97	129.7
50.97	125.66
50.97	121.62
51.16	74.75
51.45	148.69
52.1	99.8
52.68	110.71
53.22	137.78
58.45	170.51];
data.LN_Leloutre2016(:,1) = data.LN_Leloutre2016(:,1)/10; % convert mm to cm
data.LN_Leloutre2016(:,2) = data.LN_Leloutre2016(:,2)*10; % convert eggs/spwaning to total reproduction period (year)
units.LN_Leloutre2016   = {'cm', '#'};  label.LN_Leloutre2016 = {'standard length', 'fecundity'};  
temp.LN_Leloutre2016    = C2K(18);  units.temp.LN_Leloutre2016 = 'K'; label.temp.LN_Leloutre2016 = 'temperature';
% bibkey.LN_Leloutre2016 = 'Leloutre2016';

% Lt-N data from Wootton1973
data.LtN_Wootton1973 = [ ... % total length (mm), spawn size (eggs)
39.938	47.336
40.868	56.557
40.868	48.566
41.86	52.869
41.86	47.336
42.913	61.475
43.905	60.246
44.897	54.713
44.897	64.549
44.897	68.238
44.897	72.541
45.826	96.516
45.888	69.467
45.888	63.32
45.888	55.943
46.88	111.27
46.942	59.631
46.942	76.844
46.942	81.762
47.004	44.262
47.872	79.303
47.872	42.418
48.74	106.352
48.802	110.656
48.988	51.639
49.855	76.844
49.917	67.008
50.909	82.377
52.769	119.877
52.831	147.541
54.752	124.795
60.579	280.328
73.471	298.77];
data.LtN_Wootton1973(:,1) = data.LtN_Wootton1973(:,1)/10; % convert mm to cm
data.LtN_Wootton1973(:,2) = data.LtN_Wootton1973(:,2)*10; % % convert eggs/spwaning to total reproduction period (year)
units.LtN_Wootton1973   = {'cm', '#'};  label.LtN_Wootton1973 = {'total length', 'fecundity'};  
temp.LtN_Wootton1973    = C2K(18);  units.temp.LtN_Wootton1973 = 'K'; label.temp.LtN_Wootton1973 = 'temperature';
% bibkey.LtN_Wootton1973 = 'Wootton1973';

% % wet weight - respiration data
% data.WwJO_Andrea20 = [...  % wet mass (g), 	respiration (g O2/ d)	
% 0.91	0.023228571
% 1.27	0.019817143
% 1.335	0.018476571
% 1.464	0.019532571
% 2.65	0.022062857
% 3.15	0.029050286
% ];
% data.WwJO_Andrea20(:,2) = data.WwJO_Andrea20(:,2) ./ data.WwJO_Andrea20(:,1); 
% units.WwJO_Andrea20  = {'g', 'g O2/d'};  label.WwJO_Andrea20 = {'wet weight', 'O_2 consumption', '20 C'};  
% temp.WwJO_Andrea20   = C2K(20);  units.temp.WwJO_Andrea20 = 'K'; label.temp.WwJO_Andrea20 = 'temperature';
% % bibkey.WwJO_Andrea20 = 'Andrea2017';

data.WwJO_Brafield10 = [...  % wet mass (g), 	respiration (mg O2/ d)	
0.294	1.032
0.488	1.584
0.561	2.064
0.632	2.736
0.657	4.944
0.676	2.136
0.694	1.896
0.701	2.856
0.72	2.544
0.773	3.552
0.823	3.84
0.869	3.696
0.88	2.928
0.914	4.776
1.075	3.384
1.181	6.384
1.192	3.168
1.193	3.6
1.95	5.16
2.061	14.688
2.402	6.624
];
data.WwJO_Brafield10(:,2) = data.WwJO_Brafield10(:,2)/1000; % convert mg O2/d to gO2/d
data.WwJO_Brafield10(:,2) = data.WwJO_Brafield10(:,2)./data.WwJO_Brafield10(:,1); % convert mg O2/d to gO2/d
units.WwJO_Brafield10  = {'g', 'g O2/d'};  label.WwJO_Brafield10 = {'wet weight', 'O_2 consumption', '10 C'};  
temp.WwJO_Brafield10   = C2K(10);  units.temp.WwJO_Brafield10 = 'K'; label.temp.Brafield10 = 'temperature';
% bibkey.WwJO_Brafield10 = 'Brafield1976';

% temperature - respiration
data.TJO = [...  % temperature (C), 	respiration (mg O2/kg/h)	
5   108.69
12  174.87
20  335.49
];
units.TJO  = {'C', 'mg O2/kg/h'};  label.TJO = {'temperature', 'O_2 consumption'};  
% bibkey.TJO = 'Brafield1976';

% DO - respiration
data.DOJO = [...  % DO (sat %), 	respiration (mg O2/kg/h)	
79.3	340
59	    363
48.9	182
38.8	221
28.7	208
18.6	117
12.8	67
];
temp.DOJO = C2K(15); units.temp.DOJO = 'K'; label.temp.DOJO = 'temperature';
units.DOJO  = {'%', 'mg O2/kg/h'};  label.DOJO = {'oxygen saturation', 'O_2 consumption'};  
% bibkey.TJO = 'Wassenaar2024';

% Wootton, 1980
% wet weight - food consumption at 19°C
data.WwJX_Wootton19 = [...  % wet mass (mg), 	max consumption (mg/ d)	
116.64	7.502857143
144.83	7.978571429
57.73	7.864285714
66.54	9.851428571
92.18	11.63714286
212.64	10.00142857
232.63	11.66857143
280.52	12.38142857
305.98	13.25428571
370.47	12.69428571
353.67	12.97285714
429.58	13.28714286
497.18	14.23571429
328.75	14.76285714
396.45	16.22857143
362.85	16.78714286
365.85	17.74142857
380.02	18.37571429
391.32	18.73285714
389.31	22.78571429
439.76	22.22714286
562.15	16.1
613.37	19.43428571
779.19	19.90142857
1190.19	24.40714286
1256.59	33.50285714
];
% data.WwJX_Wootton19(:,2) = data.WwJX_Wootton19(:,2).*data.WwJX_Wootton19(:,1)./10; % convert g/g/d to g/d
data.WwJX_Wootton19 = data.WwJX_Wootton19/1000; % convert mg to g
units.WwJX_Wootton19  = {'g', 'g/d'};  label.WwJX_Wootton19 = {'wet weight', 'feeding rate', '19 C'};  
temp.WwJX_Wootton19   = C2K(19);  units.temp.WwJX_Wootton19 = 'K'; label.temp.WwJX_Wootton19 = 'temperature';
% bibkey.WwJX = 'Wootton1980';

% wet weight - food consumption at 12.5°C 
data.WwJX_Wootton12 = [...  % wet mass (mg), 	max consumption (mg/ d)	
222.34	2.172857143
225.14	2.132857143
121.55	3.887142857
307.14	4.87
194.75	4.717142857
169.62	5.512857143
195	5.988571429
217.53	6.305714286
268.15	6.58
487.41	7.401428571
484.71	7.958571429
375.26	8.481428571
308.19	10.23428571
463.01	11.89428571
642.29	9.34
704.66	12.27714286
595.39	13.67428571
525.23	14.03571429
598.41	14.78714286
651.68	14.26714286
763.97	13.90285714
957.87	14.28857143
994.51	14.92142857
1145.91	13.48285714
936.14	18.10428571
1350.13	23.52428571
1309.16	29.44714286
1598.3	28.47571429
];
% data.WwJX_Wootton12(:,2) = data.WwJX_Wootton12(:,2).*data.WwJX_Wootton12(:,1)./10; % convert g/g/d to g/d
data.WwJX_Wootton12 = data.WwJX_Wootton12/1000; % convert mg to g
units.WwJX_Wootton12  = {'g', 'g/d'};  label.WwJX_Wootton12 = {'wet weight', 'feeding rate', '12.5 C'};  
temp.WwJX_Wootton12   = C2K(12.5);  units.temp.WwJX_Wootton12 = 'K'; label.temp.WwJX_Wootton12 = 'temperature';
% bibkey.WwJX = 'Wootton1980';

% wet weight - food consumption at 7°C
data.WwJX_Wootton7 = [...  % wet mass (mg), 	max consumption (mg/ d)	
191.98	4.915714286
181.49	8.731428571
344.26	7.967142857
248.31	5.667142857
481.56	6.21
235.03	9.602857143
341.68	9.12
527.52	11.37285714
602.97	9.381428571
766	9.968571429
595.02	11.80571429
735.44	11.67857143
780.82	13.86142857
1306.2	14.38714286
1362.88	16.92714286
1929.76	14.19
2069.08	22.72571429
];
% data.WwJX_Wootton7(:,2) = data.WwJX_Wootton7(:,2).*data.WwJX_Wootton7(:,1)./10; % convert g/g/d to g/d
data.WwJX_Wootton7 = data.WwJX_Wootton7/1000; % convert mg to g
units.WwJX_Wootton7  = {'g', 'g/d'};  label.WwJX_Wootton7 = {'wet weight', 'feeding rate', '7 C'};  
temp.WwJX_Wootton7   = C2K(7);  units.temp.WwJX_Wootton7 = 'K'; label.temp.WwJX_Wootton7 = 'temperature';
% bibkey.WwJX = 'Wootton1980';


% % wet weight - food consumption at 23°C (optimum temp)
% data.WwJX = [...  % wet mass (g), 	max consumption (g / g/ d)	
% 0.796	0.447
% 0.798	0.745
% 0.897	0.445
% 0.9	    0.551
% 0.901	0.376
% 1.098	0.631
% 1.498	0.467
% 1.599	0.378
% 1.694	0.418
% 1.709	0.422
% 2.099	0.341
% 2.8 	0.384
% 2.896	0.345
% 3   	0.367
% 3.098	0.308];
% % data.WwJX(:,2) = data.WwJX(:,2).*data.WwJX(:,1)./10; % convert g/g/d to g/d
% units.WwJX  = {'g', 'g/g/d'};  label.WwJX = {'wet weight', 'feeding rate', '23 C'};  
% temp.WwJX   = C2K(23);  units.temp.WwJX = 'K'; label.temp.WwJX = 'temperature';
% % bibkey.WwJX = 'Hovel2015';

% % temperature - proportion of Cmax
% data.CmaxT = [...  % temperature (°C), 	proportion of Cmax	
% 7.6	0.328
% 13	0.457
% 18	0.755
% 23	1
% 25	0.096];
% % data.CmaxT(:,1) = C2K(data.CmaxT(:,1)); % convert temp from C to K
% units.CmaxT = {'C', '-'};  label.CmaxT = {'temperature', 'proportion of Cmax'};  
% % bibkey.WwJX = 'Hovel2015';

% temperature - age of hatching
data.Tah = [...  % temperature (°C), 	age of hatching (d)	
7.96	38.73
7.97	15.06
7.98	42.84
10	12.97
11.94	19.77
12.01	16.06
13.98	10.53
15.02	9.02
15.99	8.96
17.03	8
18.02	7.45
18.79	6.62
18.96	7.06
19.42	8.44
19.48	7.91
19.95	6.65
24.44	5
24.89	6.12
24.89	6.77
24.96	5
25.95	5.34];
% data.Tah(:,1) = C2K(data.Tah(:,1)); % convert temp from C to K
units.Tah = {'C', 'd'};  label.Tah = {'temperature', 'age at hatching'};  
% bibkey.WwJX = 'Wootton1976';

% DO - age of hatching
data.DOah = [...  % DO saturation (%), 	age of hatching (d)	
100	7.04
80	7.04
60	7.04
40	8.04
20	9.04
];
units.DOah = {'%', 'd'};  label.DOah = {'DO saturation', 'age at hatching'};  
temp.DOah   = C2K(19.5);  units.temp.DOah = 'K'; label.temp.DOah = 'temperature';
% bibkey.WwJX = 'Wootton1976';

% DO - Embryo Survival
DOS = [...  % time (days), survival fraction (/) for different DO saturation: 100, 80, 60, 40 and 20 %	
1.041666667	1	1	    1	1	    1
2.041666667	1	0.99	1	1	    0.98
3.041666667	1	0.99	1	1	    0.98
4.041666667	1	0.99	1	0.99	0.96
5.041666667	1	0.99	1	0.99	0.95
6.041666667	1	0.99	1	0.98	0.94
7.041666667	1	0.99	1	0.96	0.83
8.041666667	1	0.99	1	0.96	0.76
9.041666667	1	0.99	1	0.96	0.61];

data.DOS1 = DOS(:,[1,2]);
units.DOS1 = {'d', '%'};  label.DOS1 = {'days since fertilization', 'survival fraction', '100%'};  
temp.DOS1   = C2K(19.5);  units.temp.DOS1 = 'K'; label.temp.DOS1 = 'temperature';
% bibkey.WwJX = 'Wootton1976';
%
data.DOS2 = DOS(:,[1,3]);
units.DOS2 = {'d', '%'};  label.DOS2 = {'days since fertilization', 'survival fraction', '80%'};  
temp.DOS2   = C2K(19.5);  units.temp.DOS2 = 'K'; label.temp.DOS2 = 'temperature';
% bibkey.WwJX = 'Wootton1976';

data.DOS3 = DOS(:,[1,4]);
units.DOS3 = {'d', '%'};  label.DOS3 = {'days since fertilization', 'survival fraction', '60%'};  
temp.DOS3   = C2K(19.5);  units.temp.DOS3 = 'K'; label.temp.DOS3 = 'temperature';
% bibkey.WwJX = 'Wootton1976';
%
data.DOS4 = DOS(:,[1,5]);
units.DOS4 = {'d', '%'};  label.DOS4 = {'days since fertilization', 'survival fraction', '40%'};  
temp.DOS4   = C2K(19.5);  units.temp.DOS4 = 'K'; label.temp.DOS4 = 'temperature';
% bibkey.WwJX = 'Wootton1976';
%
data.DOS5 = DOS(:,[1,6]);
units.DOS5 = {'d', '%'};  label.DOS5 = {'days since fertilization', 'survival fraction', '20%'};  
temp.DOS5   = C2K(19.5);  units.temp.DOS5 = 'K'; label.temp.DOS5 = 'temperature';
% bibkey.WwJX = 'Wootton1976';


%% set weights for all real data
weights = setweights(data, []);
weights.Lb = 3 * weights.Lb;
weights.Lh = 2 * weights.Lh;
weights.ah = 4 * weights.ah;
weights.ap = 0.5 * weights.ap;
% weights.DOS1 = 0 * weights.DOS1;
% weights.DOS2 = 0 * weights.DOS2;
% weights.DOS3 = 0 * weights.DOS3;
% weights.DOS4 = 0 * weights.DOS4;
% weights.DOS5 = 0 * weights.DOS5;
% weights.WwJO_Andrea20 = 0 * weights.WwJO_Andrea20;
% weights.WwJO_Brafield10 = 0 * weights.WwJO_Brafield10;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Group plots
set1 = {'tL_WrigWoot2004'}; subtitle1 = {'Data from WrigWoot2004'};
set2 = {'tL_JoneHyne1950'}; subtitle2 = {'Data from JoneHyne1950'};
set3 = {'tL_Leloutre2016'}; subtitle3 = {'Data from Leloutre2016'};
set4 = {'LN_Leloutre2016'}; subtitle4 = {'Data from Leloutre2016'};
set5 = {'LtN_Wootton1973'}; subtitle5 = {'Data from Wootton1973'};
set6 = {'WwJO_Brafield10'}; subtitle6 = {'JO at 10°C'};
% set7 = {'WwJX'}; subtitle7 = {'JX at 23°C'};
set8 = {'WwJX_Wootton19', 'WwJX_Wootton12', 'WwJX_Wootton7'}; 
subtitle8 = {'JX at 19°C', 'JX at 12.5°C', 'JX at 7°C'};

set9 = {'tL7_f1','tL7_f2','tL19_f1','tL19_f2'}; subtitle9 = {'t-L at 7°C','t-L at 19°C'};

set10 = {'tL15_f1','tL15_f2','tL15_f3','tL15_f4'}; subtitle10 = {'t-L at 15°C'};
set11 = {'DOJO'}; subtitle11 = {'DOJO'};
set12 = {'DOah'}; subtitle12 = {'DOah'};
set13 = {'DOS1', 'DOS2', 'DOS3', 'DOS4', 'DOS5'}; subtitle13 = {'survival at different %DO'}; %
metaData.grp.sets = {set1,set2,set3,set9,set10,set4,set5,set6,set8,set12,set13};
metaData.grp.subtitle = {subtitle1,subtitle2,subtitle3,subtitle9,...
    subtitle10,subtitle4,subtitle5,subtitle6,subtitle8,subtitle12,subtitle13};

%% Facts
F1 = 'SL = 0.906 * TL + 0.200, r^2 = 0.996';
metaData.bibkey.F1 = 'WrigWoot2004'; 
F2 = 'egg weight increases with feeding rate; weight at maturity increases with feeding rate';
metaData.bibkey.F2 = 'Woot1973'; 
metaData.facts = struct('F1',F1,'F2',F2);

%% Links
metaData.links.id_CoL = '3FCKT'; % Cat of Life
metaData.links.id_ITIS = '166365'; % ITIS
metaData.links.id_EoL = '223856'; % Ency of Life
metaData.links.id_Wiki = 'Gasterosteus_aculeatus'; % Wikipedia
metaData.links.id_ADW = 'Gasterosteus_aculeatus'; % ADW
metaData.links.id_Taxo = '44801'; % Taxonomicon
metaData.links.id_WoRMS = '126505'; % WoRMS
metaData.links.id_fishbase = 'Gasterosteus-aculeatus'; % fishbase


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Gasterosteus_aculeatus}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'WrigWoot2004'; type = 'Article'; bib = [ ... 
'author = {H. A. Wright and R. J. Wootton and I. Barber}, ' ... 
'year = {2004}, ' ...
'title = {Interpopulation variation in early growth of threespine sticklebacks (\emph{Gasterosteus aculeatus}) under laboratory conditions}, ' ...
'journal = {Can. J. Fish. Aquat. Sci.}, ' ...
'volume = {61}, ' ...
'pages = {1832-1838}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Woot1973'; type = 'Article'; bib = [ ... 
'author = {R. J. Wootton}, ' ... 
'year = {1973}, ' ...
'title = {The effect of size of food ration on egg production in the female three-spined stickleback \emph{Gasterosteus aculeatus} {L}.}, ' ...
'journal = {J. Fish Biol.}, ' ...
'volume = {5}, ' ...
'pages = {89--96}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'JoneHyne1950'; type = 'Article'; bib = [ ... 
'author = {Jones, J. W. and Hynes, H. B. N.}, ' ... 
'year = {1950}, ' ...
'title = {The age and growth of \emph{Gasterosteus aculeatus}, \emph{Pygosteus pungitius} and \emph{Spinachia vulgaris} as shown by their otoliths}, ' ...
'journal = {J. Anim. Ecol.}, ' ...
'volume = {19}, ' ...
'pages = {59--73}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Paep1983'; type = 'Book'; bib = [ ... 
'author = {Paepke, H.-J.}, ' ... 
'year = {1983}, ' ...
'title = {Die  Stichlinge.}, ' ...
'series = {Die Neue Brehm-Buecherei}, ' ...
'volume = {10}, ' ...
'publisher = {A. Ziemsen Verlag DDR, WittenBerg Lutherstadt}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
% bibkey = 'AnAge'; type = 'Misc'; bib = ...
% 'howpublished = {\url{http://genomics.senescence.info/species/entry.php?species=Gasterosteus_aculeatus}}';
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ADW'; type = 'Misc'; bib = ...
'howpublished = {\url{http://animaldiversity.ummz.umich.edu/accounts/Gasterosteus_aculeatus/}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

