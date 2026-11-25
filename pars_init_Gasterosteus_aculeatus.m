function [par, metaPar, txtPar] = pars_init_Gasterosteus_aculeatus(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 3218.1584;  free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.66543;      free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.23033;  free.kap_X = 1;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.028902;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.88028;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 169.521;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5250.356;   free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 6.441e-01; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 4.492e+00; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 4.977e+02; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 8.625e-09;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.DOa = 1.3682;     free.DOa   = 1;   units.DOa = '-';          label.DOa = 'alpha'; 
par.DOb = 5.702;      free.DOb   = 1;   units.DOb = '-';          label.DOb = 'beta'; 
par.DOcrit = 11.9747;  free.DOcrit = 1;   units.DOcrit = '-';       label.DOcrit = 'critical oxygen saturation'; 
par.E_Hh = 3.724e-01; free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.E_Hp_m = 5.211e+02; free.E_Hp_m = 1;   units.E_Hp_m = 'J';       label.E_Hp_m = 'maturity at puberty for males'; 
par.K = 1.203;        free.K     = 1;   units.K = '-';            label.K = 'correction factor'; 
par.T_AH = 45756.7079;  free.T_AH  = 1;   units.T_AH = 'K';         label.T_AH = 'Arrhenius temperature for upper boundary'; 
par.T_AL = 6252.2967;  free.T_AL  = 1;   units.T_AL = 'K';         label.T_AL = 'Arrhenius temperature for lower boundary'; 
par.T_H = 298.15;     free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'upper boundary'; 
par.T_L = 276.15;     free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'lower boundary'; 
par.del_M = 0.20134;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient for standard length'; 
par.del_Me = 0.15387;  free.del_Me = 1;   units.del_Me = '-';       label.del_Me = 'shape coefficient for standard lengthof embryos'; 
par.del_T = 0.22509;  free.del_T = 1;   units.del_T = '-';        label.del_T = 'shape coefficient for total length'; 
par.f = 0.91117;      free.f     = 1;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_1 = 1.5103;     free.f_1   = 1;   units.f_1 = '-';          label.f_1 = 'scaled functional response for WWp_f1'; 
par.f_2 = 2.1428;     free.f_2   = 1;   units.f_2 = '-';          label.f_2 = 'scaled functional response for WWp_f2'; 
par.f_DOJO = 0.99906;  free.f_DOJO = 1;   units.f_DOJO = '-';       label.f_DOJO = 'scaled functional response for DOJO data'; 
par.f_Leloutre2016 = 1.0644;  free.f_Leloutre2016 = 1;   units.f_Leloutre2016 = '-';  label.f_Leloutre2016 = 'scaled functional response for Leloutre2016 data'; 
par.f_LtN_Wootton1973 = 0.64267;  free.f_LtN_Wootton1973 = 1;   units.f_LtN_Wootton1973 = '-';  label.f_LtN_Wootton1973 = 'scaled functional response for Wootton1973 data'; 
par.f_TJO = 1.8041;   free.f_TJO = 1;   units.f_TJO = '-';        label.f_TJO = 'scaled functional response for TJO data'; 
par.f_WwJO_Andrea20 = 1.0435;  free.f_WwJO_Andrea20 = 1;   units.f_WwJO_Andrea20 = '-';  label.f_WwJO_Andrea20 = 'scaled functional response for WwJO data from Andrea2017'; 
par.f_WwJO_Brafield10 = 0.81582;  free.f_WwJO_Brafield10 = 1;   units.f_WwJO_Brafield10 = '-';  label.f_WwJO_Brafield10 = 'scaled functional response for WwJO data from Brafield1976'; 
par.f_WwJX = 1;       free.f_WwJX = 0;   units.f_WwJX = '-';       label.f_WwJX = 'scaled functional response for WwJX data'; 
par.f_WwJX_Wootton = 1;  free.f_WwJX_Wootton = 0;   units.f_WwJX_Wootton = '-';  label.f_WwJX_Wootton = 'scaled functional response for WwJX data'; 
par.f_tL15_f1 = 0.33887;  free.f_tL15_f1 = 1;   units.f_tL15_f1 = '-';    label.f_tL15_f1 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL15_f2 = 0.55569;  free.f_tL15_f2 = 1;   units.f_tL15_f2 = '-';    label.f_tL15_f2 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL15_f3 = 0.91805;  free.f_tL15_f3 = 1;   units.f_tL15_f3 = '-';    label.f_tL15_f3 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL15_f4 = 1.1195;  free.f_tL15_f4 = 1;   units.f_tL15_f4 = '-';    label.f_tL15_f4 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL19_f1 = 0.84715;  free.f_tL19_f1 = 1;   units.f_tL19_f1 = '-';    label.f_tL19_f1 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL19_f2 = 0.9454;  free.f_tL19_f2 = 1;   units.f_tL19_f2 = '-';    label.f_tL19_f2 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL7_f1 = 0.6105;  free.f_tL7_f1 = 1;   units.f_tL7_f1 = '-';     label.f_tL7_f1 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL7_f2 = 2.0675;  free.f_tL7_f2 = 1;   units.f_tL7_f2 = '-';     label.f_tL7_f2 = 'scaled functional response for Leloutre2016 data'; 
par.f_tL_JoneHyne1950 = 1.2068;  free.f_tL_JoneHyne1950 = 1;   units.f_tL_JoneHyne1950 = '-';  label.f_tL_JoneHyne1950 = 'scaled functional response for JoneHyne1950 data'; 
par.f_tL_WrigWoot2004 = 1.1755;  free.f_tL_WrigWoot2004 = 1;   units.f_tL_WrigWoot2004 = '-';  label.f_tL_WrigWoot2004 = 'scaled functional response for WrigWoot2004 data'; 
par.memb = 0.063006;  free.memb  = 1;   units.memb = '-';         label.memb = 'embryo mortality'; 
par.t0 = 40.2747;     free.t0    = 1;   units.t0 = 'd';           label.t0 = 'time since hatch at start if experiment'; 
par.z_m = 0.59971;    free.z_m   = 1;   units.z_m = '-';          label.z_m = 'zoom factor for males'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
