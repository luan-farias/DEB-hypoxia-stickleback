function [prdData, info] = predict_Gasterosteus_aculeatus(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  if (E_Hh >= E_Hb || T_ref > T_H || f_tL7_f1 > f_tL7_f2 || f_tL19_f1 > f_tL19_f2...
          || f_tL15_f1 > f_tL15_f2 || f_tL15_f2 > f_tL15_f3 || f_tL15_f3 > f_tL15_f4...
          || f_tL15_f4 > 2)
    info = 0; prdData = []; return
  end

  parsT=[T_A, T_L, T_H, T_AL, T_AH];

  % compute temperature correction factors
  TC_ah = tempcorr(temp.ah, T_ref, parsT);
  TC_hb = tempcorr(temp.thb, T_ref, parsT);
  TC_aj = tempcorr(temp.aj, T_ref, parsT);
  TC_ap = tempcorr(temp.ap, T_ref, parsT);
  TC_am = tempcorr(temp.am, T_ref, parsT);
  TC_Ri = tempcorr(temp.Ri, T_ref, parsT);
  TC_tL = tempcorr(temp.tL_WrigWoot2004, T_ref, parsT);
  TC_LN = tempcorr(temp.LN_Leloutre2016, T_ref, parsT);
  TC_JO_10 = tempcorr(temp.WwJO_Brafield10, T_ref, parsT);
  % TC_JO_20 = tempcorr(temp.WwJO_Andrea20, T_ref, parsT);
  TC_TJO = tempcorr(C2K(TJO(:,1)), T_ref, parsT);
  % TC_JX = tempcorr(temp.WwJX, T_ref, parsT);
  TC_JX_19 = tempcorr(temp.WwJX_Wootton19, T_ref, parsT);
  TC_JX_12 = tempcorr(temp.WwJX_Wootton12, T_ref, parsT);
  TC_JX_7 = tempcorr(temp.WwJX_Wootton7, T_ref, parsT);
  % TC_Cmax = tempcorr(C2K(CmaxT(:,1)),T_ref, parsT);
  TC_Tah = tempcorr(C2K(Tah(:,1)), T_ref, [T_A, T_L, T_AL]);
  TC_tL7 = tempcorr(temp.tL7_f1, T_ref, parsT);
  TC_tL19 = tempcorr(temp.tL19_f1, T_ref, parsT);
  TC_tL15 = tempcorr(temp.tL15_f1, T_ref, parsT);
  TC_DOJO = tempcorr(temp.DOJO, T_ref, parsT);
  TC_DOah = tempcorr(temp.DOah, T_ref, parsT);
  TC_DOS1 = tempcorr(temp.DOS1, T_ref, parsT);
  TC_DOS2 = tempcorr(temp.DOS1, T_ref, parsT);
  TC_DOS3 = tempcorr(temp.DOS1, T_ref, parsT);
  TC_DOS4 = tempcorr(temp.DOS1, T_ref, parsT);
  TC_DOS5 = tempcorr(temp.DOS1, T_ref, parsT);

  
  % zero-variate data

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];               % compose parameter vector
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f); % -, scaled times & lengths at f
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  % E_0 = J_E_Am * U_E0 * mu_E;        % J, initial reserve 
  Wd_0 = U_E0 * p_Am * w_E/ mu_E; % g, in initial dry weight

  % hatch   
  [U_H aUL] = ode45(@dget_aul, [0; U_Hh], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(end,1);             % d, age at hatch at f and T
  aT_h = a_h/ TC_ah;             % d, age at hatch at f and T
  L_h = aUL(end,3);                  % cm, structural length at f
  Lw_h = L_h/ del_M;               % cm, S-V length at hatch at f

  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_Me;               % cm, length at birth at f
  % Lw_bt = L_b/ del_T;               % cm, total length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
  % aT_b = tau_b/ k_M/ TC_ab;         % d, age at birth at f and T
  a_b = tau_b/ k_M;         % d, age at birth at f and T
  
  t_hb = a_b - a_h;
  tT_hb = t_hb /TC_hb;

  % % metam
  % tT_j = (tau_j - tau_b)/ k_M/ TC_aj;      % d, time since birth at metam
  a_j = ((tau_j/k_M) - a_h)/ TC_aj;      % d, time since birth at metam
  L_j = L_m * l_j;                 % cm, structural length at metam
  Lw_j = L_j/ del_M;               % cm, total length at metam
  s_M = l_j/ l_b;                  % -, acceleration factor

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, standard length at puberty at f
  Ww_p = L_p^3 * (1 + f * w);       % g, wet weight at puberty 
  aT_p = tau_p/ k_M /TC_ap;         % d, age at puberty at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate standard length at f
  % Lw_it = L_i/ del_T;               % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);             % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % different f
  pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp];

  % f_1
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_1); % -, scaled times & lengths at f
  % puberty 
  L_p_1 = L_m * l_p;                  % cm, structural length at puberty at f
  Ww_pf_1 = L_p_1^3 * (1 + f_1 * w);       % g, wet weight at puberty 
  % s_M = l_j/l_b; 
  % ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0];
  % [t, ELHR] = ode45(@dget_ELHR_abj, [0, 10*7], ELHR_0, [], pars_abj, TC_LN, f_1, s_M);        
  % EN_f1 = ELHR(end,4)/E_0;
  %
  % f_2
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_2); % -, scaled times & lengths at f
  % puberty 
  L_p_2 = L_m * l_p;                  % cm, structural length at puberty at f
  Ww_pf_2 = L_p_2^3 * (1 + f_2 * w);       % g, wet weight at puberty 
  % s_M = l_j/l_b; 
  % ELHR_0 = [E_m*L_p_2^3, L_p_2, E_Hp, 0];
  % [t, ELHR] = ode45(@dget_ELHR_abj, [0, 10*7], ELHR_0, [], pars_abj, TC_LN, f_2, s_M);
  % EN_f2 = ELHR(end,4)/E_0;
  %
  % f_3
  % [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_3); % -, scaled times & lengths at f
  % puberty 
  % L_p_3 = L_m * l_p;                  % cm, structural length at puberty at f
  % Ww_pf_3 = L_p_3^3 * (1 + f_3 * w);       % g, wet weight at puberty 
  % s_M = l_j/l_b; 
  % ELHR_0 = [E_m*L_p_3^3, L_p_3, E_Hp, 0];
  % [t, ELHR] = ode45(@dget_ELHR_abj, [0, 10*7], ELHR_0, [], pars_abj, TC_LN, f_3, s_M);        
  % EN_f3 = ELHR(end,4)/E_0;                  
  %

  % males
  p_Am_m = z_m * p_M/ kap;       % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;             % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ kap/ E_m_m ;        % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m / E_G;  % mol/mol, reserve capacity
  L_m_m = v/ k_M/ g_m;           % cm, maximum length

  % J_E_Am_m = p_Am_m/ mu_E;        % mol/d.cm^2, {J_EAm}, max surface-spec assimilation flux
  % ome     = m_Em * w_E * p.d_V/ p.d_E/ w_V; % -, \omega, contribution of ash free dry mass of reserve to total ash free dry biomass
  % w       = ome;                   % -, just for consistency with the past

  U_Hb = E_Hb/ p_Am_m;                  % cm^2 d, scaled maturity at birth
  V_Hb = U_Hb/ (1 - kap);               % cm^2 d, scaled maturity at birth
  v_Hb = V_Hb * g_m^2 * k_M^3/ v^2;     % -, scaled maturity density at birth

  U_Hj = E_Hj/ p_Am_m;                  % cm^2 d, scaled maturity at metam
  V_Hj = U_Hj/ (1 - kap);               % cm^2 d, scaled maturity at metam
  v_Hj = V_Hj * g_m^2 * k_M^3/ v^2;     % -, scaled maturity density metam

  U_Hp = E_Hp_m/ p_Am_m;                  % cm^2 d, scaled maturity at puberty
  V_Hp = U_Hp/ (1 - kap);               % cm^2 d, scaled maturity at puberty
  v_Hp = V_Hp * g_m^2 * k_M^3/ v^2;     % -, scaled maturity density at puberty

  pars_tj_m = [g_m; k; l_T; v_Hb; v_Hj; v_Hp];               % compose parameter vector
  [~, ~, ~, ~, l_p_m, ~, l_i_m, ~, ~, info] = get_tj(pars_tj_m, f); % -, scaled times & lengths at f
  
  L_p_m = L_m_m * l_p_m;                  % cm, structural length at puberty at f
  Lw_p_m = L_p_m/ del_M;                % cm, standard length at puberty at f
  
  L_i_m = L_m_m * l_i_m; 
  Lw_i_m = L_i_m/ del_M;

  Ww_p_m = L_p_m^3 * (1 + f * w);       % g, wet weight at puberty 

  % pack to output
  prdData.ah = aT_h;
  prdData.aj = a_j;
  prdData.thb = tT_hb;
  % prdData.ab = aT_b;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lh = Lw_h;
  prdData.Lj = Lw_j;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  % prdData.Lit = Lw_it;
  prdData.Wd0 = Wd_0;
  prdData.Wwb = Ww_b;
  prdData.Wwp = Ww_p;
  prdData.Wwp_m = Ww_p_m;
  prdData.Wwp_f1 = Ww_pf_1;
  prdData.Wwp_f2 = Ww_pf_2;
  % prdData.Wwp_f3 = Ww_pf_3;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  % prdData.N_f1 = EN_f1;
  % prdData.N_f2 = EN_f2;
  % prdData.N_f3 = EN_f3;
  prdData.Lp_m = Lw_p_m;
  prdData.Li_m = Lw_i_m;

  % uni-variate data
  
  % time-length
  % WrigWoot2004
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL_WrigWoot2004);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL;      

  kT_M = k_M * TC_tL; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  L_b = L_m * l_b;
  L_j = L_m * l_j;
  L_p = L_m * l_p;
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL_WrigWoot2004((tL_WrigWoot2004(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL_WrigWoot2004((tL_WrigWoot2004(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_WrigWoot2004 = L./ del_M_bp; 

  % JoneHyne1950
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL_JoneHyne1950);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL;      

  kT_M = k_M * TC_tL; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  L_b = L_m * l_b; 
  L_j = L_m * l_j;
  L_p = L_m * l_p;
  L_i = L_m * l_i;

  tL_JoneHyne1950(:,1) = tL_JoneHyne1950(:,1)+t0;
  L_bj = L_b * exp(tL_JoneHyne1950((tL_JoneHyne1950(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL_JoneHyne1950((tL_JoneHyne1950(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_JoneHyne1950 = L./ del_M_bp; 

  % Leloutre2016_1
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_Leloutre2016);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL;      

  kT_M = k_M * TC_tL; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  L_b = L_m * l_b; 
  L_j = L_m * l_j;
  L_p = L_m * l_p;
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL_Leloutre2016((tL_Leloutre2016(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL_Leloutre2016((tL_Leloutre2016(:,1) >= tT_j),1) - tT_j)); % cm, struc length
  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_Leloutre2016 = L./ del_M_bp; 

  % 
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  L_b = L_m * l_b; 

  % Leloutre2016 (tL7_f1)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL7_f1);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL7;      

  kT_M = k_M * TC_tL7; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b; 
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL7_f1((tL7_f1(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL7_f1((tL7_f1(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL7_f1((tL7_f1(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL7_f1 = L./ del_M_bp; 

  % Leloutre2016 (tL7_f2)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL7_f2);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL7;      

  kT_M = k_M * TC_tL7; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL7_f2((tL7_f2(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL7_f2((tL7_f2(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL7_f2((tL7_f2(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL7_f2 = L./ del_M_bp;

  % Leloutre2016 (tL19_f1)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL19_f1);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL19;      

  kT_M = k_M * TC_tL19; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL19_f1((tL19_f1(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL19_f1((tL19_f1(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL19_f1((tL19_f1(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL19_f1 = L./ del_M_bp;

  % Leloutre2016 (tL19_f2)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL19_f2);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL19;      

  kT_M = k_M * TC_tL19; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL19_f2((tL19_f2(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL19_f2((tL19_f2(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL19_f2((tL19_f2(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL19_f2 = L./ del_M_bp;

  % Leloutre2016 (tL15_f1)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL15_f1);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL15;      

  kT_M = k_M * TC_tL15; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL15_f1((tL15_f1(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL15_f1((tL15_f1(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL15_f1((tL15_f1(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL15_f1 = L./ del_M_bp;

  % Leloutre2016 (tL15_f2)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL15_f2);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL15;      

  kT_M = k_M * TC_tL15; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL15_f2((tL15_f2(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL15_f2((tL15_f2(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL15_f2((tL15_f2(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL15_f2 = L./ del_M_bp;

  % Leloutre2016 (tL15_f3)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL15_f3);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL15;      

  kT_M = k_M * TC_tL15; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL15_f3((tL15_f3(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL15_f3((tL15_f3(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL15_f3((tL15_f3(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL15_f3 = L./ del_M_bp;

  % Leloutre2016 (tL15_f4)
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_tL15_f4);
  tT_j = (tau_j - tau_b)/ k_M/ TC_tL15;      

  kT_M = k_M * TC_tL15; 
  rT_j = rho_j * kT_M; 
  rT_B = rho_B * kT_M;        
  
  % L_b = L_m * l_b;
  L_p = L_m * l_p;
  % L_j = L_m * l_j; 
  L_i = L_m * l_i;

  L_bj = L_b * exp(tL15_f4((tL15_f4(:,1) <= tT_j),1) * rT_j/ 3);  % cm, struc length
  % L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tL15_f4((tL15_f4(:,1) > tT_j),1) - tT_j)); % cm, struc length
  L_jm = L_i - (L_i - L_bj(end)) * exp( - rT_B * (tL15_f4((tL15_f4(:,1) >= tT_j),1) - tT_j)); % cm, struc length

  L = [L_bj; L_jm];
  w_b = (L_p - L)./(L_p-L_b); w_p = 1-w_b;
  del_M_bp = w_b.*del_Me + w_p .* del_M;
  ELw_tL15_f4 = L./ del_M_bp;

  
  % length-fecundity
  EN_L_Leloutre2016 = TC_LN * 365 * reprod_rate_j(LN_Leloutre2016(:,1)*del_M, f_Leloutre2016, pars_R); % #, fecundity (approximate spawning interval)
  EN_Lt_Wootton1973 = TC_LN * 365 * reprod_rate_j(LtN_Wootton1973(:,1)*del_T, f_LtN_Wootton1973, pars_R); % #, fecundity (approximate spawning interval)

  % respiration 
  pars_p = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose pars
  p_ref = p_Am * L_m^2; % J/d, max assimilation power at max size
  % 
  % % Andrea2017 20°C 
  % [~, ~, ~, l_j20, l_p20, l_b20, l_i20, ~, ~, info20] = get_tj(pars_tj, f_WwJO_Andrea20); % -, scaled times & lengths at f
  % % 
  % L_20 = (WwJO_Andrea20(:,1)/ (1 + f_WwJO_Andrea20 * w)) .^ (1/3);             % cm, structural length
  % l_20 = L_20/ L_m; 
  % sM_20 = min(l_j20, l_20)/l_b20; 
  % % scaled powers: powers per max assimilation {p_Am} L_m^2
  % [pACSJGRD_Andrea20, info] = scaled_power_j(L_20, f_WwJO_Andrea20, pars_p, l_b20, l_j20, l_p20);  % J/d, powers
  % if info == 0
  %   prdData = []; return
  % end
  % 
  % pACSJGRD_Andrea20(:,1) = 0 * pACSJGRD_Andrea20(:,1);
  % pACSJGRD_Andrea20 = p_ref * pACSJGRD_Andrea20;
  % 
  % J_M_Andrea20 = - (n_M\n_O) * eta_O * pACSJGRD_Andrea20(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  % EJT_O_Andrea20 = - J_M_Andrea20(3,:)' * TC_JO_20 * 32 ./WwJO_Andrea20(:,1);      % g O2/g/d, O2 consumption  
  % 
  % Brafield1976 10°C 
  [~, ~, ~, l_j10, l_p10, l_b10, l_i10, ~, ~, info10] = get_tj(pars_tj, f_WwJO_Brafield10); % -, scaled times & lengths at f
  % 
  L_10 = (WwJO_Brafield10(:,1)/ (1 + f_WwJO_Brafield10 * w)) .^ (1/3);             % cm, structural length
  l_10 = L_10/ L_m; 
  sM_10 = min(l_j10, l_10)/l_b10; 
  % scaled powers: powers per max assimilation {p_Am} L_m^2
  [pACSJGRD_Brafield10, info] = scaled_power_j(L_10, f_WwJO_Brafield10, pars_p, l_b10, l_j10, l_p10);  % J/d, powers
  if info == 0
    prdData = []; return
  end

  pACSJGRD_Brafield10(:,1) = 0 * pACSJGRD_Brafield10(:,1);
  pACSJGRD_Brafield10 = p_ref * pACSJGRD_Brafield10;
  
  J_M_Brafield10 = - (n_M\n_O) * eta_O * pACSJGRD_Brafield10(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJT_O_Brafield10 = - J_M_Brafield10(3,:)' * TC_JO_10 * 32 ./WwJO_Brafield10(:,1);      % g O2/g/d, O2 consumption  
  

  % Temperature dependant respiration rate 
  [~, ~, ~, l_j, l_p, l_b, l_i, ~, ~, info] = get_tj(pars_tj, f_TJO); % -, scaled times & lengths at f
  % 
  Ww = 1.41; Lw = 5.03;
  L = Lw * del_M;             % cm, structural length
  % scaled powers: powers per max assimilation {p_Am} L_m^2
  [pACSJGRD_TJO, info] = scaled_power_j(L, f_TJO, pars_p, l_b, l_j, l_p);  % J/d, powers
  if info == 0
    prdData = []; return
  end

  pACSJGRD_TJO(:,1) = 0 * pACSJGRD_TJO(:,1);
  pACSJGRD_TJO = p_ref * pACSJGRD_TJO;
  
  J_M_TJO = - (n_M\n_O) * eta_O * pACSJGRD_TJO(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJT_O = - J_M_TJO(3,:)' .* TC_TJO * 32e6 /24/Ww;      % mg O2/kg/h, O2 consumption  


  % DO saturation dependant respiration rate 
  [~, ~, ~, l_j, l_p, l_b, l_i, ~, ~, info] = get_tj(pars_tj, f_DOJO); % -, scaled times & lengths at f
  % 
  Ww = 1.41; Lw = 5.03;
  L = Lw * del_M;             % cm, structural length
  % scaled powers: powers per max assimilation {p_Am} L_m^2
  [pACSJGRD_DOJO, info] = scaled_power_j(L, f_DOJO, pars_p, l_b, l_j, l_p);  % J/d, powers
  if info == 0
    prdData = []; return
  end

  pACSJGRD_DOJO(:,1) = 0 * pACSJGRD_DOJO(:,1);
  pACSJGRD_DOJO = p_ref * pACSJGRD_DOJO;

  if DOJO(:,1) <= DOcrit
      DOC_DOJO = 0;
  else
      DOC_DOJO = DOa .* ((DOJO(:,1) - DOcrit)./(DOJO(:,1) - DOcrit + DOb)); 
  end

  % parsDO=[DOcrit, DOa, DOb];
  % DOC_DOJO = docorr(DOJO(:,1),parsDO)';

  J_M_DOJO = - (n_M\n_O) * eta_O * pACSJGRD_DOJO(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJDO_O = - J_M_DOJO(3,:)' * TC_DOJO .* DOC_DOJO * 32e6 /24/Ww;      % mg O2/kg/h, O2 consumption  

  % feeding
  % % Hovel2015, 23°C
  % [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_WwJX); % -, scaled times & lengths at f
  % s_M = l_j/l_b;
  % L_X = (WwJX(:,1)/ (1 + f_WwJX * w)) .^ (1/3);             % cm, structural length
  % EJ_X = p_Am/kap_X * s_M * f_WwJX * L_X.^2 * (w_X/mu_X)* TC_JX ./WwJX(:,1); % g/d/d, ingestion rate
  % 
  % L_Cmax  = (1.7/ (1 + f_WwJX * w)) .^ (1/3);      
  % ECmax = p_Am * s_M * f_WwJX * L_Cmax^2 .* TC_Cmax ;  % g/d, assimilation rate
  % ECmax = ECmax./max(ECmax) ;  % g/d, assimilation rate
  
  % Wootton1980
  [tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_WwJX_Wootton); % -, scaled times & lengths at f
  s_M = l_j/l_b;
  % 19°C 
  L_X = (WwJX_Wootton19(:,1)/ (1 + f_WwJX_Wootton * w)) .^ (1/3);             % cm, structural length
  EJ_X_19 = p_Am/kap_X * s_M * f_WwJX_Wootton * L_X.^2 * (w_X/mu_X)* TC_JX_19; % g/d, ingestion rate
  % 12.5°C
  L_X = (WwJX_Wootton12(:,1)/ (1 + f_WwJX_Wootton * w)) .^ (1/3);             % cm, structural length
  EJ_X_12 = p_Am/kap_X * s_M * f_WwJX_Wootton * L_X.^2 * (w_X/mu_X)* TC_JX_12; % g/d, ingestion rate
  % 7°C 
  L_X = (WwJX_Wootton7(:,1)/ (1 + f_WwJX_Wootton * w)) .^ (1/3);             % cm, structural length
  EJ_X_7 = p_Am/kap_X * s_M * f_WwJX_Wootton * L_X.^2 * (w_X/mu_X)* TC_JX_7; % g/d, ingestion rate

  
  % temperature - age at hatching
  Ea_h = a_h./ TC_Tah;
  
  % DO - age at hatching
  if DOah(:,1) <= DOcrit
      DOC_DOah = 0;
  else
      DOC_DOah = DOa .* ((DOah(:,1) - DOcrit)./(DOah(:,1) - DOcrit + DOb)); 
  end

  % parsDO=[DOcrit, DOa, DOb];
  % DOC_DOah = docorr(DOah(:,1),parsDO)';
  EDOa_h = a_h./ TC_DOah ./ DOC_DOah;

  % time - embryo survival: Eq (6.2) of Kooy2010
  pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp, memb];
  ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0, 1];
  tspan = [0; DOS1(:,1)];
  [t, ELHR] = ode45(@dget_ELHR_abj, tspan, ELHR_0, [], pars_abj, TC_DOS1, f, s_M);  
  
  for DO = [100 80 60 40 20]
   if DO <= DOcrit
      DOC_DOS(DO) = 0;
   else
      DOC_DOS(DO) = (K * (DO - DOcrit)) / (1 + K * (DO - DOcrit));
   end
  end

  % DO=100
  EDOS1 = ELHR(:,5) * DOC_DOS(100); 
  EDOS1(1) =[];
  
  % DO=80
  % pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp, mu_emb80];
  % ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0, 1];
  % tspan = [0; DOS2(:,1)];
  % [t, ELHR] = ode45(@dget_ELHR_abj, tspan, ELHR_0, [], pars_abj, TC_DOS2, f, s_M);  
  EDOS2 = ELHR(:,5)* DOC_DOS(80); 
  EDOS2(1) =[];
  
  %DO=60 
  % pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp, mu_emb60];
  % ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0, 1];
  % tspan = [0; DOS3(:,1)];
  % [t, ELHR] = ode45(@dget_ELHR_abj, tspan, ELHR_0, [], pars_abj, TC_DOS3, f, s_M);  
  EDOS3 = ELHR(:,5) * DOC_DOS(60); 
  EDOS3(1) =[];

  % DO=40
  % pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp, mu_emb40];
  % ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0, 1];
  % tspan = [0; DOS4(:,1)];
  % [t, ELHR] = ode45(@dget_ELHR_abj, tspan, ELHR_0, [], pars_abj, TC_DOS4, f, s_M);  
  EDOS4 = ELHR(:,5) * DOC_DOS(40); 
  EDOS4(1) =[];  
  
  % DO=20
  % pars_abj = [p_Am, v, p_M, k_J, kap, kap_G, E_G, E_Hb, E_Hj, E_Hp, mu_emb20];
  % ELHR_0 = [E_m*L_p_1^3, L_p_1, E_Hp, 0, 1];
  % tspan = [0; DOS5(:,1)];
  % [t, ELHR] = ode45(@dget_ELHR_abj, tspan, ELHR_0, [], pars_abj, TC_DOS5, f, s_M);  
  EDOS5 = ELHR(:,5) * DOC_DOS(20); 
  EDOS5(1) =[];

  % pack to output
  prdData.tL_WrigWoot2004 = ELw_WrigWoot2004;
  prdData.tL_JoneHyne1950 = ELw_JoneHyne1950;
  prdData.tL_Leloutre2016 = ELw_Leloutre2016;
  prdData.tL7_f1 = ELw_tL7_f1;
  prdData.tL7_f2 = ELw_tL7_f2;
  prdData.tL19_f1 = ELw_tL19_f1;
  prdData.tL19_f2 = ELw_tL19_f2;
  prdData.tL15_f1 = ELw_tL15_f1;
  prdData.tL15_f2 = ELw_tL15_f2;
  prdData.tL15_f3 = ELw_tL15_f3;
  prdData.tL15_f4 = ELw_tL15_f4;

  prdData.LN_Leloutre2016 = EN_L_Leloutre2016;
  prdData.LtN_Wootton1973 = EN_Lt_Wootton1973;
  % prdData.WwJO_Andrea20 = EJT_O_Andrea20;
  prdData.WwJO_Brafield10 = EJT_O_Brafield10;
  prdData.TJO = EJT_O;
  prdData.DOJO = EJDO_O;
  % prdData.WwJX = EJ_X;
  prdData.WwJX_Wootton19 = EJ_X_19;
  prdData.WwJX_Wootton12 = EJ_X_12;
  prdData.WwJX_Wootton7 = EJ_X_7;  
  % prdData.CmaxT = ECmax;
  prdData.Tah = Ea_h;
  prdData.DOah = EDOa_h;
  prdData.DOS1 = EDOS1;
  prdData.DOS2 = EDOS2;
  prdData.DOS3 = EDOS3;
  prdData.DOS4 = EDOS4;
  prdData.DOS5 = EDOS5;

 %% subfunctions
 function dELHR = dget_ELHR_abj(t, ELHR, p, tTC, tf, s_M)
  % Define changes in the state variables for abj model
  % t: time
  % ELHR: 4-vector with state variables
  %         E , J, reserve energy
  %         L , cm, structural length
  %         E_H , J , cumulated energy inversted into maturity (E_H in Kooijman 2010)
  %         E_R , J, reproduction buffer (E_R in Kooijman 2010)
  %         
  % dELHR: 4-vector with change in E, L, H, R

  % unpack state variables

  E = ELHR(1); L = ELHR(2); E_H = ELHR(3); S = min(1,ELHR(5));

  % unpack par
  p_Am = p(1); v = p(2); p_M = p(3); k_J = p(4); 
  kap = p(5); kap_G = p(6); 
  E_G = p(7); E_Hb = p(8); E_Hj = p(9); E_Hp = p(10);
  memb = p(11);

  % TC = spline1(t, tTC);  % C, temperature at t
  % f_t = spline1(t, tf);  % -, scaled functional response at t
  TC = tTC;
  f_t  = tf;

  % temp correction
  pT_Am = TC * p_Am ;
  vT = TC * v;  
  pT_M = TC * p_M;
  kT_J = TC * k_J; 

  % Fluxes
  pA = (pT_Am * s_M * f_t * L^2) * (E_H >= E_Hb);

  if  kap * E * s_M * vT >= pT_M * L^4 % section 4.1.5 comments to Kooy2010
    r = (E * s_M * vT/ L - pT_M * L^3/ kap)/ (E + E_G * L^3/ kap); % d^-1, specific growth rate  
  else 
    r = (E * s_M * vT/ L - pT_M * L^3/ kap)/ (E + kap_G * E_G * L^3/ kap); % d^-1, specific growth rate                                      
  end
  pC  = E * (s_M * vT/ L - r); % J/d, mobilized energy flux
  % generate derivatives
  dE    = pA - pC;  % J/d, change in energy in reserve
  dL    = r * L / 3;    % cm^3/d, change in structural volume
  dE_H  = ((1 - kap) * pC - kT_J * E_H) * (E_H<E_Hp);     % J/d, change in cumulated energy invested in maturation
  dE_R  = ((1 - kap) * pC - kT_J * E_Hp) * (E_H >= E_Hp); % J/d, change in reproduction buffer
  dS = - S * memb; % survival rate embryos

  % pack derivatives
  dELHR = [dE; dL; dE_H; dE_R; dS];

