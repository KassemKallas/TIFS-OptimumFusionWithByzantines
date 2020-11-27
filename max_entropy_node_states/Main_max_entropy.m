
sim_param.PH1 = 0.5;

sim_param.N = 20; % TOTALE NODI
% sim_param.K = 9; % ONESTI
% sim_param.M = sim_param.N-sim_param.K; % DISONESTI
% sim_param.alfa = (sim_param.N-sim_param.K)/sim_param.N;
sim_param.T = 4;

sim_param.L = sim_param.N/2;


sim_param.possible_system_states =  dec2bin(0:2^sim_param.T -1,sim_param.T)-'0';

sim_param.epsilon = 0.1;
sim_param.delta= 1 - sim_param.epsilon;


%===============Initialization of Varshney and LLR===================
sim_param.gammas = 0:sim_param.T;
sim_param.Pd_Hp = 1-sim_param.epsilon;%detection at honest
sim_param.Pfa_Hp = sim_param.epsilon;%fa at honest
sim_param.Pd_Bp = 1-sim_param.epsilon;%detection at Byzantine

sim_param.Pfa_Bp = sim_param.epsilon;%false alarm at byzantine
sim_param.Nsoglie_LLR =  10;
%==============End Initialization of Varshney and LLR===============

sim_param.Nprove = 50000;

fprintf('N = %f\n',sim_param.N);
% fprintf('K = %f\n',sim_param.K);
% fprintf('M = %f\n',sim_param.M);
% fprintf('alfa = %f\n',sim_param.alfa);

fprintf('T = %f\n',sim_param.T);

for pmal_dec = 1:1:6
    
    sim_param.Pmal = (pmal_dec+4)/10;
    
    fprintf('Pmal Real = %f\n',sim_param.Pmal);
    
    %[results] = function_fixed_states(sim_param);
    [results] = function_max_entropy_optimized(sim_param);
    
    %==========RESULTS COMPUTATIONS OF VARSHNEY AND LLR========================================
    PD_IDB = results.PD_IDB;
    PFA_IDB = results.PFA_IDB;
    P_ISO_H= results.P_ISO_H; %ISO VALUE OH HONEST AT VARSHNEY SCHEME
    P_ISO_B= results.P_ISO_B;%ISO VALUE OF BYZANTINES AT VARSHNEY
    PD_IDB_LLR = results.PD_IDB_LLR;
    PFA_IDB_LLR = results.PFA_IDB_LLR;
    P_ISO_H_LLR = results.P_ISO_H_LLR;
    P_ISO_B_LLR = results.P_ISO_B_LLR;
    PFA = results.PFA;
    PD = results.PD;
    PFAr = results.PFAr; % Varshney
    PDr = results.PDr;
    PFAr_LLR = results.PFAr_LLR; % LLR
    PDr_LLR = results.PDr_LLR;
    
    Hard_N_13_T_4(pmal_dec) = min(PFAr + 1-PDr);
    Soft_N_13_T_4(pmal_dec) = min(PFAr_LLR+1-PDr_LLR);
    %==========END OF RESULTS COMPUTATIONS OF VARSHNEY AND LLR========================================
    perr_majority(pmal_dec) = results.error_majority;
%         
%     MB_max_entropy_N_16_T_4(pmal_dec,:)=  results.p_err_MB; % this is only for max entropy game result MB
     
max_entropy_N_13_T_4(pmal_dec,:)=  results.p_err; % this is only for max entropy game result KK


    
end

 %save('MB_max_entropy_N_16_T_4.mat','MB_max_entropy_N_16_T_4');
 save('max_entropy_N_13_T_4.mat','max_entropy_N_13_T_4');

 
% M_all = results.M_all;
 %R_mat_array = results.R_mat_array;
 
% SS = results.SS;
 
 
 