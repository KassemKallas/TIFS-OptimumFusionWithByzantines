clear all;
clc;
%======================== T 4 N 16 ======================================

%-------------------------------------------------------------------------
sim_param.PH1 = 0.5;



sim_param.T = 4;

sim_param.N = 20;
sim_param.K1 = 11;
sim_param.alfa = (sim_param.N-sim_param.K1)/sim_param.N;

sim_param.epsilon = 0.1;

sim_param.Nprove = 50000;
sim_param.delta= 1 - sim_param.epsilon;

sim_param.possible_system_states =  dec2bin(0:2^sim_param.T -1,sim_param.T)-'0';
sim_param.L = sim_param.N/2;



%=====================Initialize Varshney and LLR=========================


        sim_param.gammas = 0:sim_param.T;
        sim_param.Pd_Hp = 1-sim_param.epsilon;%detection at honest
        sim_param.Pfa_Hp = sim_param.epsilon;%fa at honest
        sim_param.Pd_Bp = 1-sim_param.epsilon;%detection at Byzantine
        sim_param.Pfa_Bp = sim_param.epsilon;%false alarm at byzantines
        sim_param.Nsoglie_LLR =  10;
        
%=====================END Initialize Varshney and LLR=========================
        

fprintf('Alfa = %f\n',sim_param.alfa);
fprintf('T = %f\n',sim_param.T);
fprintf('N = %f\n',sim_param.N);
fprintf('K honests = %f\n',sim_param.K1);
fprintf('M byzantines = %f\n',(sim_param.N-sim_param.K1));
fprintf('epsilon = %f\n',sim_param.epsilon);

xx=0;
for pmal_dec = 1:0.01:6
xx=xx+1;    
sim_param.Pmal = (pmal_dec+4)/10;
fprintf('Pmal = %f\n',sim_param.Pmal);


[results] = function_independent_states(sim_param);
%============= Results Varshney and LLR====================================
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

            % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
            
            Varsh(xx)= min(PFAr + 1-PDr);
            SOFT(xx) =min(PFAr_LLR+1-PDr_LLR);
           Majority(xx)= results.error_majority;

%============= End Results Varshney and LLR====================================

        Independent(xx,:) =  results.error_eq4;


end
    

save('Independent.mat','Independent');
save('Varsh.mat','Varsh');
save('SOFT.mat','SOFT');
save('Majority.mat','Majority');



