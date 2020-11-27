function[results] = function_fixed_states(sim_param)

Nprove = sim_param.Nprove;
PH1 = sim_param.PH1;
PH0 = 1-PH1;
L=sim_param.L;

Pmal = sim_param.Pmal;
%================ Varsh and LLR===================================
Pd_Hp = sim_param.Pd_Hp;
Pfa_Hp = sim_param.Pfa_Hp;
Pd_Bp = sim_param.Pd_Bp;
Pfa_Bp = sim_param.Pfa_Bp;
gammas = sim_param.gammas;


Pd_H = Pd_Hp;
Pfa_H = Pfa_Hp;
Pd_B = Pmal*(1-Pd_Bp)+(1-Pmal)*Pd_Bp;
Pfa_B = Pmal*(1-Pfa_Bp)+(1-Pmal)*Pfa_Bp;


Nsoglie_LLR = sim_param.Nsoglie_LLR;
Nerr_h = 0;
Nerr_b = 0;
N0=0;
N1=0;
Nerr_hr = zeros(length(gammas),1);
Nerr_br = zeros(length(gammas),1);
Nerr_hr_LLR = zeros(Nsoglie_LLR,1);
Nerr_br_LLR = zeros(Nsoglie_LLR,1);
Nerr_H = zeros(length(gammas),1);
Nerr_B = zeros(length(gammas),1);
Nerr_H_LLR = zeros(Nsoglie_LLR,1);
Nerr_B_LLR = zeros(Nsoglie_LLR,1);
%=================================================================


N = sim_param.N;
M = sim_param.M;
K = sim_param.K;
alfa = sim_param.alfa;
T = sim_param.T;

possible_states = sim_param.possible_system_states;



epsilon = sim_param.epsilon;
delta = sim_param.delta;

delta_Byz = (1-delta)*(1-Pmal) + (delta)*Pmal;





for np = 1:Nprove
    if rem(np,500) == 0
        fprintf('Simulation %d su %d\n',np,Nprove);
    end;
    rd = rand(1,T);
    P = zeros(1,T);
    P(rd < PH1) = 1; % generating the system state
    
    
    UH = zeros(K,T); % to save the decisions at honests
    UB = zeros(M,T); % to save the decisions at Byzantines
    D = zeros(1,T); % here this is used to save the decisions of the majority rule.
    LLRs_OUT = zeros(N,T); % this belongs to state of the art schemes
    for t = 1:T
        if P(t) == 1 % if the system state was 1
            UH(:,t) = 1; % the reports of honest were first equal to one
            GH = rand(K,1);
            UH(GH < epsilon,t) = 0; % switch there report for honest
            UB(:,t) = 1;
            GB = rand(M,1);
            UB(GB < delta_Byz,t) = 0;
        else % if the system state was zero
            GH = rand(K,1);
            UH(GH < epsilon,t) = 1;
            GB = rand(M,1);
            UB(GB < delta_Byz,t) = 1;
        end;
        U_ALL = [UB(:,t);UH(:,t)];
        
        %=========================  SETUP VARSHNEY AND LLR ===========================
        Prob_err = Pmal*M/N;
        
        U_ALL_2 = [UH(:,t);UB(:,t)];
        Num_ones = length(find(U_ALL_2 == 1));
        Num_zeros = length(find(U_ALL_2 == 0));
        
        
        
        alpha = M/N;
        P1 = (1-alpha)*Pfa_H+alpha*Pfa_B;
        P2 = (1-alpha)*Pd_H+alpha*Pd_B;
        PUH0 = ((1-Prob_err)*P1+Prob_err*(1-P1))^Num_ones * ((1-Prob_err)*(1-P1)+Prob_err*P1)^Num_zeros;
        PUH1 = ((1-Prob_err)*P2+Prob_err*(1-P2))^Num_ones * ((1-Prob_err)*(1-P2)+Prob_err*P2)^Num_zeros;
        for dec = 1:N
            if U_ALL_2(dec) == 0
                PUH0d = PUH0/(((1-Prob_err)*(1-P1)+Prob_err*P1));
                PUH1d = PUH1/((1-Prob_err)*(1-P2)+Prob_err*P2);
                Px0U = (1-Prob_err)*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = Prob_err*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            else
                PUH0d = PUH0/((1-Prob_err)*P1+Prob_err*(1-P1));
                PUH1d = PUH1/((1-Prob_err)*P2+Prob_err*(1-P2));
                Px0U = Prob_err*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = (1-Prob_err)*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            end;
            LLRs_OUT(dec,t) = abs(log(Px0U/Px1U));
        end;
        %=========================END SETUP VARSHNEY AND LLR ==========================================
        
        
        
        if sum([UH(:,t);UB(:,t)]) >= L % obtain the majority rule result
            D(t) = 1;
        else
            D(t) = 0;
        end;
        
        R_matrix(:,t) = U_ALL;
        
    end
    
    
    %=====================Decoding Using Varshney and LLR=======================================
    REL = sum(LLRs_OUT,2);
    if  (max(REL)- min(REL)) > 0 %Nsoglie_LLR > 0
        
        for i=1:Nsoglie_LLR
            SOGLIA_LLRs(i) = min(REL) + ((max(REL)-min(REL))/Nsoglie_LLR)*i;
        end
    else
        SOGLIA_LLRs = mean(REL)*ones(1,Nsoglie_LLR);
    end;
    %disp(SOGLIA_LLRs);
    for is = 1:Nsoglie_LLR
        SOGLIA_LLR = SOGLIA_LLRs(is);
        Nerr_H_LLR(is) = Nerr_H_LLR(is) + length(find(REL(1:K) < SOGLIA_LLR));
        Nerr_B_LLR(is) = Nerr_B_LLR(is) + length(find(REL(K+1:N) >= SOGLIA_LLR));
    end;
    
    
    Dall_H = repmat(D,K,1);
    Errs_H = xor(Dall_H,UH);
    Dall_B = repmat(D,M,1);
    Errs_B = xor(Dall_B,UB);
    eta_H = sum(Errs_H,2);
    eta_B = sum(Errs_B,2);
    for ig = 1:length(gammas)
        gamma = gammas(ig);
        Nerr_H(ig) = Nerr_H(ig) + length(find(eta_H > gamma));
        Nerr_B(ig) = Nerr_B(ig) + length(find(eta_B <= gamma));
    end;
    
    
    
    indx = find(P == 1);
    Nerr_h = Nerr_h + length(find(D(indx) == 0));
    indx = find(P == 0);
    Nerr_b = Nerr_b + length(find(D(indx) == 1));
    % variabili che servono nel caso non simmetrico
    N0 = N0 + numel(find(P==0));
    N1 = N1 + numel(find(P==1));
    
    %Valutazione delle prestazioni dopo rimozione
    for ig = 1:length(gammas)
        Dr = 0*D;
        gamma = gammas(ig);
        indxH = find(eta_H <= gamma);
        indxB = find(eta_B <= gamma);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;
        for t = 1:T
            %Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))/2
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;
        end
        indx = find(P == 1);
        Nerr_hr(ig) = Nerr_hr(ig) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br(ig) = Nerr_br(ig) + length(find(Dr(indx) == 1));
    end;
    %Valutazione delle prestazioni dopo rimozione nel caso LLR
    for is = 1:Nsoglie_LLR
        Dr = 0*D;
        SOGLIA_LLR = SOGLIA_LLRs(is);
        indxH = find(REL(1:K) >= SOGLIA_LLR);
        indxB = find(REL(K+1:N) >= SOGLIA_LLR);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;
        for t = 1:T
            %%%%Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))/2
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;
        end
        indx = find(P == 1);
        Nerr_hr_LLR(is) = Nerr_hr_LLR(is) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br_LLR(is) = Nerr_br_LLR(is) + length(find(Dr(indx) == 1));
    end;
    %============================END DECODING VARSHNEY AND LLR==============================================
    
    
    
    
    
    %===================     decoding max entropy         ==================================================
    all_idx_combinations = permute_matrix(N,M);
    for i=1:2^T
        
        state = possible_states(i,:);
        
%         for x=1:size(all_idx_combinations,1)
%             
%             current_permutation = all_idx_combinations(x,:);
%             
%             for xx=1:length(current_permutation)
%                 permuted_matrix(xx,:) = R_matrix(current_permutation(xx),:);
%             end
%             
%             
%             %========= From here to compute T-neq and neq========================
%             nonzeros = bsxfun(@minus,permuted_matrix,state);
%             
%             neq_state = T - sum(nonzeros~=0,2);
%             
%             Tminusneq = bsxfun(@minus,T,neq_state);
%             
%             list_H = 1:M;
%             M_Hon = permuted_matrix;   % obtain a matrix with honest part only
%             M_Hon(list_H,:)=[];
%             
%             list_B = M+1:N;
%             M_Byz = permuted_matrix;   % obtain a matrix with byzantine part only
%             M_Byz(list_B,:)=[];
%             
%             
%             
%             nonzeros_H = bsxfun(@minus,M_Hon,state);
%             
%             neq_H = T - sum(nonzeros_H~=0,2);   %values for honest part
%             
%             Tminusneq_H = bsxfun(@minus,T,neq_H);
%             
%             
%             
%             nonzeros_B = bsxfun(@minus,M_Byz,state);
%             
%             neq_B = T - sum(nonzeros_B~=0,2);
%             
%             Tminusneq_B = bsxfun(@minus,T,neq_B);
%             
%             %========= Finish computing T-neq and neq========================
%             
%             for Pmal_guess_dec = 1:1:6
%                 
%                 Pmal_guess = (Pmal_guess_dec + 4)/10;
%                 
%                 delta_Byz_guess = (1-delta)*(1-Pmal_guess) + (delta)*Pmal_guess;
%                 
%                 fB(x,Pmal_guess_dec)= prod(((1-delta_Byz_guess).^neq_B).*(delta_Byz_guess.^Tminusneq_B));
%                 
%                 fH(x)= prod(((1-epsilon).^neq_H).*(epsilon.^Tminusneq_H));
%                 
%                 factor(x, Pmal_guess_dec) = fB(x,Pmal_guess_dec)*fH(x);
%                 
%             end
%             
%         end
%         states_score(i,:)=log(sum(factor,1));
        
        nonz = bsxfun(@minus,R_matrix,state);
        neq = T - sum(nonz~=0,2);
        
        for Pmal_guess_dec = 1:1:6
            Pmal_guess = (Pmal_guess_dec + 4)/10;
            delta_Byz_guess = (1-delta)*(1-Pmal_guess) + (delta)*Pmal_guess;
            % states_score(i,Pmal_guess_dec) = log(fnk(M,N,N,epsilon, delta_Byz_guess, T,neq,1)); %recursive implementation of the function
            states_score(i,Pmal_guess_dec) = FNKmatrix(M, N, epsilon, delta_Byz_guess, T, neq); % iterative implementation (convenient)
        end
        
    end
    
    
    %========================= End of Decoding Max Entropy Case===========================================================
    
    
    
    %=========================Final Decision Max entropy =================================================================
    [val_dec, idx_dec] = max(states_score,[],1); % max and index over each column
    
    for iii = 1:length(val_dec)
        dec_fusion(iii,:) = possible_states(idx_dec(iii),:); % each row in this matrix corresponds to the decision at FC with different guessed Pmal
        % dec_fusion matrix should be of size 6xT
        
        err_dec(np,iii) = numel(find(P~=dec_fusion(iii,:))); % computing the number of decision errors for each pmal guess over all the runs
        
    end
    %=========================Final Decision Max entropy =================================================================
    
    total_nb_trials(np) = length(P);
    
    diff_majority = P-D;
    
    err_nb_majority (np) = nnz(diff_majority); % number of error decision for majority rule
    
end

p_err(1,:) = sum(err_dec,1)/sum(total_nb_trials); % this row contains the results for different pmal guess under one real pmal then, its size
%is 1x6


results.p_err = p_err;


%====================Decision and Error Computation of Varshney and LLR =========================
results.PD = 1-Nerr_h/Nprove/T;
results.PFA = Nerr_b/Nprove/T;
for ig = 1:length(gammas)
    results.PDr(ig) = 1-Nerr_hr(ig)/Nprove/T;
    results.PFAr(ig) = Nerr_br(ig)/Nprove/T;
    results.PD_IDB(ig) = 1-Nerr_H(ig)/K/Nprove;
    results.PFA_IDB(ig) = Nerr_B(ig)/M/Nprove;
    results.P_ISO_H(ig) = Nerr_H(ig)/K/Nprove;
    results.P_ISO_B(ig) = 1 - Nerr_B(ig)/M/Nprove;
end;
for is = 1:Nsoglie_LLR
    results.PD_IDB_LLR(is) = 1-Nerr_H_LLR(is)/K/Nprove; % 1 - P_ISO^H
    results.PFA_IDB_LLR(is) = Nerr_B_LLR(is)/M/Nprove;  % 1 - P_ISO^B = P_NONISO^B
    results.P_ISO_H_LLR(is) = Nerr_H_LLR(is)/K/Nprove; % P_ISO^H
    results.P_ISO_B_LLR(is) = 1 - Nerr_B_LLR(is)/M/Nprove;  % P_ISO^B
    results.PDr_LLR(is) = 1-Nerr_hr_LLR(is)/Nprove/T;
    results.PFAr_LLR(is) = Nerr_br_LLR(is)/Nprove/T;
end;
%====================End Decision and Error Computation of Varshney and LLR ======================
error_majority = sum(err_nb_majority)/sum(total_nb_trials);

results.error_majority = error_majority;

end


