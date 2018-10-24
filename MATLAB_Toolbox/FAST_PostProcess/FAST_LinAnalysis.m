% FAST_LinAnalysis
% - Load and analyse a number of FAST linearization files

%% Setup Paths
Outdir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW';
OutfileBase = 'DTU_10MW_RWT';
nlin = 36;



%% Load Files
clear linout
for ind = 1:nlin
    % Load linearization
    linfile = [Outdir filesep OutfileBase '.' num2str(ind) '.lin'];
    linout(ind) = ReadFASTLinear(linfile);
    FileNameVec{ind} = [OutfileBase '.' num2str(ind) '.lin'];
end


%% Multiblade Coordinate Transform
% This takes the mixed rotating and non-rotating reference frame outputs
% from the FAST .lin files and transforms them into an entirely nonrotating
% reference frame state space. The idea here "normalizing" the outputs to
% be consistent, and thus enable linear analysis and control design.
%
% All nonrotating frame state space matrices are found using the methods
% that result in equations (29) through (33) of NREL_MBC.pdf (NREL report
% on MBC
%
% Note that this is done more extensively in mbc3.m, but mbc3.m needs some
% cleanup/updating for use with OpenFAST

for ind = 1:nlin
% some setup for cleaner code (poor memory management, though)
%     psi = linout(ind).Azimuth;
%     x_size = length(linout(ind).x_rotFrame)/2;
%     u_size = length(linout(ind).u_rotFrame);
%     y_size = length(linout(ind).y_rotFrame);
%     Omega = linout(ind).RotSpeed;
%     Omega_dot = 0;
%     A = linout(ind).A;
%     B = linout(ind).B;
%     C = linout(ind).C;
%         C1 = C(:,1:x_size);
%         C2 = C(:,x_size+1:end);
%     D = linout(ind).D;
    
% % Constructing setup matrices & variables
%     t1_til = [1 cos(psi) sin(psi); 1 cos(psi+2*pi/3) sin(psi+2*pi/3); 1 cos(psi+4*pi/3) sin(psi+4*pi/3)];
%     T1 = [];
%     for ind = 1:x_size/3
%         T1 = blkdiag(T1, t1_til);
%     end
%     T1 = T1 * diag(linout(ind).x_rotFrame(1:end/2));
%     T1 = T1 + eye(x_size) + -1.*diag(linout(ind).x_rotFrame(1:x_size));
%     
%     t2_til = [0 -sin(psi) cos(psi); 0 -sin(psi+2*pi/3) cos(psi+2*pi/3); 0 -sin(psi+4*pi/3) cos(psi+4*pi/3)];
%     T2 = [];
%     for ind = 1:x_size/3
%         T2 = blkdiag(T2, t2_til);
%     end
%     T2 = T2 * diag(linout(ind).x_rotFrame(1:x_size));
%     
%     t3_til = [0 -cos(psi) -sin(psi); 0 -cos(psi+2*pi/3) -sin(psi+2*pi/3); 0 -cos(psi+4*pi/3) -sin(psi+4*pi/3)];
%     T3 = [];
%     for ind = 1:x_size/3
%         T3 = blkdiag(T3, t2_til);
%     end
%     T3 = T3 * diag(linout(ind).x_rotFrame(1:x_size));
%     
%     T1c = [];
%     for ind = 1:u_size/3
%         T1c = blkdiag(T1c, t1_til);
%     end
%     T1c = T1c * diag(linout(ind).u_rotFrame);
%     T1c = T1c + eye(u_size) + -1.*diag(linout(ind).u_rotFrame);
%     
%     T1o = [];
%     for ind = 1:y_size/3
%         T1o = blkdiag(T1o, t1_til);
%     end
%     T1o = T1o * diag(linout(ind).y_rotFrame);
%     T1o = T1o + eye(y_size) + -1.*diag(linout(ind).y_rotFrame);

    FileNames = FileNameVec(ind)
    GetMats_f8
    mbc3
    
    linout(ind).A_MBC = MBC_A;
    linout(ind).B_MBC = MBC_B;
    linout(ind).C_MBC = MBC_C;
    linout(ind).D_MBC = MBC_D;
% % Nonrotating state matrices
%     linout(j).A_MBC = blkdiag(inv(T1), inv(T1)) * (A * [T1 zeros(size(T1)); Omega*T2 T1] - [Omega*T2 zeros(size(T2)); Omega^2*T3 + Omega_dot*T2 2*Omega*T2]);
%     linout(j).B_MBC = blkdiag(inv(T1), inv(T1)) * B * T1c;
%     linout(j).C_MBC = inv(T1o) * [C1*T1 + Omega*C2*T2 C2*T1];
%     linout(j).D_MBC = inv(T1o)*D*T1c;
end

%% Some analysis


for ind = 1:nlin

%     % Some organization
    A_MBC = linout(ind).A_MBC;
    B_MBC = linout(ind).B_MBC;
    C_MBC = linout(ind).C_MBC;
    D_MBC = linout(ind).D_MBC;

    sys = ss(A_MBC, B_MBC, C_MBC, D_MBC);
    
    linout(ind).A_MBC = A_MBC;
    linout(ind).B_MBC = B_MBC;
    linout(ind).C_MBC = C_MBC;
    linout(ind).D_MBC = D_MBC;
    
    % Find uncrontrollable modes
    lambdas = eig(A_MBC);
    lamrec = [];
    for k = 1:length(lambdas)
        ctr = PBH(A_MBC,B_MBC,lambdas(k));
        lamrec(k) = ctr;
    end
    linout(ind).lambdas = lambdas;
    linout(ind).lamPBHres = lamrec;
    
    % Remove unnecessary states and inputs
     rmstates = [5];                                    % 5 = ED Variable speed generator DOF
     rminputs = [1:7];                                  % 1-7 = All states besides GenTq and ColBldPitch
     A_MBC(rmstates,:) = []; A_MBC(:,rmstates) = [];
     B_MBC(rmstates,:) = []; 
     B_MBC(:,rminputs) = [];
     C_MBC(:,rmstates) = [];
     D_MBC(:,rminputs) = [];
     sys_rm = ss(A_MBC, B_MBC, C_MBC, D_MBC);
     linout(ind).x_desc_rm = linout(ind).x_desc;
     linout(ind).x_desc_rm(rmstates) = [];

    % Save removed states
    linout(ind).A_MBC_rm = A_MBC;
    linout(ind).B_MBC_rm = B_MBC;
    linout(ind).C_MBC_rm = C_MBC;
    linout(ind).D_MBC_rm = D_MBC;
    % Try to calculate the controllability Gramian and SVD
    try
        % Gramian
        linout(ind).Gc = gram(sys_rm,'c');
        % SVD
        [U,S,V] = svd(linout(ind).Gc);
        linout(ind).SVD.U = U;
        linout(ind).SVD.S = diag(S);
        linout(ind).SVD.V = V;
        
        % Processing 
        linout(ind).ctrldir = vecnorm(S*V');   % this is funky and probably wrong
        
%         % ~~THIS IS A WEIRD WAY TO LOOK AT THINGS~~
%         linout(i).SVD.Xwork = [];
%         xvec = zeros(length(linout(i).x_desc_rm),1);
%         for k = 1:length(linout(i).x_desc_rm)
%             xvec(k) = 1;
%             linout(i).SVD.Xwork = [linout(i).SVD.Xwork xvec'/linout(i).Gc*xvec];
%         end
%             
        
        %%
        % Need to look at U*S*V'*V to look at first ~6 singular values
    catch
        display(['Controllability Gramian could not be calculated for linearization point ', num2str(ind)]);
        
        % Fill the structure if its not controllable
        linout(ind).Gc = [];
        linout(ind).SVD.U = [];
        linout(ind).SVD.S = [];
        linout(ind).SVD.V = [];
        
        % Processing
        linout(ind).ctrldir = [];
        linout(ind).SVD.Xwork = [];
    end
    
    % Eigenvalues
    linout(ind).eigs = eig(A_MBC);
    
end

% contd = [];
% for j= 1:length(linout)
%     contd = [contd linout(j).contd']
% end
%% Some averaging analysis - CLEAN THIS UP AND SEPARATE IT INTO OTHER FUNCTIONS

Asum = [zeros(size(linout(1).A_MBC_rm))];
Bsum = [zeros(size(linout(1).B_MBC_rm))];
Csum = [zeros(size(linout(1).C_MBC_rm))];
Dsum = [zeros(size(linout(1).D_MBC_rm))];
for ind = 1:nlin
   Asum = Asum + linout(ind).A_MBC_rm;
   Bsum = Bsum + linout(ind).B_MBC_rm;
   Csum = Csum + linout(ind).C_MBC_rm;
   Dsum = Dsum + linout(ind).D_MBC_rm;
end

linavg.A = Asum./nlin;
linavg.B = Bsum./nlin;
linavg.C = Csum./nlin;
linavg.D = Dsum./nlin;

sys_avg = ss(linavg.A,linavg.B,linavg.C,linavg.D);
linavg.x_desc = linout(1).x_desc_rm;

% Gramian
linavg.Gc = gram(sys_avg,'c');
% SVD
[U,S,V] = svd(linavg.Gc);
linavg.SVD.U = U;
linavg.SVD.S = diag(S);
linavg.SVD.V = V;

% % This is a strange way to look at things 
% linavg.SVD.Xwork = [];
% xvec = zeros(length(linout(i).x_desc_rm),1);
% for k = 1:length(linout(i).x_desc_rm)
%     xvec(k) = 1;
%     linout(i).SVD.Xwork = [linout(i).SVD.Xwork xvec'/linout(i).Gc*xvec];
% end

%% Make some plots
figure(1), hold on
for ind = 1:nlin
    plot(linout(ind).SVD.S)
end
title('Singular Values')

% figure(2), hold on
% for i = 1:nlin
%     plot(linout(i).SVD.Xwork./max(linout(i).SVD.Xwork))
% end
% title('Work per state')

% figure(3), hold on
% for i = 1:nlin
%     plot(linout(i).SVD.Xwork./max(linout(i).SVD.Xwork))
% end
% title('Work per state')

% figure(4), hold on
% for i = 1:4
%     plot(abs(linout(i).SVD.V(i,:)))
% end
% title('Directions')

figure(5),
m = [];
for ind = 1%:4
    v = abs(linavg.SVD.V(ind,:))';
    m = [m v];
end
cats = linavg.x_desc;
% for i = 1:length(linout(1).x_desc_rm)
%     cats(i) = {cats{i}(4:6)}
% end
c = categorical(cats);
barh(c',m')
% legend('1','2','3','4')

% % First six singular values look interesting, lets think about their
% % directions
% figure, hold on
% for i = 1:nlin
%     for j = 1:6
%         Vt = [linout(i).SVD.V]';
%         vbar = Vt(i,:);
%         ubar = linout(i).SVD.U(:,i);
%         dir = vbar.*ubar';
%         plot(dir,'x-','MarkerSize',10);
%     end
% end
% legend('\sigma_1', '\sigma_2', '\sigma_3',  '\sigma_4',  '\sigma_5', '\sigma_6')
% 
% linout(1).x_desc_rm





