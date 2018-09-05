% FAST_LinAnalysis
% - Load and analyse a number of FAST linearization files

%% Setup Paths
Outdir = '/Users/nabbas/Documents/TurbineModels/DTU_10MW';
OutfileBase = 'DTU_10MW_RWT';
nlin = 18;



%% Load Files
for i = 1:nlin
    linfile = [Outdir filesep OutfileBase '.' num2str(i) '.lin'];
    linout(i) = ReadFASTLinear(linfile);
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

for i = 1:length(nlin)
% some setup for cleaner code (poor memory management, though)
    psi = linout(i).Azimuth;
    x_size = length(linout(i).x_rotFrame)/2;
    u_size = length(linout(i).u_rotFrame);
    y_size = length(linout(i).y_rotFrame);
    Omega = linout(i).RotSpeed;
    Omega_dot = 0;
    A = linout(i).A;
    B = linout(i).B;
    C = linout(i).C;
        C1 = C(:,1:x_size);
        C2 = C(x_size+1:end);
    D = linout(i).D;
    
% Constructing setup matrices & variables
    t1_til = [1 cos(psi) sin(psi); 1 cos(psi+2*pi/3) sin(psi+2*pi/3); 1 cos(psi+4*pi/3) sin(psi+4*pi/3)];
    T1 = [];
    for ind = 1:x_size/3
        T1 = blkdiag(T1, t1_til);
    end
    T1 = T1 * diag(linout(i).x_rotFrame(1:end/2));
    T1 = T1 + eye(x_size) + -1.*diag(linout(i).x_rotFrame(1:x_size));
    
    t2_til = [0 -sin(psi) cos(psi); 0 -sin(psi+2*pi/3) cos(psi+2*pi/3); 0 -sin(psi+4*pi/3) cos(psi+4*pi/3)];
    T2 = [];
    for ind = 1:x_size/3
        T2 = blkdiag(T2, t2_til);
    end
    T2 = T2 * diag(linout(i).x_rotFrame(1:x_size));
    
    t3_til = [0 -cos(psi) -sin(psi); 0 -cos(psi+2*pi/3) -sin(psi+2*pi/3); 0 -cos(psi+4*pi/3) -sin(psi+4*pi/3)];
    T3 = [];
    for ind = 1:x_size/3
        T3 = blkdiag(T3, t2_til);
    end
    T3 = T3 * diag(linout(i).x_rotFrame(1:x_size));
    
    T1c = [];
    for ind = 1:u_size/3
        T1c = blkdiag(T1c, t1_til);
    end
    T1c = T1c * diag(linout(i).u_rotFrame);
    T1c = T1c + eye(u_size) + -1.*diag(linout(i).u_rotFrame);
    
    T1o = [];
    for ind = 1:y_size/3
        T1o = blkdiag(T1o, t1_til);
    end
    T1o = T1o * diag(linout(i).y_rotFrame);
    T1o = T1o + eye(y_size) + -1.*diag(linout(i).y_rotFrame);

    
% Nonrotating state matrices
    linout(i).Anr = blkdiag(inv(T1), inv(T1)) * (A * [T1 zeros(size(T1)); Omega*T2 T1] - [Omega*T2 zeros(size(T2)); Omega^2*T3 + Omega_dot*T2 2*Omega*T2]);
    linout(i).Bnr = blkdiag(inv(T1), inv(T1)) * B * T1c;
    linout(i).Cnr = inv(T1o) * [C1*T1 + Omega*C2*T2 C2*T1];
    linout(i).Dnr = inv(T1o)*D*T1c;
    
end
