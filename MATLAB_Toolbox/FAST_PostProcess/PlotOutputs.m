%% PlotOutputs.m
% Plots desired outputs from defined OpenFAST .out file or simulink s-function output 

%% Load fAST data

% simulink or fast?
type = 2;                   % 1 for fast, 2 for simulink

if type == 1
    % Load FAST output file
    file = ['/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU_10MW_RWT.out'];
    fastout = LoadFastOut(file);
elseif type == 2
    fastout = simout;
else
    error('Make sure your simulation type is defined properly') 
end

fastplot = fastout;


%% Define desired FAST outputs to plot
outstr = {'Wind1VelX';
    'BldPitch1';
             'GenTq';
            'RotSpeed';
            'GenPwr';
};
        
%% Plot Outputs
for j = 1:length(outstr)
    figure(100+j)
    plot(fastplot.Time, fastplot.(outstr{j}))
    xlabel('Time (s)')
    title(outstr{j})
    hold on
end


