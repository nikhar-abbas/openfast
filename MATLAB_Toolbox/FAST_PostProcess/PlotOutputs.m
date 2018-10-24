%% PlotOutputs.m
% Plots desired outputs from defined OpenFAST .out file or simulink s-function output 

%% Load fAST data

% simulink or fast?
type = 1;                   % 1 for fast outfile, 2 for simulink 
% save figures?
sfig = 0;                   % 1 to save figures
%save data?
sdata = 0;                  % 1 to save data

if type == 1
    % Load FAST output file
    file = '/Users/nabbas/Documents/TurbineModels/DTU_10MW/DTU_10MW_RWT.out';
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
            'YawBrTDxt';
};
        
%% Plot Outputs
for j = 1:length(outstr)
    fid = figure(100+j);
    plot(fastplot.Time, fastplot.(outstr{j}))
    xlabel('Time (s)')
    title(outstr{j})
    hold on
    
    % Save figures
    runtype = 'VSPI_NoFilt_step_fixed';
    if sfig == 1
        root = '/Volumes/GoogleDrive/My Drive/FOWT_ResearchCollab/Figures/GenericControl_v0/';        
        saveas(fid, [root, outstr{j}, '_', runtype], 'png')
    end
    % Save Data
    if sdata == 1
        root = '/Volumes/GoogleDrive/My Drive/FOWT_ResearchCollab/simout/GenericControl_v0/';
        save([root,runtype], 'fastout')
    end
end


