function fastout = LoadFastOut(FAST_OutFile)

% Load Model Output
% Loads FAST8 model output into a MATLAB structure to be post processed

% Open Fast Output File
% root = 'C:\Users\Nikhar\Documents\Research\TurbineModels\NREL5MW-OC3Spar\00-Base'; %File Location
% file = 'Model.out'; %File Name
% FAST_OutFile = [root filesep file];
fid = fopen(FAST_OutFile, 'r');
if fid == -1, error('Error loading file'), end

% Define Headers
n_rec = 0;                              % record keeper to keep while loop running until the header line has been found    
ind = 0;
while n_rec == 0 
    ind = ind+1;
    tline = fgetl(fid);
    w_str = strtrim(tline);
    w_str = strsplit(w_str); 
    if strcmpi(w_str{1},'time')
        n_rec = 1;
        headers = w_str;                % define headers
        
        tline = fgetl(fid);
        check = strtrim(tline);         
        if check(1) == '('              % find units line
            units = strtrim(tline);
            units = strsplit(units);    % define units
            ind = ind+1;                % index headerlines
        end
    end
end

frewind(fid)

dstr = cell(1,length(headers));                 %load data
dstr(:) = {'%f '};

dat = textscan(fid, strjoin(dstr), 'headerlines',ind);

%Put data in structure
fastout.headers = headers;
fastout.units = units;
for i = 1:length(headers)
    try
    fastout.(headers{i}) = dat{i};
    catch
        warning(['Outlist Parameter ' headers{i} ' was not loaded from the fast.out file.'])
    end
end

end