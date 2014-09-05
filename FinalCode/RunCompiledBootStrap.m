function RunCompiledBootStrap(command_line_arg)
% /Applications/MATLAB_R2012b.app/bin/mcc -R -nodisplay -R -nojvm -R -singleCompThread -m -v -w enable -d ../CompiledCode RunCompiledBootStrap.m 

message = ['Command line string: ',command_line_arg];
disp(message)
% Parse the comand line for the argument
params_file = sscanf(command_line_arg,'%s');
fprintf('Parameters file name: %s\n',params_file);
% read the file
fid = fopen(params_file,'r');
line = fgetl(fid)
% split the line up
equalIndex = strfind(line,'=')
DataPath = strtrim(sscanf(line(equalIndex+1:end),'%s'))
fprintf(1,'Using data: %s\n',DataPath);

% load the data
%load(DataPath)
Results = CycleOverVoxelsProcessBootstrap(DataPath);