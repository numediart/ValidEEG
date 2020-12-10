%% Fashionable validation framework for ERP and oscillatory source reconstruction using Fieldtrip

mfilename = 'validation_framework.m';
cd(fileparts(which(mfilename)))

%read config file
fname = 'config.json';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
config = jsondecode(str);