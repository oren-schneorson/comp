function [opts] = gen_opts(fpath_table, date_vars, InputFormat, var_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%{
***************************************************************************
Variables types mapping
***************************************************************************
%}

if isempty(var_type)
comp_dir = '/mnt/D/data/comp';
c_dir = fullfile(comp_dir, 'var_type');
d=dir(c_dir);
d([d.isdir]) = [];
var_type = containers.Map;

for f=d'
    fpath=fullfile(c_dir, f.name);
    opts=detectImportOptions(fpath);
    opts.VariableNamesLine = 1;
    opts.DataLines = [2, inf];
    opts.VariableNames= {'var'};
    opts.Delimiter= ',';
    vars = readtable(fpath,  opts);
    for v=vars.var'
        var_type(v{1})=f.name;
    end
end

end


opts=detectImportOptions(fpath_table);
for v=1:numel(opts.VariableNames)
    %opts.VariableNames{v}
    opts.VariableTypes{v} = var_type(opts.VariableNames{v}); 
end

if isempty(InputFormat)
    InputFormat = 'yyyy-MM-dd';
end
if isempty(date_vars)
    InputFormat = {'datadate'};
end

%opts=setvaropts(opts, date_vars, 'InputFormat', InputFormat);



end

