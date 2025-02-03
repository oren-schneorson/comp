%{
***************************************************************************
Variables types mapping
***************************************************************************
%}

if ~exist('username', 'var') || isempty(username)
    username = getenv("USER");
end
if ~exist('lib_data', 'var') || strcmp(lib_data, '')
    lib_data = fullfile('/media', username, 'D', 'data');
end

if ~exist('comp_dir', 'var') || strcmp(comp_dir, '')
    comp_dir = fullfile(lib_data, 'comp');
end

c_dir = fullfile(lib_comp, 'var_type');
d=dir(c_dir);
d = d(~[d.isdir]);
var_type = containers.Map;
type_var = containers.Map;

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
    type_var(f.name) = vars;
end
