

fname = [glob, 'sechistory.csv'];

fpath = fullfile(comp_dir, fname);
opts=detectImportOptions(fpath);

for v=1:numel(opts.VariableNames)
    opts.VariableTypes{v} = var_type(opts.VariableNames{v}); 
end

sechistory = readtable(fpath, opts);
