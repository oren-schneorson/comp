%{
script to replace compute_merton_vars.m

cleaning of data should be done in advance
all values should be in USD
all fundq vars should be in millions (already are...?)
mcap should be in millions

# fundq: holes sholud be filled in advance, or perhaps
I don't need to fill holes, because I'm using retime anyway

# fundq: edges cropped
# leave only unique observations
# save fundq as fundq_proc


secd_proc contains only the relevant variables

firms with missing data should be cleaned out in advance
and tracked



%}
clear; clc

% set root
username = getenv("USER");
root = mfilename('fullpath');
if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    % change this line if running on a different system
    root = fullfile('/home', username, 'Documents/MATLAB/comp');
    cd(root)
    
else
    root = fileparts(root);
end

cd(root)


comp = 'comp';
fundq = 'g_fundq';

matlab_dir = fullfile('/home', username, 'Documents', 'MATLAB');
lib_data = fullfile('/media', username, 'D', 'data');
temp_lib_data = fullfile('/media', username, 'D', 'data_temp4backup');


comp_dir = fullfile(lib_data, comp);
% comp_dir = fullfile(temp_lib_data, comp);
fundq_dir = fullfile(comp_dir, fundq);
fig_dir = fullfile(root, 'gfx', 'fundq');

addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))




fx_daily_fpath = fullfile(lib_data, 'fx', 'fx_daily.csv');
fx_daily = readtable(fx_daily_fpath);

%%
clc

plot_figure = false;

redo_all = true;

query_mode = false;
% query = '015532'; % na
query = '015532'; % global

fpath = fullfile(comp_dir, 'liab_structure.xlsx');
liab_structure = readtable(fpath);
liab_structure = liab_structure(liab_structure.usgaap==1,:);
components_lctq = {'apq', 'lcoq', 'dlcq', 'txpq'};
components_ltq = {'lctq', 'loq', 'dlttq', 'txditcq'};
%
fundq_vars = {'apq', 'lcoq', 'dlcq', 'txpq', 'dlttq', 'loq', 'txditcq', 'xaccq', 'dd1q', 'npq'};
%fundq_vars = {'apq', 'lcoq', 'dlcq', 'txpq', 'dlttq', 'loq', 'txditcq'};


d = dir(fundq_dir);
d = d(~[d.isdir]);
d = d(~cellfun(@(c) contains(c, '#'), {d.name}));

fpath = 'gvks_neg_liab.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvks_neg_liab = readtable(fpath, opts);
gvks_neg_liab = gvks_neg_liab{:,1};


for i_d = 1:numel(d)
    
    clc
    i_d/numel(d)*100

    fname = d(i_d).name;
    gvkey = strrep(fname, '.csv', '');
    
    if ~ismember(gvkey, gvks_neg_liab)
        %continue
    end
    
    if ~strcmp(gvkey, query) && query_mode
        continue
    end
    
    fpath_save = fullfile(comp_dir, [fundq, '_proc'], fname);
    if exist(fpath_save, 'file') == 2 && ~redo_all && ~query_mode
        continue
    end
    
    T = load_fundq(gvkey, fx_daily, fundq_dir);
    if isempty(T)
        continue
    end
    writetable(T, fpath_save)
    
    if plot_figure
        T = T(T.datadate.Year < 2021,:);
        clf
        subplot(1,2,1)
        hold on
        area(T.datadate, T{:, components_ltq}, 'LineStyle', 'none')
        plot(T.datadate, T.ltq, 'r--', 'LineWidth', 1.5)
        legend([components_ltq, {'ltq'}], 'Location', 'northwest')
        title('ltq and its components')
        
        if sum(T.ltq<0) >0
            vline(T.datadate(T.ltq<0))
        end
        
        subplot(1,2,2)
        hold on
        area(T.datadate, T{:, components_lctq}, 'LineStyle', 'none')
        plot(T.datadate, T.lctq, 'r--', 'LineWidth', 1.5)
        legend([components_lctq, {'lctq'}], 'Location', 'northwest')
        title('lctq and its components')
        if sum(T.lctq<0)>0
            vline(T.datadate(T.lctq<0))
        end
        
        fpath_fig = fullfile(fig_dir, [gvkey, '.png']);
        export_fig(fpath_fig)
        
        
    end
    
    
    
    

end
return





