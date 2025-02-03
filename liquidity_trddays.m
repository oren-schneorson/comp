clear; clc


%{
The first step in implementing the Merton DD model is to estimate sigma_E from
either historical stock returns data or from option-implied volatility data.
%}

addpath('/home/oren/Documents/MATLAB/altmany-export_fig-410f0ad')
addpath('/home/oren/Documents/MATLAB/my_functions')

glob = '';
comp_dir = '/media/oren/D/data/comp';
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
root = '/home/oren/Documents/MATLAB/comp';
cd(root)


load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

%%

clc

redo_all = false;
plot_figure = true;

query_mode = false;
%query = '015634';
query = '025056';


fpath_company = fullfile(comp_dir, [glob, 'company.csv']);
company = readtable(fpath_company);

    
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq']);
liquidity_dir = fullfile(comp_dir, [glob, 'liquidity']);

fig_dir = fullfile(comp_dir, 'figs', [glob, 'liquidity']);

%{
fpath = '/home/oren/PycharmProjects/comp/gvkey_set_lctq.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvks = readtable(fpath, opts);
gvks = gvks{:,1};
%}

d = dir(secd_proc_dir);
d = d(~[d.isdir]);
d = d(~contains({d.name}, '#'));

n = numel(d);
loop_times = NaN(n, 1);
loop_times(1) = 0;

%%
for i_file= 1523:n
    clc
    tic % measure time
    pcent_ = i_file/n*100;
    time_left = (n-i_file-1)*mean(loop_times, 'omitnan')/60/60/24;
    time_left = datestr(time_left, 'HH:MM:SS.FFF');
    
    fprintf('%.2f%%, time left: %s\n', i_file/n*100, time_left)
    
    fname = d(i_file).name;
    gvkey = strrep(fname, '.csv', ''); % get gvkey
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');
    
    fpath_secd_proc = fullfile(secd_proc_dir, fname);    
    fpath_liquidity = fullfile(liquidity_dir, fname);    
    
    if ~strcmp(gvkey, query) && query_mode
        continue
    elseif exist(fpath_secd_proc, 'file') ~= 2
        continue
    elseif exist(fpath_liquidity, 'file') == 2 && ~redo_all
        continue
    end
    
    idx_company = company.gvkey==str2double(gvkey);
    %{
    loc = company{idx_company, 'loc'};
    if ~ismember(loc, {'USA', 'CAN'}) && numel(glob) == 0
        %continue
    end
    %}

          
    fpath_median = fullfile(...
        liquidity_dir, 'median_trddays');
    fpath_mean = fullfile(...
        liquidity_dir, 'mean_trddays');

    if ~exist(fpath_median, 'dir')
        mkdir(fpath_median)
    end
    if ~exist(fpath_mean, 'dir')
        mkdir(fpath_mean)
    end
    
    fpath_mean = fullfile(fpath_mean, fname);
    fpath_median = fullfile(fpath_median, fname);

    if exist(fpath_median, 'file') ~= 2 || true

    T = readtable(fpath_secd_proc);

    liq_mean_i = varfun(@(date_)...
        mean(days252bus(...
        date_(1:end-1),...
        date_(2:end)),...
        'omitnan'),...
        T(:, {'datadate'}));
    liq_median_i = varfun(@(date_)...
        median(days252bus(...
        date_(1:end-1),...
        date_(2:end)),...
        'omitnan'),...
        T(:, {'datadate'}));

    writetable(liq_median_i, fpath_median);
    writetable(liq_mean_i, fpath_mean);

    else
        liq_median_i = readtable(fpath_median);
        liq_mean_i = readtable(fpath_mean);

    end

    if isempty(liq_median_i)
        continue
    end

    
    loop_times(i_file) = toc;
        
    if query_mode
        return
    end

end


sound(end_sound, end_sound_Fs);

return





















