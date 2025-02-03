clear; clc


%{
The first step in implementing the Merton DD model is to estimate sigma_E from
either historical stock returns data or from option-implied volatility data.
%}

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


addpath(fullfile('/home', username, 'Documents', 'MATLAB', 'altmany-export_fig-410f0ad'))
addpath(fullfile('/home', username, 'Documents', 'MATLAB', 'my_functions'))


glob = '';

lib_data = fullfile('/media', username, 'D', 'data');
comp_dir = fullfile(lib_data, 'comp');
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);


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


fpath = fullfile(comp_dir, [glob, 'company.csv']);
company = readtable(fpath);

    
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq']);

fig_dir = fullfile(comp_dir, 'figs', [glob, 'liquidity']);

%{
fpath = fullfile('/home', username, 'PycharmProjects', 'comp', 'gvkey_set_lctq.csv');

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
lo

for i_file= 1:n
    clc
    tic % measure time
    pcent_ = i_file/n*100;
    time_left = (n-i_file-1)*mean(loop_times, 'omitnan')/60/60/24;
    time_left = datestr(time_left, 'HH:MM:SS.FFF');
    
    fname = d.name(i_file);
    gvkey = strrep(fname, '.csv', ''); % get gvkey
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');
    
    if ~strcmp(gvkey, query) && query_mode
        continue
    end
    
    idx_company = company.gvkey==gvkey;
    %{
    loc = company{idx_company, 'loc'};
    if ~ismember(loc, {'USA', 'CAN'}) && numel(glob) == 0
        %continue
    end
    %}

    
        
    fpath_secd_proc = fullfile(secd_proc_dir, fname);    
        
    if exist(fpath_merton_liquidity, 'file') && ~redo_all && ~query_mode
    	fprintf('%s already done\n', gvkey)
        continue
    end
    
    liq_mean = table;
    for Ticker = Tickers'
        fname = [char(Ticker), '.csv'];
       
        fpath_CDS = fullfile(by_ticker_DocClause, fname);
        
        if exist(fpath_CDS, 'file') ~= 2
            continue
        end
        
        fpath_median = fullfile(...
            by_ticker_DocClause_liquidity, 'median');
        if ~exist(fpath_median, 'dir')
            mkdir(fpath_median)
        end
        
        fpath_mean = fullfile(...
            by_ticker_DocClause_liquidity, 'mean');
        if ~exist(fpath_mean, 'dir')
            mkdir(fpath_mean)
        end
        fpath_mean = fullfile(fpath_mean, fname);
        fpath_median = fullfile(fpath_median, fname);
        
        if ~exist(fpath_median, 'file') == 2
            
        T = readtable(fpath_CDS);
        [~, idx] = unique(T(:, {'Date', 'Tenor'}), 'rows');

        liq_mean_i = varfun(@(date_)...
            mean(days252bus(...
            date_(1:end-1),...
            date_(2:end)),...
            'omitnan'),...
            T(idx, {'Tenor', 'Date'}),...
            'GroupingVariables', 'Tenor');
        liq_median_i = varfun(@(date_)...
            median(days252bus(...
            date_(1:end-1),...
            date_(2:end)),...
            'omitnan'),...
            T(idx, {'Tenor', 'Date'}),...
            'GroupingVariables', 'Tenor');
        
        writetable(liq_median_i, fpath_median);
        writetable(liq_mean_i, fpath);
        
        else
            liq_median_i = readtable(fpath_median);
            liq_mean_i = readtable(fpath_mean);
        
        end
        
        if isempty(liq_median_i)
            continue
        end

    if plot_figure

        liq_median_i.Tenor_y = cellfun(@(c) steps_dict(c),...
            liq_median_i.Tenor);
        liq_mean_i.Tenor_y = cellfun(@(c) steps_dict(c),...
            liq_mean_i.Tenor);
        liq_mean_i = sortrows(liq_mean_i, 'Tenor_y');
        liq_median_i = sortrows(liq_median_i, 'Tenor_y');
        clf
        plot(liq_mean_i.Tenor_y, liq_mean_i.Fun_Date)
        xlabel('Tenor')
        ylabel('Avg. bus. days between trades')
        title(Ticker)
        fname = [char(Ticker), '.png'];
        fpath_fig = fullfile(fig_dir_CDS, fname);
        export_fig(fpath_fig)

    end
    end
   

    loop_times(i_file) = toc;
        
    if query_mode
        return
    end

end



fpath = fullfile(comp_dir, [glob, 'errors'], 'merton_error.csv');
writetable(liq_mean, fpath);



sound(end_sound, end_sound_Fs);

return





















