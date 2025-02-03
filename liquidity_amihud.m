clear; clc


%{
The first step in implementing the Merton DD model is to estimate sigma_E from
either historical stock returns data or from option-implied volatility data.
%}

addpath('/home/oren/Documents/MATLAB/altmany-export_fig-410f0ad')
addpath('/home/oren/Documents/MATLAB/my_functions')

% set root
root = mfilename('fullpath');
if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    % change this line if running on a different system
    root = '/home/oren/Documents/MATLAB/comp';
    cd(root)
    
else
    root = fileparts(root);
end

cd(root)


glob = '';
comp_dir = '/media/oren/D/data/comp';
liquidity_dir = fullfile(comp_dir, [glob, 'liquidity'], 'amihud');
secd_dir = fullfile(comp_dir, [glob, 'secd']);
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
sechistory_dir = fullfile(comp_dir, [glob, 'sechistory']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq']);
fig_dir = fullfile(comp_dir, 'figs', [glob, 'liquidity']);


fpath = fullfile(comp_dir, [glob, 'company.csv']);
company = readtable(fpath);
fpath = fullfile(comp_dir, [glob, 'security.csv']);
security = readtable(fpath);

run('var_types.m')


load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

DateFormat = 'yyyy-MM-dd';
min_year = 1900;

%%

clc

%locs = {'USA', 'CAN'};
locs = {'USA', };

redo_all = false;
plot_figure = true;

query_mode = false;
%query = '015634';
query = '025056';

%
fpath = '/home/oren/PycharmProjects/comp/gvkey_set_lctq.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvks = readtable(fpath, opts);
gvks = gvks{:,1};
%}


gvkeys = dir(secd_proc_dir);
gvkeys = gvkeys(~[gvkeys.isdir]);
gvkeys = gvkeys(~contains({gvkeys.name}, '#'));
%}

ngvkeys = numel(gvkeys);
loop_times = NaN(ngvkeys, 1);
loop_times(1) = 0;

liquidity_tab = table;


for gvkey_i= 1:ngvkeys
    clc
    tic % measure time
    
    fname = gvkeys(gvkey_i).name;
    gvkey = strrep(fname, '.csv', ''); % get gvkey
    if ~ismember(gvkey, gvks)
        %continue
    end
    
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');
    gvkey = str2double(gvkey);
    
    % get data on all securities of company by gvkey
    fpath_sechistory = fullfile(sechistory_dir, fname);
    company_i = company(ismember(company.gvkey, gvkey), :);
    security_i = security(ismember(security.gvkey, gvkey), :);
    sechistory_i = readtable(fpath_sechistory);
    
    loc = company_i{:, 'loc'};
    
    if ~strcmp(gvkey, query) && query_mode
        continue
    end
    
    if ~ismember(loc, locs) && numel(glob) == 0
        continue
    end

    fpath_secd = fullfile(secd_dir, fname);    
    fpath_secd_proc = fullfile(secd_proc_dir, fname);
    fpath_liquidity = fullfile(liquidity_dir, fname);
        
    if exist(fpath_liquidity, 'file') == 2 && ~redo_all && ~query_mode
    	fprintf('%s already done\n', gvkey)
        continue
    end
    
    opts = gen_opts(fpath_secd_proc, {'datadate'}, DateFormat, var_type);
    secd_proc = readtable(fpath_secd_proc, opts);
    
    
    
    if numel(unique(secd_proc.datadate)) ~= size(secd_proc,1)
        error('datadate not unique, %d: %s', gvkey_i, gvkey)
    end
    
    
    if numel(glob) == 0 % north america
        secd_proc.qunit = ones(size(secd_proc, 1), 1);
    end

    secd_proc = secd_proc(secd_proc.datadate.Year>=min_year, :);
    secd_proc = secd_proc(~isnat(secd_proc.datadate), :);
    
    if size(secd_proc, 1) < 10 || all(isnan(secd_proc.retd))
        continue
    end
    
    secd_proc.retd = (secd_proc.retd-1)*100;
    amihud = abs(secd_proc.retd)./secd_proc.volume;
    
    amihud = mean(amihud(~isinf(amihud)), 'omitnan');
    gamma = -cov(secd_proc.retd(2:end), secd_proc.retd(1:end-1));
    gamma = gamma(2);
    mcap_mean = mean(secd_proc.mcap, 'omitnan');
    retd_mean = mean(secd_proc.retd, 'omitnan');
    retd_std = std(secd_proc.retd, 'omitnan');
    volume_mean = mean(secd_proc.volume, 'omitnan');

    
    vars_company = {'loc','gsector','gind','gsubind'};
    
    liquidity_tab_i = table(gvkey, amihud, gamma, mcap_mean, volume_mean, retd_mean, retd_std,...
        'VariableNames', {'gvkey', 'amihud', 'gamma', 'mcap_mean', 'volume_mean', 'retd_mean', 'retd_std'});
    liquidity_tab_i = [liquidity_tab_i, company_i(:, vars_company)];
    liquidity_tab = [liquidity_tab; liquidity_tab_i];
    
    loop_times(gvkey_i) = toc;
    clc
    [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', strrep(fname, '.csv', '')]);
    time_left_


    if query_mode
        return
    end

    continue
    
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
   

    loop_times(gvkey_i) = toc;
        
    if query_mode
        return
    end

end






sound(end_sound, end_sound_Fs);

return

%%

clc
clf
hold on
%{
% gamma
idx = true(size(liquidity_tab.gamma));
idx = idx & liquidity_tab.mcap_mean > 1e9;
idx = idx & liquidity_tab.gamma < prctile(liquidity_tab.gamma, 98);
idx = idx & liquidity_tab.gamma > prctile(liquidity_tab.gamma, 2);
scatter(liquidity_tab.mcap_mean(idx), liquidity_tab.gamma(idx))
%}

% amihud
idx = true(size(liquidity_tab.amihud));
idx = idx & liquidity_tab.mcap_mean > 1e9;
idx = idx & liquidity_tab.amihud < prctile(liquidity_tab.amihud, 50);
idx = idx & liquidity_tab.amihud > prctile(liquidity_tab.amihud, 2);
scatter(liquidity_tab.mcap_mean(idx)/1e6, liquidity_tab.amihud(idx))

ax = gca;
ax.XLim(2)
x = liquidity_tab.mcap_mean(idx)/1e6;
    func = @(x, a, b) a.*exp(-b*x);
a = max(liquidity_tab.amihud(idx));
b = .001;
scatter(x, func(x, a, b));

J = sum((func(x, a, b)-liquidity_tab.amihud(idx)).^2)


% negative gamma is ok, it means liquid...
%idx = idx & liquidity_tab.gamma < 0;
%liquidity_tab(idx,:)



















