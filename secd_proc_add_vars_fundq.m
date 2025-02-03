
clear; clc

%{
***************************************************************************
Preliminaries
***************************************************************************
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


plot_figure = true;
redo_all = true;
missing_data = false;
% when both false, continues redo_all from where stopped.

query_mode = false;
query = '015325';

min_year = 1900;

min_num_obs = 10;

load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

clear y Fs

comp = 'comp';
glob = '';
PycharmProjects_dir = fullfile('/home', username, 'PycharmProjects');
matlab_dir = fullfile('/home', username, 'Documents', 'MATLAB');
lib_data = fullfile('/media', username, 'D', 'data');

comp_dir = fullfile(lib_data, comp);
sechistory_dir = fullfile(comp_dir, [glob, 'sechistory']);
secd_dir = fullfile(comp_dir, [glob, 'secd']);
secm_dir = fullfile(comp_dir, [glob, 'secm']);
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
secm_proc_dir = fullfile(comp_dir, [glob, 'secm_proc']);
fundq_dir = fullfile(comp_dir, [glob, 'fundq']);
fig_dir = fullfile(comp_dir, 'figs', [glob, 'secd_proc']);


addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))



run('var_types.m') % Variables types mapping


run('load_company.m')
run('load_security.m')
run('load_sechistory.m')
run('load_currency.m')

%load(fullfile(lib_data, 'ir', 'TB.mat'))
%run('TBill.m') % risk free rate table, TB
%run('identification.m') % load company/security specific variables (dicts)
fpath = fullfile(lib_data, 'ir', 'DTB');

[R_f, R_f_metadata] = load_data_files('DTB3', fpath, 'FRED');

fpath = fullfile('/home', username, 'PycharmProjects', 'wrds_get_data', 'xpressfeed', 'metadata.csv');
opts = detectImportOptions(fpath);
opts.VariableTypes = repmat({'char'}, size(opts.VariableTypes));

compustat_metadata = readtable(fpath, opts);




units = containers.Map;
units('thousand') = 1e3;
units('million') = 1e6;
units('billion') = 1e9;
units_keys = units.keys();

for key=units.keys()
    key = char(key);
    units(lower(key)) = units(key);
    units(upper(key)) = units(key);
    units([upper(key(1)), key(2:end)]) = units(key);

    units([key, 's']) = units(key);
    units([upper(key), 'S']) = units(key);
    units([upper(key(1)), key(2:end), 's']) = units(key);
end



%{
This loop processes comapny level data and saves it to a separate library
from the raw data. It computes market value of each firm by adding the
market value of each security that participates in earning report (epf =
'Y'). ret is defined as the change in total market value of the firm:
ret = MCAP1/MCAP0

The code filters for each legilation date only the firms
that have at least 30 unique observations in terms of ret.

%}

epf = containers.Map; % for tracking extrapolation of cshoc
epf('g_') = 'epf';
epf('') = 'mkvalincl';
%%
clc

vars_fundq_proc2add = {'datadate', 'rdq' 'atq', 'ltq'};


gvkey_i_start = 1;
%{
fpath = fullfile(PycharmProjects_dir, 'comp', 'gvkey_set_lctq.csv')
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvkeys = readtable(fpath, opts);
gvkeys = gvkeys{:,1};
%}


gvkeys = dir(secd_proc_dir);
gvkeys = gvkeys(~[gvkeys.isdir]);
gvkeys = gvkeys(~contains({gvkeys.name}, '#'));

gvkeys=cellfun(@(c) c(1:end-4), {gvkeys.name},...
    'UniformOutput', false);
gvkeys = gvkeys(~cellfun(@isempty, gvkeys));
%}

% fpath = fullfile(PycharmProjects_dir, 'comp', 'redo_gvkeys_secd_proc.csv')
fpath = fullfile(PycharmProjects_dir, 'comp', 'gvkeys_no_fundq.csv')
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
redo_gvkeys = readtable(fpath, opts);
redo_gvkeys = redo_gvkeys{:,1};
gvkeys = redo_gvkeys;

%return
%}

loop_times = NaN(size(gvkeys));
loop_times(1) = 0;
ngvkeys = numel(gvkeys);

%for gvk_i=flip(1:ngvkeys)
%(ngvkeys-1)/2+1:ngvkeys
% 1:(ngvkeys-1)/2
%for gvkey_i=ngvkeys/2+1:ngvkeys


for gvkey_i=gvkey_i_start:ngvkeys
    
    gvkey_i_start=gvkey_i;
    tic % measure time
    gvkey = gvkeys{gvkey_i}; % get gvkey
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');

    if ~strcmp(gvkey, query) && query_mode
        continue
    end
    
    % get data on all securities of company by gvkey
    fpath_sechistory = fullfile(sechistory_dir, sprintf('%s.csv', gvkey));
    company_i = company(ismember(company.gvkey, gvkey), :);
    security_i = security(ismember(security.gvkey, gvkey), :);
    sechistory_i = sechistory(ismember(sechistory.gvkey, gvkey), :);    
    
    Country = unique(company_i.loc);
    
    
    %{
    % 007228, 003221 is an excpetion. Better data in na file...
    if ismember(Country, {'USA', 'CAN'}) && ~strcmp(lib, 'lib_na')
        continue
    elseif ~ismember(Country, {'USA', 'CAN'}) && strcmp(lib, 'lib_na')
        continue
    end
    %}
    
    fpath_secd = fullfile(secd_dir, sprintf('%s.csv', gvkey));
    fpath_secm = fullfile(secm_dir, sprintf('%s.csv', gvkey));
    fpath_fundq = fullfile(fundq_dir, sprintf('%s.csv', gvkey));
    fpath_secd_proc = fullfile(secd_proc_dir, sprintf('%s.csv', gvkey));
    fpath_secm_proc = fullfile(secm_proc_dir, sprintf('%s.csv', gvkey));
    fpath_fig = fullfile(fig_dir, [gvkey, '.png']);
    
    if exist(fpath_secd_proc, 'file') ~=2
        fprintf('%s not present in secd_proc.\n', gvkey)
        continue
    end

    
    if exist(fpath_secd_proc, 'file')==2 && ~redo_all && ~query_mode    
        fprintf('%s already present in secd_proc.\n', gvkey)
        continue
    end
    
    finfo = dir(fpath_secd_proc);
    if finfo.bytes <= 292
        continue
    end
      
    % table secd is read from csv with options
    DateFormat = 'yyyy-MM-dd';   
    
    % fundq information:
    opts = gen_opts(fpath_fundq, {'datadate', 'rdq'}, DateFormat, var_type);
    fundq=readtable(fpath_fundq, opts);

    if isempty(fundq)
        continue
    end
    
    opts = gen_opts(fpath_secd, {'datadate'}, DateFormat, var_type);
    secd=readtable(fpath_secd, opts);
    opts = gen_opts(fpath_secd_proc, {'datadate'}, DateFormat, var_type);
    secd_proc=readtable(fpath_secd_proc, opts);
    
    if numel(glob) == 0 % north america
        secd.qunit = ones(size(secd, 1), 1);
    end
    
    if any(ismember(secd_proc.Properties.VariableNames, {'atq', 'ltq'}))
        continue
    end
    

    if ~isempty(fundq)
    
        vars_fundq = {'datadate', 'rdq' 'atq', 'ltq', 'cshoq'};
        fundq = fundq(:, vars_fundq);
        
        % TODO: change this to work automatically with metadata
        fundq.atq = fundq.atq * 1e6;
        fundq.ltq = fundq.ltq * 1e6;
        fundq.ltq = fundq.cshoq * 1e6;
        %fundq.book_value = fundq.atq-fundq.ltq;
        %fundq = removevars(fundq, {'atq', 'ltq'});

        % if date not reported, assume one month lag
        idx = isnat(fundq.rdq);
        fundq.rdq(idx) = eomdate(fundq.datadate(idx) + calmonths(1));
        
        % keep first report for any datadate
        fundq = sortrows(fundq, {'datadate', 'rdq'}, 'ascend');
        [~, idx] = unique(fundq.datadate);
        fundq = fundq(idx, :);
        fundq.datadate = fundq.rdq;
        fundq = removevars(fundq, 'rdq');
        fundq = sortrows(fundq, 'datadate');

        % retime fundq (robust merge, Keys: datadate)
        fundq = table2timetable(fundq);
        fundq = retime(fundq, 'daily', 'previous');
        fundq = timetable2table(fundq);

        secd = outerjoin(secd, fundq,...
            'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');        
        secd_proc = outerjoin(secd_proc, fundq(:, vars_fundq_proc2add([1, 3:4])),...
            'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');        
              
        idx = isnan(secd.cshoc);
        if ~all(idx)
            % complete cshoc from cshocq only at tails
            idx_ = idx & secd.datadate < min(secd.datadate(~idx));
            idx_ = idx_ & secd.datadate > max(secd.datadate(~idx));
        else
            idx_ = idx;
        end

        if sum(idx) ~= sum(idx_)        
            secd.cshoc(idx_) = secd.cshoq(idx_);
            clear idx idx_
        end
        secd = removevars(secd, 'cshoq');
        
    end


    secd = secd(secd.datadate.Year>=min_year, :);
    secd = secd(~isnat(secd.datadate), :);
    
    if size(secd, 1) < 2 || all(isnan(secd.prccd))
        continue
    end
    
    if all(isnan(secd.trfd))
        secd.trfd = ones(size(secd,1), 1);
    end
    
    secd.iid = cellfun(@(c) pad(c, 2, 'left', '0'), secd.iid, 'UniformOutput', false);
    
    if ismember('gvkey', secd.Properties.VariableNames)
        secd = removevars(secd, 'gvkey');
    end

    secd = secd(~isnan(secd.prccd),:);
    secd.prccd = secd.prccd./secd.qunit;
    secd = sortrows(secd, {'iid', 'datadate'});
    iids = unique(secd.iid); % list of unique iids
    
    if numel(iids) == 1
        %idx_secd_proc = secd_proc.retd == 1;
        secd.intraday_retd = secd.prccd./secd.prcod;
        idx_secd = ~isnan(secd.intraday_retd) & ~isinf(secd.intraday_retd);
        secd_proc = outerjoin(secd_proc, secd(idx_secd, {'datadate', 'intraday_retd'}),...
            'Keys', {'datadate'}, 'MergeKeys', true, 'Type', 'left');
        
        idx_secd_proc = ...
            ~isnan(secd_proc.intraday_retd) &...
            ~isinf(secd_proc.intraday_retd);
        idx_secd_proc = idx_secd_proc & (secd_proc.retd == 1 | secd_proc.retd == 0 | isnan(secd_proc.retd) | isinf(secd_proc.retd));
        
        secd_proc.retd(idx_secd_proc) = secd_proc.intraday_retd(idx_secd_proc);
        secd_proc = removevars(secd_proc, 'intraday_retd');
        secd_proc.retd(idx_secd_proc)

    end
    
    writetable(secd_proc, fpath_secd_proc)

   
    [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', gvkey]);
    clc
    time_left_
    if query_mode
        return
    end




end


sound(end_sound, end_sound_Fs);

