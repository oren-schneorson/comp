
clear; clc

%{
***************************************************************************
Preliminaries
***************************************************************************
%}

plot_figure = true;
redo_all = false;
missing_data = false;
% when both false, continues redo_all from where stopped.

query_mode = true;
query = '023133';

min_year = 1900;
DateFormat = 'yyyy-MM-dd';            



min_num_obs = 10;

load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

clear y Fs

comp = 'comp';
glob = '';
matlab_dir = '/home/oren/Documents/MATLAB';
data_dir = '/media/oren/D/data';
root = fullfile(matlab_dir, comp);
comp_dir = fullfile(data_dir, comp);
sechistory_dir = fullfile(comp_dir, [glob, 'sechistory']);

secd_dir = fullfile(comp_dir, [glob, 'secd']);
secm_dir = fullfile(comp_dir, [glob, 'secm']);

secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
secm_proc_dir = fullfile(comp_dir, [glob, 'secm_proc']);

fundq_dir = fullfile(comp_dir, [glob, 'fundq']);

fig_dir = fullfile(comp_dir, 'figs', [glob, 'secm_proc']);


addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))


cd (root)
addpath('/home/oren/Documents/MATLAB/my_functions')

run('var_types.m') % Variables types mapping


run('load_company.m')
run('load_security.m')
run('load_sechistory.m')
run('load_currency.m')

load('/media/oren/D/data/ir/TB.mat')
%run('TBill.m') % risk free rate table, TB
%run('identification.m') % load company/security specific variables (dicts)

fpath = '/home/oren/PycharmProjects/wrds_get_data/xpressfeed/metadata.csv';
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

epf = containers.Map; % for tracking extrapolation of cshom
epf('g_') = 'epf';
epf('') = 'mkvalincl';
%%
clc


fpath = '/home/oren/PycharmProjects/comp/gvkey_set_lctq.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvkeys = readtable(fpath, opts);
gvkeys = gvkeys{:,1};


%
gvkeys = dir(secm_dir);
gvkeys = gvkeys(~[gvkeys.isdir]);
gvkeys = gvkeys(~contains({gvkeys.name}, '#'));

gvkeys=cellfun(@(c) c(1:end-4), {gvkeys.name},...
    'UniformOutput', false);
gvkeys = gvkeys(~cellfun(@isempty, gvkeys));
%}

%{
fpath = '/home/oren/PycharmProjects/comp/redo_gvkeys_secm_proc.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
redo_gvkeys = readtable(fpath, opts);
redo_gvkeys = redo_gvkeys{:,1};
gvkeys = redo_gvkeys;
%}


loop_times = NaN(size(gvkeys));
loop_times(1) = 0;
ngvkeys = numel(gvkeys);
gvkey_i_start = (ngvkeys-1)/2+1;
gvkey_i_start = 1;

%for gvkey_i=gvkey_i_start:(ngvkeys-1)/2
for gvkey_i=flip(gvkey_i_start:ngvkeys)
%for gvkey_i=gvkey_i_start:ngvkeys
    
    tic % measure time
    gvkey_i_start=gvkey_i;
    gvkey = gvkeys{gvkey_i}; % get gvkey
    gvkey = pad(num2str(gvkey), 6, 'Left', '0');

    if ~strcmp(gvkey, query) && query_mode
        continue
    end
      
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
    fpath_secm_proc = fullfile(secm_proc_dir, sprintf('%s.csv', gvkey));
    fpath_fig = fullfile(fig_dir, [gvkey, '.png']);
    
    if exist(fpath_secm, 'file') ~=2
        fprintf('%s not present in secm.', gvkey)
        continue
    end
    
    if exist(fpath_secm_proc, 'file')==2 && ~redo_all && ~query_mode    
        loop_times(gvkey_i) = toc;
        %
        %current_time = now;
        %current_time = datetime(current_time,'ConvertFrom','datenum');
        clc
        [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', gvkey]);
        time_left_        
        continue
    end
    
    finfo = dir(fpath_secm);
    if finfo.bytes <= 292
        continue
    end
    
    if numel(glob) == 0 % north america
        finfo = dir(fpath_secm);
        if finfo.bytes <= 292
            continue
        end
    end
    
    % get data on all securities of company gvk
    fpath_sechistory = fullfile(sechistory_dir, sprintf('%s.csv', gvkey));
    company_i = company(ismember(company.gvkey, gvkey), :);
    security_i = security(ismember(security.gvkey, gvkey), :);
    sechistory_i = sechistory(ismember(sechistory.gvkey, gvkey), :);       
    Country = unique(company_i.loc);

    % fundq information:
    opts = gen_opts(fpath_fundq, {'datadate', 'rdq'}, DateFormat, var_type);
    fundq=readtable(fpath_fundq, opts);
    
    if isempty(fundq)
        loop_times(gvkey_i) = toc;
        [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', gvkey]);
        clc
        time_left_        
        continue        
    end
    
    % table secm is read from csv with options
    opts = gen_opts(fpath_secm, {'datadate'}, DateFormat, var_type);    
    secm=readtable(fpath_secm, opts);
    
    % make sure empty cshom is read as NaN...
    secm.cshom(secm.cshom==0) = NaN;
    
    % revtq, 'niq', 'nopiq'
    %vars_fundq = {'datadate', 'rdq', 'cshoq', 'atq', 'ltq', 'seqq', 'mkvaltq', 'cshtrq', 'saleq'};
    vars_fundq = {'datadate', 'rdq', 'cshoq', 'atq', 'ltq', 'mkvaltq'};
    fundq = fundq(:, vars_fundq);
    fundq.cshoq = fundq.cshoq*1e6;
    fundq.mkvaltq = fundq.mkvaltq*1e6;
    fundq.atq = fundq.atq*1e6;
    fundq.ltq = fundq.ltq*1e6;

    idx = isnat(fundq.rdq);
    fundq.rdq(idx) = eomdate(fundq.datadate(idx) + calmonths(1));
    fundq = sortrows(fundq, {'datadate', 'rdq'}, 'ascend');
    [~, idx] = unique(fundq.datadate);
    fundq = fundq(idx, :); % keep first report for any datadate
    fundq.datadate = fundq.rdq;
    fundq = removevars(fundq, 'rdq');
    fundq = sortrows(fundq, 'datadate');

    % retime fundq (robust merge, Keys: datadate)
    fundq = table2timetable(fundq);
    fundq = retime(fundq, 'daily', 'previous');
    fundq = timetable2table(fundq);

    fundq = renamevars(fundq, 'cshoq', 'cshoq_fundq');
    secm = renamevars(secm, 'cshoq', 'cshoq_secm');

    secm = outerjoin(secm, fundq(:, {'datadate', 'cshoq_fundq'}),...
        'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');

    idx = isnan(secm.cshoq_secm);
    if sum(idx) ~= size(secm, 1)
        idx_ = idx & secm.datadate < min(secm.datadate(~idx));
        idx_ = idx_ & secm.datadate > max(secm.datadate(~idx));
    else
        idx_ = idx;
    end

    secm.cshoq_secm(idx_) = secm.cshoq_fundq(idx_);
    secm = removevars(secm, 'cshoq_fundq');
    secm = renamevars(secm, 'cshoq_secm', 'cshoq');
    clear idx idx_

    idx = isnan(secm.cshom);
    if sum(idx) ~= size(secm, 1)
        idx_ = idx & secm.datadate < min(secm.datadate(~idx));
        idx_ = idx_ & secm.datadate > max(secm.datadate(~idx));
    else
        idx_ = idx;
    end
    

    secm.cshom(idx_) = secm.cshoq(idx_);
    secm = removevars(secm, 'cshoq');
    clear idx idx_
    

    if numel(glob) == 0 % north america
        secm.qunit = ones(size(secm, 1), 1);
    end
    
    secm = secm(secm.datadate.Year>=min_year, :);
    secm = secm(~isnat(secm.datadate), :);
    
    if size(secm, 1) < 24 || all(isnan(secm.prccm))
        continue
    end
    
    %{
    if numel(glob) == 0
        DateFormat = 'yyyy-MM-dd';
        opts = gen_opts(fpath_secm, {'datadate'}, DateFormat, var_type);
        secm=readtable(fpath_secm, opts);
        %secm=readtable(fpath_secm);
        secm = secm(secm.datadate.Year>=min_year, :);
        secm = secm(~isnat(secm.datadate), :);

        secm.year = secm.datadate.Year;
        secm.month = secm.datadate.Month;
        secm.year = secm.datadate.Year;
        secm.month = secm.datadate.Month;

        Keys = {'iid', 'year', 'month'};
        secm = outerjoin(secm, secm(:, [Keys, {'mkvalincl'}]), 'Keys', Keys, 'MergeKeys', true);
    end
    %}
    
    if all(isnan(secm.trfm))
        secm.trfm = ones(size(secm,1), 1);
    end
    secm.iid = cellfun(@(c) pad(c, 2, 'left', '0'), secm.iid, 'UniformOutput', false);
    if ismember('gvkey', secm.Properties.VariableNames)
        secm = removevars(secm, 'gvkey');
    end
    
    secm = secm(~isnan(secm.prccm),:);
    secm.cshom(secm.cshom == 0) = NaN;
    secm.prccm = secm.prccm./secm.qunit;
    secm = sortrows(secm, {'iid', 'datadate'});
    iids = unique(secm.iid); % list of unique iids
    
    % remove letters from iids
    iids_no_letter = cellfun(@(c) erase(c, 'C'), iids, 'UniformOutput', false);
    iids_no_letter = cellfun(@(c) erase(c, 'W'), iids_no_letter, 'UniformOutput', false);
    iids_to_remove = {};
    %iids_to_remove = [iids_to_remove; iids(str2double(iids_no_letter)>=90)'];
    if strcmp(Country, 'CAN')
        iids_to_remove = [iids_to_remove, iids(~contains(iids, 'C', 'IgnoreCase', true))];
    elseif ~strcmp(Country, 'CAN')
        iids_to_remove = [iids_to_remove, iids(contains(iids, 'C', 'IgnoreCase', true))];
    end
    
    
    secm = secm(~ismember(secm.iid, iids_to_remove), :);
    iids = unique(secm.iid); % list of unique iids
        
    if isempty(secm)
        continue
    end
    
    for iid_i=1:numel(iids)
        iid = iids{iid_i};
        idx_iid = ismember(secm.iid, iid);
        if all(isnan(secm.prccm(idx_iid))) || numel(unique(secm.prccm(idx_iid))) == 1
            % remove iids with defunct information
            try
                iids_to_remove = [iids_to_remove, {iid}];
            catch
                iids_to_remove = [iids_to_remove; {iid}]';
            end
            continue
        end
        
        % fill curcdm (previous)
        secm.curcdm(idx_iid) = fillmissing(secm.curcdm(idx_iid), 'previous');
        % fill epf (previous)
        secm.(epf(glob))(idx_iid) = fillmissing(secm.(epf(glob))(idx_iid), 'previous');
        
        % complete first observation for trfm
        if isnan(secm.trfm(find(idx_iid, 1))) && ~isnan(secm.trfm(find(idx_iid, 1)+1))
            secm.trfm(find(idx_iid, 1)) = secm.trfm(find(idx_iid, 1)+1);
            % note: secm rows sorted by 1. iid; 2. datadate
        end
        
        % it is incorrect to extrapolate cshoc in this way.
        % it includes dividends into mcap... too low mcap
        
        % extrapolate cshom via ajexm
        if any(~isnan(secm.cshom(idx_iid))) && any(isnan(secm.cshom(idx_iid)))
            fprintf('\nExtrapolating cshom via ajexm\n')
            % first and last non-NaN cshom
            aux_date=[min(secm.datadate(~isnan(secm.cshom) & idx_iid)),...
                 max(secm.datadate(~isnan(secm.cshom) & idx_iid))];

             % only extrapolte backwards
            secm.cshom(secm.datadate<aux_date(1) & idx_iid)=...
                 secm.cshom(secm.datadate == aux_date(1) & idx_iid).*...
                (secm.ajexm(secm.datadate == aux_date(1) & idx_iid)./...
                 secm.ajexm(secm.datadate <  aux_date(1) & idx_iid));
             %{
              % never used. likely incorrect.
             ajpm = secm.ajpm(secm.datadate<aux_date(1) & idx_iid);
             ajpm = [ajpm(1); ajpm];
             dajpm = ajpm(2:end)./ajpm(1:end-1);
             
             % undo dividend payment effect of ajexm
            secm.cshom(secm.datadate<aux_date(1) & idx_iid)=...
                secm.cshom(secm.datadate<aux_date(1) & idx_iid)./dajpm;
                 %}

        end
                
    end
    
    secm = secm(~ismember(secm.iid, iids_to_remove), :);
    iids = unique(secm.iid); % list of unique iids

    if numel(iids) == 1
        secm.(epf(glob)) = repmat({'Y'}, size(secm, 1), 1);
        % compute mcap on the single iid (may need to extrapolate cshom
        % treat this iid as the primary issue
    elseif numel(iids) > 1 %&& strcmp(lib, 'lib_g')
        % merge epf from sechistory if lib==g
        % head(secm(idx_iid & ~ismember(secm.iid, pri.itemvalue), :))
        for i=find(ismember(sechistory_i.item, 'EPF'))'
            idx = ismember(secm.iid, sechistory_i.iid(i)) &...
                secm.datadate >= sechistory_i.effdate(i);
            if ~isnat(sechistory_i.thrudate(i))
                idx = idx & secm.datadate <= sechistory_i.thrudate(i);
            end
            
            % make sure to not override data with sechistory (e.g. there's
            % an error in sechistory in 248384
            %idx = idx & ~ismember(secm.epf, {'Y', 'N'});
            
            % I found inconsistent behaviour. Is this WRDS issue or
            % something with S&P more generally? when I download 281740
            % data from WRDS I get inconsistent data with sechistory, and
            % in this case sechistory is correct!
            % 030024
            secm.(epf(glob))(idx) = repmat(sechistory_i.itemvalue(i), sum(idx), 1);
            %secm.epf(idx) = repmat(sechistory.itemvalue(i), sum(idx), 1);

        end

        idx = ~strcmp(secm.(epf(glob)), 'N');
        secm = secm(idx, :);
    
    end
    
    [g, g_] = findgroups(secm.iid);
    VariableTypes = varfun(@class,secm,'OutputFormat','cell');
    chars = strcmp(VariableTypes, 'cell');
    doubles = strcmp(VariableTypes, 'double');
    
    vars_retime = {'ajexm', 'cshom', 'prccm', 'trfm'}; % vars to operate retime on
    vars_all_other = setdiff(secm.Properties.VariableNames, vars_retime); % all other variables (excluding datadate)
    vars_all_other = setdiff(vars_all_other, {'datadate'}); % exclude datadate
    vars_fill_missing = vars_all_other(~contains(vars_all_other, 'ret')); % vars to fill with different method
    vars_fill_missing = setdiff(vars_fill_missing, {'GroupCount'});
    
    T3 = timetable;
    

    for i=1:max(g)
        idx = g == i;
        
        if days(max(diff(secm.datadate(idx)))) > inf
            T3 = [T3; table2timetable(secm(idx, :))];
            % don't retime if there are gaps of at least two weeks in iid
            fprintf('big gap')
            temp = secm.datadate(idx);
            plot(temp(2:end), days(diff(secm.datadate(idx))))
            
            continue
        end
        
        dates_iid = secm.datadate(idx);
        
        
        all_dates_iid_might_miss = secm.datadate(...
            secm.datadate > min(secm.datadate(idx)) &...
            secm.datadate < max(secm.datadate(idx)));
        newDates = union(dates_iid, all_dates_iid_might_miss);

        T1 = table2timetable(secm(idx, :));
        T2 = retime(T1(:, vars_retime), newDates, 'previous');
       
        T2 = outerjoin(T2, T1(:, vars_all_other), 'Keys', {'datadate'});
        for v=vars_fill_missing
            T2.(v{:}) = fillmissing(T2.(v{:}), 'nearest');
        end
        T3 = [T3; T2];
        
    end
    
    secm = timetable2table(T3);
    clear T1 T2 T3
    
    if isempty(secm)
        continue
    end
    
    
    % join fx    
    secm = innerjoin(secm, fx_daily,...
    'LeftKeys', {'datadate', 'curcdm'},...
    'RightKeys', {'datadate', 'curcdq'});
    
    % Convert closing price to USD
    secm.prccm_usd = secm.prccm ./ secm.currency_val;
    secm.volume_usd_iid = secm.cshtrm .* secm.prccm_usd;

    secm.mcap_iid=abs(secm.prccm_usd.*secm.cshom); % market value in USD, unadjusted at secuirty level
    secm.mcap_iid(secm.mcap_iid==0) = 1;
    secm.mcap_iid(isnan(secm.mcap_iid)) = 1;
        
    % compute market capitalization
    mcap = varfun(@(x) sum(x, 'omitnan'), secm, 'InputVariables', {'mcap_iid'},...
        'GroupingVariables', 'datadate'); % market value at company level
    mcap = renamevars(mcap, {'Fun_mcap_iid'}, {'mcap'});
    secm = innerjoin(secm, mcap, 'Keys', {'datadate'});
    secm.weight_iid = secm.mcap_iid./secm.mcap;
    

    % Compute adjusted price
    % ((prccm/ ajexm )* trfm )[ current ] /( prccm/ ajexm )* trfm ))[ prior time period ]-1)*100)
    secm.prcadj = secm.prccm./secm.ajexm; % get adj price for the purpose of rets
    secm.prcadj_trfm = secm.prcadj.*secm.trfm;
    
    % declare ret columns
    secm.ret_iid = ones(size(secm,1),1);
    secm.retm_iid = ones(size(secm,1),1);
    %secm.retm_iid_w = NaN(size(secm,1),1);
    
    % by iid, compute daily ret (retm_iid) and multiply by weight
    [gs, grps] = findgroups(secm.iid);
    for iid_i=1:numel(grps)
        iid = grps{iid_i};
        idx_iid = find(gs==iid_i);
        secm.ret_iid(idx_iid(2:end))=...
            secm.prcadj_trfm(idx_iid(2:end))./...
            secm.prcadj_trfm(idx_iid(1:end-1));
        
        num_bus_days = days252bus(...
            secm.datadate(idx_iid(1:end-1)),...
            secm.datadate(idx_iid(2:end)))';
        num_bus_days(num_bus_days==0)=1;  % avoid division by 0
        secm.retm_iid(idx_iid(2:end)) = secm.ret_iid(idx_iid(2:end)).^(...
            1./num_bus_days');

    end
    
    secm.retm_iid_w = secm.retm_iid .* secm.weight_iid;
    secm.volume_usd_iid_w = secm.volume_usd_iid .* secm.weight_iid;
    
    retm = varfun(@(x) sum(x, 'omitnan'), secm,...
        'InputVariables', {'retm_iid_w'},...
        'GroupingVariables', {'datadate'});
    retm = renamevars(retm, 'Fun_retm_iid_w', 'retm'); 
    
    % compute trading volume in USD
    volume = varfun(@(x) sum(x, 'omitnan'), secm,...
        'InputVariables', {'volume_usd_iid_w'},...
        'GroupingVariables', {'datadate'});
    volume = renamevars(volume, 'Fun_volume_usd_iid_w', 'volume'); 


    Keys = {'datadate', 'GroupCount'};
    output = join(mcap, retm, 'Keys', Keys);
    output = join(output, volume, 'Keys', Keys);
    %{
    output = outerjoin(output,...
        fundq(:, {'datadate', 'mkvaltq', 'atq', 'ltq', 'seqq', 'cshoq', 'cshtrq', 'saleq'}),...
        'Keys', Keys, 'MergeKeys', true, 'type', 'left');
    %}
    Keys = {'datadate'};
    output = outerjoin(output,...
        fundq(:, {'datadate', 'atq', 'ltq', 'mkvaltq'}),...
        'Keys', Keys, 'MergeKeys', true, 'type', 'left');

    output.mcap(output.mcap==1) = NaN;
    writetable(output, fpath_secm_proc)
    
    
    if plot_figure
        clf
        set(gcf, 'Position', [50 50 1700 1500])
        set(gcf, 'Color', 'w')
        if ~isempty(secm)
            [~, idx] = unique(secm.datadate);

        subplot(2,3,1)
        hold on
        plot(output.datadate, output.mcap, '-', 'LineWidth', 2)
        plot_group_by(secm, {'iid'}, 'datadate', 'mcap_iid')
        plot(output.datadate, output.mkvaltq, '--', 'LineWidth', 1.5)
        iids = unique(secm.iid);
        legend([iids', {'mcap', 'mkvaltq'}])
        title('Market cap.')
        ylabel('USD')
        end
        
        subplot(2,3,2)
        plot_group_by(secm, {'iid'}, 'datadate', 'prccm')
        xlabel('Date')
        ylabel('Local currency')
        title('prccm')

        subplot(2,3,3)
        plot_group_by(secm, {'iid'}, 'datadate', 'trfm')
        title('trfm')
        
        subplot(2,3,4)
        plot_group_by(secm, {'iid'}, 'datadate', 'cshom')
        ylabel('# of shares')
        title('cshom')
        
        
        subplot(2,3,5)
        plot_group_by(secm, {'iid'}, 'datadate', 'ajexm')
        title('ajexm')

        sgtitle([company_i.conm, ' gvkey: ', gvkey])
        
        export_fig(fpath_fig)
        
        
    end

    
    %secm = table2timetable(secm);
    

    loop_times(gvkey_i) = toc;
    %
    %current_time = now;
    %current_time = datetime(current_time,'ConvertFrom','datenum');
    clc
    [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', gvkey]);
    time_left_
    if query_mode
        return
    end




end


sound(end_sound, end_sound_Fs);

