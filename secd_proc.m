
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

query_mode = false;
query = '015325';

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
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
secm_dir = fullfile(comp_dir, [glob, 'secm']);
fundq_dir = fullfile(comp_dir, [glob, 'fundq']);
fig_dir = fullfile(comp_dir, 'figs', [glob, 'secd_proc']);


addpath(fullfile(matlab_dir, 'my_functions'))
addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad'))


cd (root)
addpath('/home/oren/Documents/MATLAB/my_functions')

run('var_types.m') % Variables types mapping


run('load_company.m')
run('load_security.m')
run('load_sechistory.m')
run('load_currency.m')

%load('/media/oren/D/data/ir/TB.mat')
%run('TBill.m') % risk free rate table, TB
%run('identification.m') % load company/security specific variables (dicts)
fpath = '/media/oren/D/data/ir/DTB';

[R_f, R_f_metadata] = load_data_files('DTB3', fpath, 'FRED');

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

epf = containers.Map; % for tracking extrapolation of cshoc
epf('g_') = 'epf';
epf('') = 'mkvalincl';
%%
clc

run('var_types.m')

gvkey_i_start = 1;
fpath = '/home/oren/PycharmProjects/comp/gvkey_set_lctq.csv';
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvkeys = readtable(fpath, opts);
gvkeys = gvkeys{:,1};


%
gvkeys = dir(secd_dir);
gvkeys = gvkeys(~[gvkeys.isdir]);
gvkeys = gvkeys(~contains({gvkeys.name}, '#'));

gvkeys=cellfun(@(c) c(1:end-4), {gvkeys.name},...
    'UniformOutput', false);
gvkeys = gvkeys(~cellfun(@isempty, gvkeys));
%}

%{
fpath = '/home/oren/PycharmProjects/comp/redo_gvkeys_secd_proc.csv';
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


% 1:34514 V
for gvkey_i=flip(1:ngvkeys)
    
    gvkey_i_start=gvkey_i;
    tic % measure time
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
    fpath_secd_proc = fullfile(secd_proc_dir, sprintf('%s.csv', gvkey));
    fpath_fig = fullfile(fig_dir, [gvkey, '.png']);
    
    if exist(fpath_secd, 'file') ~=2
        fprintf('%s not present in secd.', gvkey)
        continue
    end
    
   
    if exist(fpath_secd_proc, 'file')==2 && ~redo_all && ~query_mode    
        continue
    end
    
    finfo = dir(fpath_secd);
    if finfo.bytes <= 292
        continue
    end
    
    if numel(glob) == 0 % north america
        finfo = dir(fpath_secm);
        if finfo.bytes <= 292
            continue
        end
    end
    
    
    % get data on all securities of company by gvkey
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
    
    % table secd is read from csv with options
    opts = gen_opts(fpath_secd, {'datadate'}, DateFormat, var_type);
    secd = readtable(fpath_secd, opts);
    
    % make sure empty cshom is read as NaN...
    secd.cshoc(secd.cshoc==0) = NaN;

    % revtq, 'niq', 'nopiq'
    %vars_fundq = {'datadate', 'rdq', 'cshoq', 'atq', 'ltq', 'seqq', 'mkvaltq', 'cshtrq', 'saleq'};
    vars_fundq = {'datadate', 'rdq', 'cshoq', 'atq', 'ltq', 'mkvaltq'};
    fundq = fundq(:, vars_fundq);
    fundq.cshoq = fundq.cshoq*1e6;
    fundq.mkvaltq = fundq.mkvaltq*1e6;
    fundq.atq = fundq.atq*1e6;
    fundq.ltq = fundq.ltq*1e6;

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

    secd = outerjoin(secd, fundq(:, {'datadate', 'cshoq'}),...
        'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');

    idx = isnan(secd.cshoc);
    if sum(idx) ~= size(secd, 1)
        idx_ = idx & secd.datadate < min(secd.datadate(~idx));
        idx_ = idx_ & secd.datadate > max(secd.datadate(~idx));
    else
        idx_ = idx;
    end

    secd.cshoc(idx_) = secd.cshoq(idx_);
    secd = removevars(secd, 'cshoq');
    clear idx idx_


        
    

    if numel(glob) == 0 % north america
        secd.qunit = ones(size(secd, 1), 1);
    end

    secd = secd(secd.datadate.Year>=min_year, :);
    secd = secd(~isnat(secd.datadate), :);
    
    if size(secd, 1) < 252*2 || all(isnan(secd.prccd))
        continue
    end
    
    if numel(glob) == 0
        DateFormat = 'yyyy-MM-dd';
        opts = gen_opts(fpath_secm, {'datadate'}, DateFormat, var_type);
        secm=readtable(fpath_secm, opts);
        %secm=readtable(fpath_secm);
        secm = secm(secm.datadate.Year>=min_year, :);
        secm = secm(~isnat(secm.datadate), :);

        secd.year = secd.datadate.Year;
        secd.month = secd.datadate.Month;
        secm.year = secm.datadate.Year;
        secm.month = secm.datadate.Month;

        Keys = {'iid', 'year', 'month'};
        secd = outerjoin(secd, secm(:, [Keys, {'mkvalincl'}]), 'Keys', Keys, 'MergeKeys', true);
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
    
    
    secd = secd(~ismember(secd.iid, iids_to_remove), :);
    iids = unique(secd.iid); % list of unique iids
        
    if isempty(secd)
        continue
    end
    
    for iid_i=1:numel(iids)
        iid = iids{iid_i};
        idx_iid = ismember(secd.iid, iid);
        if all(isnan(secd.prccd(idx_iid))) || numel(unique(secd.prccd(idx_iid))) == 1
            % remove iids with defunct information
            try
                iids_to_remove = [iids_to_remove, {iid}];
            catch
                iids_to_remove = [iids_to_remove; {iid}]';
            end
            continue
        end
        
        % fill curcdd (previous)
        secd.curcdd(idx_iid) = fillmissing(secd.curcdd(idx_iid), 'previous');
        % fill epf (previous)
        secd.(epf(glob))(idx_iid) = fillmissing(secd.(epf(glob))(idx_iid), 'previous');
        
        % complete first observation for trfd
        if isnan(secd.trfd(find(idx_iid, 1))) && ~isnan(secd.trfd(find(idx_iid, 1)+1))
            secd.trfd(find(idx_iid, 1)) = secd.trfd(find(idx_iid, 1)+1);
            % note: secd rows sorted by 1. iid; 2. datadate
        end
        
        % extrapolate cshoc via ajexdi
        if any(~isnan(secd.cshoc( idx_iid))) && any(isnan(secd.cshoc( idx_iid)))
            fprintf('Extrapolating cshoc via ajexdi\n')
            aux_date=[min(secd.datadate(~isnan(secd.cshoc)& idx_iid)),...
                 max(secd.datadate(~isnan(secd.cshoc) & idx_iid))];

            secd.cshoc(secd.datadate<aux_date(1) & idx_iid)=...
                 secd.cshoc( secd.datadate==aux_date(1) & idx_iid).*...
                (secd.ajexdi(secd.datadate==aux_date(1) & idx_iid)./...
                 secd.ajexdi(secd.datadate< aux_date(1) & idx_iid));

        end
                
    end
    
    
    secd = secd(~ismember(secd.iid, iids_to_remove), :);
    iids = unique(secd.iid); % list of unique iids

    if numel(iids) == 1
        secd.(epf(glob)) = repmat({'Y'}, size(secd, 1), 1);
        % compute mcap on the single iid (may need to extrapolate cshoc
        % treat this iid as the primary issue
    elseif numel(iids) > 1 %&& strcmp(lib, 'lib_g')
        % merge epf from sechistory if lib==g
        % head(secd(idx_iid & ~ismember(secd.iid, pri.itemvalue), :))
        for i=find(ismember(sechistory_i.item, 'EPF'))'
            idx = ismember(secd.iid, sechistory_i.iid(i)) &...
                secd.datadate >= sechistory_i.effdate(i);
            if ~isnat(sechistory_i.thrudate(i))
                idx = idx & secd.datadate <= sechistory_i.thrudate(i);
            end
            
            % make sure to not override data with sechistory (e.g. there's
            % an error in sechistory in 248384
            %idx = idx & ~ismember(secd.epf, {'Y', 'N'});
            
            % I found inconsistent behaviour. Is this WRDS issue or
            % something with S&P more generally? when I download 281740
            % data from WRDS I get inconsistent data with sechistory, and
            % in this case sechistory is correct!
            % 030024
            secd.(epf(glob))(idx) = repmat(sechistory_i.itemvalue(i), sum(idx), 1);
            %secd.epf(idx) = repmat(sechistory.itemvalue(i), sum(idx), 1);

        end

        idx = ~strcmp(secd.(epf(glob)), 'N');
        secd = secd(idx, :);
    
    end
    
    
    
    
    [g, g_] = findgroups(secd.iid);
    VariableTypes = varfun(@class,secd,'OutputFormat','cell');
    chars = strcmp(VariableTypes, 'cell');
    doubles = strcmp(VariableTypes, 'double');
    
    vars_retime = {'ajexdi', 'cshoc', 'prccd', 'trfd'}; % vars to operate retime on
    vars_all_other = setdiff(secd.Properties.VariableNames, vars_retime); % all other variables (excluding datadate)
    vars_all_other = setdiff(vars_all_other, {'datadate'}); % exclude datadate
    vars_fill_missing = vars_all_other(~contains(vars_all_other, 'ret')); % vars to fill with different method
    vars_fill_missing = setdiff(vars_fill_missing, {'GroupCount'});
    
    T3 = timetable;

    for i=1:max(g)
        idx = g == i;
        
        if days(max(diff(secd.datadate(idx)))) > inf
            T3 = [T3; table2timetable(secd(idx, :))];
            % don't retime if there are gaps of at least two weeks in iid
            fprintf('big gap')
            temp = secd.datadate(idx);
            plot(temp(2:end), days(diff(secd.datadate(idx))))
            
            continue
        end
        
        dates_iid = secd.datadate(idx);
        
        all_dates_iid_might_miss = secd.datadate(...
            secd.datadate > min(secd.datadate(idx)) &...
            secd.datadate < max(secd.datadate(idx)));
        newDates = union(dates_iid, all_dates_iid_might_miss);

        T1 = table2timetable(secd(idx, :));
        T2 = retime(T1(:, vars_retime), newDates, 'previous');

        T2 = outerjoin(T2, T1(:, vars_all_other), 'Keys', {'datadate'});
        for v=vars_fill_missing
            T2.(v{:}) = fillmissing(T2.(v{:}), 'nearest');
        end
        T3 = [T3; T2];
        
    end
    
    secd = timetable2table(T3);
    clear T1 T2 T3
    
    if isempty(secd)
        continue
    end
    
    % join fx    
    secd = innerjoin(secd, fx_daily,...
    'LeftKeys', {'datadate', 'curcdd'},...
    'RightKeys', {'datadate', 'curcdq'});
    
    % Convert closing price to USD
    secd.prccd_usd = secd.prccd ./ secd.currency_val;
    secd.volume_usd_iid = secd.cshtrd .* secd.prccd_usd;
    % Note: For American companies cshtrd is for the whole company?
    % this is what the docs say, but it seems to be incorrect

    secd.mcap_iid=abs(secd.prccd_usd.*secd.cshoc); % market value in USD, unadjusted at secuirty level
 
    % compute market capitalization
    mcap = varfun(@(x) sum(x, 'omitnan'), secd,...
        'InputVariables', {'mcap_iid'},...
        'GroupingVariables', 'datadate'); % market value at company level
    mcap = renamevars(mcap, {'Fun_mcap_iid'}, {'mcap'});
    secd = innerjoin(secd, mcap, 'Keys', {'datadate'});
    secd.weight_iid = secd.mcap_iid./secd.mcap;

    % Compute adjusted price
    % ((prccd/ ajexdi )* trfd )[ current ] /( prccd/ ajexdi )* trfd ))[ prior time period ]-1)*100)
    secd.prcadj = secd.prccd./secd.ajexdi; % get adj price for the purpose of rets
    secd.prcadj_trfd = secd.prcadj.*secd.trfd;
    
    % declare ret columns
    secd.ret_iid = ones(size(secd,1),1);
    secd.retd_iid = ones(size(secd,1),1);
    %secd.retd_iid_w = NaN(size(secd,1),1);
    
    % by iid, compute daily ret (retd_iid) and multiply by weight
    [gs, grps] = findgroups(secd.iid);
    for iid_i=1:numel(grps)
        iid = grps{iid_i};
        idx_iid = find(gs==iid_i);
        secd.ret_iid(idx_iid(2:end))=...
            secd.prcadj_trfd(idx_iid(2:end))./...
            secd.prcadj_trfd(idx_iid(1:end-1));
        
        num_bus_days = days252bus(...
            secd.datadate(idx_iid(1:end-1)),...
            secd.datadate(idx_iid(2:end)))';
        num_bus_days(num_bus_days==0)=1;  % avoid division by 0
        secd.retd_iid(idx_iid(2:end)) = secd.ret_iid(idx_iid(2:end)).^(...
            1./num_bus_days');

    end
    
    secd.retd_iid_w = secd.retd_iid .* secd.weight_iid;
    secd.volume_usd_iid_w = secd.volume_usd_iid .* secd.weight_iid;
    
    retd = varfun(@(x) sum(x, 'omitnan'), secd,...
        'InputVariables', {'retd_iid_w'},...
        'GroupingVariables', {'datadate'});
    retd = renamevars(retd, 'Fun_retd_iid_w', 'retd'); 

    % compute trading volume in USD
    volume = varfun(@(x) sum(x, 'omitnan'), secd,...
        'InputVariables', {'volume_usd_iid_w'},...
        'GroupingVariables', {'datadate'});
    volume = renamevars(volume, 'Fun_volume_usd_iid_w', 'volume'); 
    
    Keys = {'datadate', 'GroupCount'};
    output = join(mcap, retd, 'Keys', Keys);
    output = join(output, volume, 'Keys', Keys);

    Keys = {'datadate'};
    %{
    output = outerjoin(output,...
        fundq(:, {'datadate', 'mkvaltq', 'atq', 'ltq', 'seqq', 'cshoq', 'cshtrq', 'saleq'}),...
        'Keys', Keys, 'MergeKeys', true, 'type', 'left');
    %}
    output = outerjoin(output,...
        fundq(:, {'datadate', 'mkvaltq', 'atq', 'ltq'}),...
        'Keys', Keys, 'MergeKeys', true, 'type', 'left');
    
    %output.amihud = output.retd./output.volume;
    writetable(output, fpath_secd_proc)
    
    
    if plot_figure
        clf
        set(gcf, 'Position', [50 50 1700 1500])
        set(gcf, 'Color', 'w')
        if ~isempty(secd)
            [~, idx] = unique(secd.datadate);

        subplot(2,3,1)
        hold on
        plot(output.datadate, output.mcap, '-', 'LineWidth', 2)
        plot_group_by(secd, {'iid'}, 'datadate', 'mcap_iid')
        plot(output.datadate, output.mkvaltq, '--', 'LineWidth', 1.5)
        iids = unique(secd.iid);
        legend([iids', {'mcap', 'mkvaltq'}])
        title('Market cap.')
        ylabel('USD')
        end
        
        subplot(2,3,2)
        plot_group_by(secd, {'iid'}, 'datadate', 'prccd')
        xlabel('Date')
        ylabel('Local currency')
        title('prccd')

        subplot(2,3,3)
        plot_group_by(secd, {'iid'}, 'datadate', 'trfd')
        title('trfd')
        
        subplot(2,3,4)
        plot_group_by(secd, {'iid'}, 'datadate', 'cshoc')
        ylabel('# of shares')
        title('cshoc')
        
        
        subplot(2,3,5)
        plot_group_by(secd, {'iid'}, 'datadate', 'ajexdi')
        title('ajexdi')

        sgtitle([company_i.conm, ' gvkey: ', gvkey])
        
        export_fig(fpath_fig)
        
        
    end

    
    %secd = table2timetable(secd);
    

    loop_times(gvkey_i) = toc;
    %
    %current_time = now;
    %current_time = datetime(current_time,'ConvertFrom','datenum');
    
    [~, ~, time_left_] = time_left(gvkey_i, loop_times, false, ['gvkey: ', gvkey]);
    clc
    time_left_
    if query_mode
        return
    end




end


sound(end_sound, end_sound_Fs);

