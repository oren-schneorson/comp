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
CDS_dir = fullfile('/media', username, 'D', 'data', 'CDS');
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

Ccy = 'USD';
DocClause = 'CR14';


fpath = fullfile(lib_data, 'ir', 'DGS1.csv');
DGS1 = readtable(fpath);
DGS1 = renamevars(DGS1, 'DATE', 'datadate');
DGS1 = table2timetable(DGS1);
DGS1 = retime(DGS1, 'daily');
DGS1 = fillmissing(DGS1, 'previous');
DGS1 = timetable2table(DGS1);

fpath = fullfile(comp_dir, [glob, 'company.csv']);
company = readtable(fpath);


fname = 'attributes.csv';
fpath = fullfile(CDS_dir, fname);
attr = readtable(fpath);
[~, idx] = unique(attr(:, {'Ticker', 'gvkey_0'}), 'rows');
attr = attr(idx, :);

clc
missing = {};
ME_table = table;

fundq_vars = {'lctq', 'lltq'};
%fundq_vars = {'ltq'};

merton_dir_fundq_vars = join(fundq_vars, '_');

if ismember({'lltq', 'lctq'}, fundq_vars)
    F_formula = @(lctq, lltq) lctq + 0.5*lltq;
elseif ismember({'ltq'}, fundq_vars)
    F_formula = @(ltq) ltq;
end
    
fig_dir = fullfile(comp_dir, 'figs', [glob, 'merton_CDS']);
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
merton_CDS_dir = fullfile(comp_dir, [glob, 'merton_CDS']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq']);


fpath = fullfile('/home', username, 'PycharmProjects', 'comp', 'gvkey_set_lctq.csv');
opts = detectImportOptions(fpath);
opts.VariableTypes = {'char'};
gvks = readtable(fpath, opts);
gvks = gvks{:,1};

d = dir(secd_proc_dir);
d = d(~[d.isdir]);
d = char({d.name}');
d = d(:, 1:end-4);
d = cellstr(d);
gvks = gvks(ismember(gvks, d));
clear d


%{
gvks = dir(secd_proc_dir);
gvks = gvks(~[gvks.isdir]);
gvks = gvks([gvks.bytes] > 30);
gvks = gvks(~cellfun(@(c) contains(c, '#'), {gvks.name}));

gvks = char({gvks.name}');
gvks = gvks(:,1:end-4);
gvks = cellstr(gvks);
%}

loop_times = NaN(size(gvks));
loop_times(1) = 0;
n = numel(gvks);
    
j = 1;
tot = 1;
% 4: 2541-->

% ?
%{
negative F
I found some misreporting:
for gvkey:
105460,
dlltq and ltq are negative in 2003Q1; I'v deleted data <= 2003Q1

066072
loq and ltq are negative in 2003-01-31;

%}
neg_liab = [797, 932, 1654, 1688, 1693, 2378, 2109, 2562, 3629];
neg_liab_gvks = {};

for gvk_i= 1:n
    clc
    fprintf('%.2f%%\n', gvk_i/numel(gvks)*100)
    gvk_i

    tic % measure time
    gvkey = gvks{gvk_i}; % get gvkey

    if gvk_i <= (j-1)*n/tot || gvk_i > j*n/tot
        %continue
    end
    if ~any(gvk_i==neg_liab)
        %continue
    end
    
    idx = attr.gvkey_0 == str2double(gvkey);
    if ~any(idx)
        continue
    elseif sum(idx) > 1
        %fprintf('%s has more than one Ticker attached to it.', gvkey)
        %return
        %continue
        
    end
    
    Tickers = attr.Ticker(idx);
    %Ticker = char(Ticker);
    
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');
    idx_company = company.gvkey==str2double(gvkey);
    loc = company{idx_company, 'loc'};
    if ~ismember(loc, {'USA', 'CAN'}) && numel(glob) == 0
        %continue
    end


    
    if ~strcmp(gvkey, query) && query_mode
        continue
    end
        
    fpath_merton = fullfile(merton_CDS_dir, [gvkey, '.csv']);
    fpath_fig = fullfile(fig_dir, [gvkey, '.png']);
    
        
    if exist(fpath_merton, 'file') && ~redo_all && ~query_mode
    	fprintf('%s already done\n', gvkey)
        continue
    end    
    
    [T, ME] = compute_merton_CDS_vars(gvkey, Tickers, glob, comp_dir, CDS_dir, fundq_vars, F_formula, DGS1, Ccy, DocClause);
 
    if ~strcmp(ME, '')
        missing = [missing; {gvkey}];

        ME_table = [ME_table; table({gvkey}, {ME}, 'VariableNames', {'gvkey', 'ME'})];
        fprintf('%s error. %.2f\n', gvkey, gvk_i/n*100)
        continue
    elseif size(T,1) < 10
        ME = 'Too few valid observations';
        ME_table = [ME_table; table({gvkey}, {ME}, 'VariableNames', {'gvkey', 'ME'})];
        fprintf('%s error. %.2f\n', gvkey, gvk_i/n*100)
        continue
        
    end
    
    
    writetable(T, fpath_merton)

    
    if plot_figure
        clf
        %set(gcf, 'Position', [50 50 1700 1500])
        set(gcf, 'Color', 'w')

        subplot(1,3,1)
        plot(T.datadate, (T.F_star-T.F)./T.F)
        yyaxis right
        plot(T.datadate, normcdf(-T.DD))
        leg = {'F^*/F-1', 'PoD'};
        legend(leg, 'Location', 'North')
        title('Difference between DD implied equity and mcap')
        xlabel('Date')
        %ylabel('mm$')

        subplot(1,3,2)
        hold on
        area(T.datadate, T{:, {'F', 'mcap'}})
        plot(T.datadate, T{:, {'F_star'}})
        %title('F, mcap and E')
        leg = {'F', 'mcap', 'F^*'};
        legend(leg, 'Location', 'North')

        subplot(1,3,3)
        plot(T.datadate, T.DD)
        leg = {'(mcap-E)/F'};
        legend(leg, 'Location', 'North')
        yyaxis right
        plot(T.datadate, normcdf(-T.DD))
        leg = {'DD', 'PoD'};
        legend(leg, 'Location', 'North')
        title('DD and PoD')
        xlabel('Date')


        sgtitle([company.conm{idx_company}, '--', [' gvkey: ', gvkey]])
        
        export_fig(fpath_fig)
        
        
        
    end
   

    loop_times(gvk_i) = toc;
    
    
    %[~, ~] = time_left(gvk_i, loop_times, true, ['gvkey: ', gvkey]);
    
    if query_mode
        ME
        return
    end

end



fpath = fullfile(comp_dir, [glob, 'errors'], 'merton_error.csv');
writetable(ME_table, fpath);



sound(end_sound, end_sound_Fs);

return





















