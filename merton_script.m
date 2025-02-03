clear; clc


%{
The first step in implementing the Merton DD model is to estimate sigma_E from
either historical stock returns data or from option-implied volatility data.
%}

addpath('/home/u70o/Documents/MATLAB/altmany-export_fig-410f0ad')
addpath('/home/u70o/Documents/MATLAB/my_functions')

lib_data = '/media/u70o/D/data_temp4backup';
glob = '';
comp_dir = fullfile(lib_data, 'comp');
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
fig_dir = fullfile(comp_dir, 'figs', [glob, 'merton']);
merton_dir = fullfile(comp_dir, [glob, 'merton']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq']);

root = '/home/u70o/Documents/MATLAB/comp';
cd(root)


load handel.mat;
end_sound = y;
end_sound_Fs = Fs;

%%

clc

redo_all = true;
redo_query = true;
plot_figure = true;

query_mode = true;
%query = '015634';
%query = '015532';
%query = '028256';
%query = '011366'; % done
%query = '009340'; % done

% missing from ltq_lctq sample
%query = '013135'; 

%query = '163999'; % missing secd && secm
%query = '011628'; % missing secd && secm
%query = '154257'; % missing secd && secm
%query = '007864'; % missing secd && secm

%query = '003874'; % missing secd, exists secm
%query = '004478'; % missing secd, exists secm



query = '157373';
%fpath_queries = '/home/u70o/PycharmProjects/comp/gvkey_lists/queries.csv';
fpath_queries = '/home/u70o/Documents/MATLAB/measuring_rollover_risk/err1.xlsx';
opts = detectImportOptions(fpath_queries);
opts.VariableTypes = {'char'};
queries = readtable(fpath_queries, opts);




fpath = '/media/u70o/D/data/ir/DGS1.csv';
DGS1 = readtable(fpath);
DGS1 = renamevars(DGS1, 'DATE', 'datadate');
DGS1 = table2timetable(DGS1);
DGS1 = retime(DGS1, 'daily');
DGS1 = fillmissing(DGS1, 'previous');
DGS1 = timetable2table(DGS1);

fpath = fullfile(comp_dir, [glob, 'company.csv']);
company = readtable(fpath);


clc
missing = {};
ME_table = table;

%fundq_vars = {'lctq', 'lltq'};
fundq_vars = {'ltq'};

merton_dir_fundq_vars = join(fundq_vars, '_');

if ismember({'lltq', 'lctq'}, fundq_vars)
    F_formula = @(lctq, lltq) lctq + 0.5*lltq;
elseif ismember({'ltq'}, fundq_vars)
    F_formula = @(ltq) ltq;
end
    


fpath = '/home/u70o/PycharmProjects/comp/gvkey_set_lctq.csv';
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

if query_mode
    gvks = queries.gvkey;
    %gvks = {'005894'};
end


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
neg_liab_gvks = {};

for gvk_i= 1:n
%for gvk_i= flip(1:n)
    tic % measure time
    gvkey = gvks{gvk_i}; % get gvkey
    gvkey = pad(gvkey, 6, 'left', '0');


    if gvk_i <= (j-1)*n/tot || gvk_i > j*n/tot
        continue
    end
    
    
    
    %gvkey = pad(num2str(gvkey), 6, 'Left', '0');
    idx_company = company.gvkey==str2double(gvkey);
    loc = company{idx_company, 'loc'};

    if sum(idx_company) == 0
        fprintf('Missing company information. %s\n', gvkey)
        continue
    end

    if ~ismember(loc, {'USA', 'CAN'}) && numel(glob) == 0
        %continue
    end
        
    fpath_secd_proc = fullfile(secd_proc_dir, [gvkey, '.csv']);
    fpath_merton = fullfile(merton_dir, [gvkey, '.csv']);
    fpath_fig = fullfile(fig_dir, [gvkey, '.png']);

    if exist(fpath_secd_proc, 'file') ~= 2
        fprintf('Missing secd_proc - %s \n', gvkey)
        continue
    end
    
        
    if exist(fpath_merton, 'file') && ~redo_all && ~query_mode
    	fprintf('%s already done\n', gvkey)
        continue
    elseif query_mode && redo_query
    	fprintf('%s already done\n', gvkey)
        continue
    end    
    
    [T, ME] = compute_merton_vars(gvkey, glob, comp_dir, fundq_vars, F_formula, DGS1);
    
    
    if ~strcmp(ME, '')
        missing = [missing; {gvkey}];

        ME_table = [ME_table; table({gvkey}, {ME}, 'VariableNames', {'gvkey', 'ME'})];
        fprintf('%s error. %.2f\n', gvkey, gvk_i/n*100)
        continue
    end

    writetable(T, fpath_merton)

    
    if plot_figure
        clf
        set(gcf, 'Position', [50 50 1700 1500])
        set(gcf, 'Color', 'w')

        subplot(1,2,1)
        hold on
        plot(T.datadate, T.V)
        plot(T.datadate, T.F)
        plot(T.datadate, T.mcap)
        xlabel('Date')
        ylabel('mm$')
        leg = {'Firm value', 'Face value of debt', 'Market cap.'};
        legend(leg, 'Location', 'North')

        subplot(1,2,2)
        idx = 60;
        plot(T.datadate(idx:end), T.DD(idx:end))
        title('Distance to default')
        leg = {'DD'};
        legend(leg, 'Location', 'North')


        sgtitle([company.conm(idx_company), ' gvkey: ', gvkey])
        
        export_fig(fpath_fig)
        
        
    end

    loop_times(gvk_i) = toc;
    
    
    %[~, ~] = time_left(gvk_i, loop_times, true, ['gvkey: ', gvkey]);
    

end

sound(end_sound, end_sound_Fs);


return
fpath = fullfile(comp_dir, [glob, 'errors'], 'merton_error.csv');
writetable(ME_table, fpath);




return





















