
%{
merge2T.m, part of a package for TASE.

1. main.m
2. preamble.m
3. load_data.m
4. merge2T.m
5. get_indices.m

merges all the data loaded to a single table, T.

%}
clc


join_TS = false;
%{
% ** Quarterly reports (qr)
% * merge security data with reported data by quarter-year
% allow time for reported data to be published for information to be
% incorporated into market prices
information_set_shift = 0; % 0 by default

unq_gvkey = unique(T.gvkey);

fundq_ = table;

for gvkey_i = 1:ngvkeys
    clc
    fprintf('retiming fundq: %.2f%%\n', gvkey_i/ngvkeys*100)
    
    fname = gvkeys(gvkey_i).name;
    [~, gvkey, ~] = fileparts(fname);
    gvkey = str2double(char(gvkey));
    idx = fundq.gvkey == gvkey;
    fundq_i = fundq(idx, :);
    fundq_i = sortrows(fundq_i, {'datadate', 'rdq'}, {'descend', 'ascend'});
    [~, idx] = unique(fundq_i(:, {'rdq'}), 'rows');
    fundq_i = fundq_i(idx, :);
    
    fundq_i = removevars(fundq_i, 'datadate');
    fundq_i = renamevars(fundq_i, 'rdq', 'datadate');
    fundq_i = table2timetable(fundq_i);
    fundq_i = retime(fundq_i, 'daily', 'previous');
    fundq_i = timetable2table(fundq_i);
    if isempty(fundq_)
        fundq_ = fundq_i;
    else
        fundq_ = [fundq_; fundq_i];
    end
end

% join on Keys
Keys = {'gvkey', 'datadate'};
vars_fundq = fundq_.Properties.VariableNames;
vars_fundq = vars_fundq(~ismember(vars_fundq, Keys));

vars_T = T.Properties.VariableNames;
vars_T = vars_T(~ismember(vars_T, Keys));
vars_T = vars_T(~ismember(vars_T, vars_fundq));
% keep market information without reports, outerjoin...
T = outerjoin(T(:, [Keys, vars_T]), fundq_(:, [Keys, vars_fundq]),...
    'Keys', Keys, 'Type', 'left', 'MergeKeys', true);

%}
%{
% drop market information without reports
T = innerjoin(T, qr(:, [vars_qr, Keys]),...
    'Keys', Keys);
%}


% Book to market equity
T.b2m = (T.atq-T.ltq)./T.mcap;
%}

% add Risk free rate (R_f)
Keys = {'datadate'};
vars_R_f = {'DTB3'};
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(vars_R_f, Keys));

T = outerjoin(T(:, vars_T), R_f(:, [Keys, vars_R_f]), 'Keys', Keys,...
    'Type', 'left', 'MergeKeys', true);


if join_TS
    


vars_CPI = {'PI_SA_12'};
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(vars_CPI, Keys));
% add CPI variables
T = outerjoin(T(:, vars_T), CPI_D(:, [Keys, vars_CPI]), 'Keys', Keys,...
    'Type', 'left', 'MergeKeys', true);

    
% add market index and returns (incl realized inflation)
vars_mkt = mkt.Properties.VariableNames;
vars_T = ~ismember(T.Properties.VariableNames,...
    setdiff(vars_mkt, Keys));
T = outerjoin(T(:, vars_T), mkt, 'Keys', Keys,...
    'Type', 'left', 'MergeKeys', true);

% add HML index (b2m factor)
vars_T = ~ismember(T.Properties.VariableNames,...
    HML.Properties.VariableNames(2:end));
T = outerjoin(T(:, vars_T), HML, 'Keys', Keys,...
    'Type', 'left', 'MergeKeys', true);

% add SMB index (size factor)
vars_T = ~ismember(T.Properties.VariableNames,...
    SMB.Properties.VariableNames(2:end));
T = outerjoin(T(:, vars_T), SMB, 'Keys', Keys,...
    'Type', 'left', 'MergeKeys', true);



% ** expected inflation, ex_PI_[M/D]:
vars_ex_PI = ex_PI_TIPS.Properties.VariableNames(...
    ~ismember(ex_PI_TIPS.Properties.VariableNames, Keys));
vars_T = ~ismember(T.Properties.VariableNames, vars_ex_PI);
T = outerjoin(T(:, vars_T), ex_PI_D(:, [Keys, vars_ex_PI]),...
    'Keys', 'datadate', 'MergeKeys', true, 'Type', 'left');

end


T = renamevars(T, 'DTB3', 'R_f');



%}
%%

if strcmp(freq, 'D')
    varorder = {'gvkey', 'datadate', 'mcap', 'atq', 'ltq', 'b2m', 'retd', 'retrd', 'xretd', 'R_f',...
        'gsector', 'ggroup', 'gind'};
elseif strcmp(freq, 'M')
    varorder = {'gvkey', 'datadate', 'mcap', 'atq', 'ltq', 'b2m', 'retm', 'retrm', 'xretm', 'R_f',...
        'gsector', 'ggroup', 'gind'};
end    
    
    
varorder = [varorder, setdiff(T.Properties.VariableNames, varorder)];
T = T(:, varorder);
T = sortrows(T, {'gvkey', 'datadate'});
