%{
***************************************************************************
RISK FREE RATE
***************************************************************************
%}
TB_dir = '/mnt/D/data/ir/TBs';
fpath = fullfile(TB_dir, 'TB.csv');
opts=detectImportOptions(fpath);
for i=1:numel(opts.VariableNames)
    opts.VariableTypes{i} = var_type(opts.VariableNames{i}); 
end
TB=readtable(fpath, opts);
TB = TB(~isnan(TB.TDRETNUA),:);
TB.BUSD2MAT = days252bus(TB.CALDT, TB.MATDT);


func = @(TDRETNUA, BUSD2MAT) TDRETNUA(BUSD2MAT==min(BUSD2MAT));
TB=rowfun(func, TB,...
    'InputVariables', {'TDRETNUA', 'BUSD2MAT'},...
    'GroupingVariables', 'CALDT',...
    'OutputVariable', 'tdyld');

TB = table2timetable(removevars(TB, {'CALDT', 'GroupCount'}), 'RowTimes', TB.CALDT);
TB = retime(TB, min(TB.Time):caldays(1):max(TB.Time), 'previous');
TB = timetable2table(TB);
TB.Properties.VariableNames{1} = 'datadate';
save(fullfile(TB_dir, 'TB.mat'), 'TB')
%TB.BUSD2MAT = days252bus(TB.Time, TB.MATDT);

%TB.MATDAT = cellfun(@(c) datetime(c(1:8), 'InputFormat', 'yyyyMMdd'), TB.KYCRSPID, 'UniformOutput', false);
%TB.pseudo_date = TB.caldt+28;
