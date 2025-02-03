
% load CPI data

clc

lib_cpi = fullfile(lib_data, 'cpi_usa');
series = 'CPIAUCSL';
[CPI_M, CPI_metadata] = load_data_files(series, lib_cpi, 'FRED');
CPI_M = renamevars(CPI_M, series, 'CPI_SA');

CPI_D = table2timetable(CPI_M);
CPI_D = retime(CPI_D, 'daily', 'previous');
CPI_D = timetable2table(CPI_D);



CPI_M.PI_SA = [NaN; CPI_M.CPI_SA(2:end)./CPI_M.CPI_SA(1:end-1)];
CPI_M.PI_SA_12 = [NaN(12,1); CPI_M.CPI_SA(13:end)./CPI_M.CPI_SA(1:end-12)];

CPI_D.PI_SA = [NaN; CPI_D.CPI_SA(2:end)./CPI_D.CPI_SA(1:end-1)];

aux=CPI_D.date-[CPI_D.date-calmonths(12)]';
[time_diff, idx] = min(abs(aux));
idx = idx';
time_diff = time_diff';

CPI_D.PI_SA_12 = CPI_D.CPI_SA./CPI_D.CPI_SA(idx);
CPI_D.PI_SA_12 = CPI_D.PI_SA_12.^(1+days(time_diff));


aux=CPI_D.date-[CPI_D.date-calmonths(1)]';
[time_diff, idx] = min(abs(aux));
idx = idx';
time_diff = time_diff';

CPI_D.PI_SA_1 = CPI_D.CPI_SA./CPI_D.CPI_SA(idx);
CPI_D.PI_SA_1 = CPI_D.PI_SA_1.^(1+days(time_diff));

CPI_D.PI_SA_lagged_1M = CPI_D.PI_SA_1(idx);


aux=CPI_D.date-[CPI_D.date-calmonths(2)]';
[time_diff, idx] = min(abs(aux));
idx = idx';
time_diff = time_diff';

CPI_D.PI_SA_2 = CPI_D.CPI_SA./CPI_D.CPI_SA(idx);
CPI_D.PI_SA_2 = CPI_D.PI_SA_2.^(1+days(time_diff));

CPI_D.PI_SA_lagged_2M = CPI_D.PI_SA_1(idx);

CPI_M.PI_SA_lagged_1M = [NaN; CPI_M.PI_SA(1:end-1)];
CPI_M.PI_SA_lagged_2M = [NaN; CPI_M.PI_SA_lagged_1M(1:end-1)];


idx = ~cellfun(@(c) isempty(regexp(c, '^PI_', 'once')), CPI_M.Properties.VariableNames);

for var = CPI_M.Properties.VariableNames(idx)
    CPI_M.(char(var)) = (CPI_M.(char(var))-1)*100;
end

idx = ~cellfun(@(c) isempty(regexp(c, '^PI_', 'once')), CPI_D.Properties.VariableNames);
for var = CPI_D.Properties.VariableNames(idx)
    CPI_D.(char(var)) = (CPI_D.(char(var))-1)*100;
end


CPI_D = renamevars(CPI_D, 'date', 'datadate');
CPI_M = renamevars(CPI_M, 'date', 'datadate');
