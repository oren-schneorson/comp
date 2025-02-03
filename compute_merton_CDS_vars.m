function [T, ME] = compute_merton_CDS_vars(gvkey, Tickers, glob, comp_dir, CDS_dir, fundq_vars, F_formula, DGS1, Ccy, DocClause)
%COMPUTE_MERTON_VARS Computes Distance to Default (DD) and Default Likelihood Indicator (DLI) for firm by gvkey
%   This function computes distance to default and default likelihood indicator from company equity and balance sheet data.
%   The main formula it solves is the Black-Scholes-Merton (BSM), equation (2) in Vassalou and Xing (2004).
%   Inputs to the formula are: market capitalization (mcap), risk free rate (DGS1), standard deviation of returns (sigma_E)
%   and a measure of the face value of debt (F). The latter could be either total libabilities (ltq) or following
%   the literature: lctq+.5*lltq. The precise way to process variables for the calculation of F is given in F_formula
%   Lastly, to compute the value of firm assets, an iterative process is used (see comment below).

%%


n_fundq_vars = numel(fundq_vars);
libs = {'lib_g', 'lib_na'};


run('var_types.m')
ME = '';

if isempty(Ccy)
    Ccy = 'USD';
end

if isempty(DocClause)
    DocClause = 'CR14';
end

if strcmp(comp_dir, '')
	comp_dir = '/media/oren/D/data/comp';
end

if strcmp(CDS_dir, '')
	CDS_dir = '/media/oren/D/data/CDS';
end



if numel(fundq_vars) == 0
	fundq_vars = {'ltq'};
	F_formula = @(x) x;

end


if ~isa(F_formula, 'function_handle')
	F_formula = @(x) sum(x, 'omitnan');
end

% libraries
secd_proc_dir = fullfile(comp_dir, [glob, 'secd_proc']);
fundq_proc_dir = fullfile(comp_dir, [glob, 'fundq_proc']);
CDS_dir = fullfile(CDS_dir, 'bank_run_index', 'data', Ccy, DocClause);

CDS = table;
for Ticker = Tickers'
    fname = [char(Ticker), '.mat'];
    fpath = fullfile(CDS_dir, fname);

    if exist(fpath, 'file') == 2
        load(fpath)
        CDS = struct2table(s_);
        clear s_
    end
end

if isempty(CDS)
    T = table;
    ME = 'No CDS data';
    return
end
    
CDS = renamevars(CDS, 'Date', 'datadate');
CDS.PoD = CDS.P(:, 4); % PoD at 1 year Tenor
CDS.PoD(CDS.PoD == 0) = NaN;
CDS.PoD = fillmissing(CDS.PoD, 'previous');
CDS.DD = -norminv(CDS.PoD); % Distance to Default



fname = [gvkey, '.csv'];
fpath = fullfile(secd_proc_dir, fname);
T = readtable(fpath);
T.mcap = T.mcap/1e6; % mcap in millions

if isempty(T)
    ME = 'secd_proc empty';
    return
end


fpath = fullfile(fundq_proc_dir, fname);
fundq = readtable(fpath);
fundq.lltq = fundq.ltq-fundq.lctq;

if isempty(fundq)
    ME = 'fundq_proc empty';
    return
end

%{
gvkey
numel(unique(fundq.datadate))-numel(fundq.datadate)

[~, idx] = unique(fundq.datadate);
idx = ismember(1:numel(fundq.datadate), idx);
fundq(ismember(fundq.datadate, fundq.datadate(~idx)),:)
%}

% keep last observation which was reported.
fundq = fundq(~isnat(fundq.datadate), :);
[~, idx] = unique(fundq.datadate, 'last');

fundq = table2timetable(fundq(idx,:));
fundq = fillmissing(fundq, 'previous');
F = rowfun(F_formula, fundq, 'InputVariables', fundq_vars);
fundq.F = F{:, 1};
clear F

fundq = retime(fundq, 'daily');
fundq = fillmissing(fundq, 'previous');
fundq = timetable2table(fundq);


T = innerjoin(T, fundq(:, {'datadate', 'F'}));

if isempty(T)
    ME = 'fundq_proc and secd_proc do not overlap.';
    return
end

%{
unused...
% conditional volatility, one year backwards
T = [T, ...
    rowfun(@(datadate) std(T.retd(...
    T.datadate < datadate &...
    T.datadate >= datetime(datadate.Year-1, datadate.Month, datadate.Day)),...
    'omitnan'), T, 'InputVariables', {'datadate'}, 'OutputVariableNames', 'sigma_E')];
T.sigma_E = T.sigma_E * sqrt(252); % annualized
%}

% risk free rate: 3. 1-Year Treasury Constant Maturity Rate
if isa(DGS1, 'char') % in case DGS1 is not already loaded
    DGS1 = readtable(DGS1);
end

T = innerjoin(T, DGS1, 'Keys', 'datadate');
T = renamevars(T, 'DGS1', 'r');
T.r = log(1+T.r/100);

if isempty(T)
    ME = 'No overlap between fundq and DGS1';
    return
end

T = innerjoin(T, CDS(:, {'datadate', 'DD'}), 'Keys', 'datadate');

if isempty(T)
    ME = 'No overlap between secd_proc and CDS';
    return
end


%
idx = ~isnan(T{:, {'F', 'r', 'DD'}});
idx = idx & T{:, {'F', 'r', 'DD'}} ~= 0;
idx = all(idx, 2);
T = T(idx,:);

if isempty(T)
    ME = 'not enough information.';
    return
end
%}



%{
***************
Functions

%}

% Value of equity: Black-Scholes-Merton
E = @(V_, F_, r_, T_, d_1_, d_2_) V_ .* normcdf(d_1_) - exp(-r_ .* T_) .* F_ .* normcdf(d_2_);

d_1 = @(V_, F_, r_, T_, sigma_V_) ( log(V_./F_) + (r_ + .5 * sigma_V_.^2) .* T_ ) ./ ( sigma_V_ .* sqrt(T_) );

d_2 = @(V_, F_, r_, T_, sigma_V_) d_1(V_, F_, r_, T_, sigma_V_) - sigma_V_ .* sqrt(T_);

F_star = @(V_, DD_, T_, sigma_V_, mu_V_) V_./exp(sigma_V_.*sqrt(T_).*DD_-mu_V_+.5*sigma_V_.^2);


% Black-Scholes-Merton Formula, for a given day
%BSM = @(V_, sigma_V_, mu_V_, t) E(V_, F_star(V_, T.DD(t), 1, sigma_V_, mu_V_), T.r(t), 1, d_1(V_, F_star(V_, T.DD(t), 1, sigma_V_, mu_V_), T.r(t), 1, sigma_V_), d_2(V_, F_star(V_, T.DD(t), 1, sigma_V_, mu_V_), T.r(t), 1, sigma_V_)) - T.mcap(t);

% here I changed F_star to F outside d_1, d_2
BSM = @(V_, sigma_V_, mu_V_, t)...
    E(V_, T.F(t), T.r(t), 1,...
    d_1(V_, F_star(V_, T.DD(t), 1, sigma_V_, mu_V_), T.r(t), 1, sigma_V_),...
    d_2(V_, F_star(V_, T.DD(t), 1, sigma_V_, mu_V_), T.r(t), 1, sigma_V_)) - T.mcap(t);


%{
ITERATIVE PROCESS
*************************************************************
Vassalou and Xing, 2004: Default Risk in Equity Returns
To calculate sigma_V we adopt an iterative procedure. We use daily data from the
past 12 months to obtain an estimate of the volatility of equity sigma_E , which is
then used as an initial value for the estimation of sigma_V . Using the Black–Scholes
formula, and for each trading day of the past 12 months, we compute V using
E as the market value of equity of that day. In this manner, we obtain daily
values for V. We then compute the standard deviation of those V’s, which is
used as the value of sigma_V , for the next iteration. This procedure is repeated until
the values of sigma_V from two consecutive iterations converge.
*************************************************************
%}

opts = optimoptions('fsolve', 'Algorithm', 'Levenberg-Marquardt', 'Display', 'none');
%opts = optimoptions('fsolve', 'Algorithm', 'Levenberg-Marquardt', 'Display', 'iter');
diff_V = 1;

E0 = T.mcap;
F0 = T.F;
V0 = E0 + F0; % initial condition

precision = 1e-7; % precision of iterative process


if size(T, 1) < 100
    ME = 'not enough information.';
    return
end
 
try
    %%
    counter = 1;
    
    while diff_V > precision
        
        %{
        Note: extreme changes in V may affect the ability of the process to
        converge, even though those changes are both in V0 and V1. e.g. if
        V is very close to zero at a certain part of the time series, then
        even minute changes in v would dramatically shift sigma_V.
        
        Censoring these changes when computing sigma_V should take care of
        this.
        
        %}
        V1 = NaN(size(T, 1), 1);
        mu_V = ones(size(T,1), 1)*0.05;
        
        dV0 = [NaN; log(V0(2:end)./V0(1:end-1))];
        
        % set extreme values as NaN
        idx = dV0 < prctile(dV0,99) & dV0 > prctile(dV0,1);
        dV0(~idx) = NaN;

        % compute sigma_V_0 based on V0   
        % unconditional volatility
        sigma_V0 = std(dV0(idx), 'omitnan')*sqrt(252)*ones(size(T,1), 1);
        %{
        % conditional volatility
        sigma_V0 = arrayfun(@(datadate) std(dV0(...
            T.datadate < datadate &...
            T.datadate >= datadate-365), 'omitnan'), T.datadate);
        sigma_V0 = sigma_V0*sqrt(252); % annualize
        sigma_V0(isnan(sigma_V0)) = std(dV0(idx), 'omitnan')*sqrt(252);
        %}
        
        for t = 1:size(T,1)
            [V1(t), ~, ~, ~] = fsolve(@(V_) BSM(V_, sigma_V0(t), mu_V(t), t), V0(t), opts);
        end
                
        diff_V = sum(abs(V1./V0-1))/numel(V1); % condition for iteratation stop
        
        dV1 = log(V1(2:end)./V1(1:end-1));
        idx = dV1 < prctile(dV1,99) & dV1 > prctile(dV1,1);
        dV1(~idx) = NaN;
        %{
        sigma_V1 = arrayfun(@(datadate) std(dV1(...
            T.datadate < datadate &...
            T.datadate >= datadate-365), 'omitnan'), T.datadate);
        sigma_V1 = sigma_V1*sqrt(252); % annualize
        sigma_V1(isnan(sigma_V1)) = std(dV1(idx), 'omitnan')*sqrt(252);
        %}
        sigma_V1 = std(dV1, 'omitnan')*sqrt(252) * ones(size(T,1), 1);

        %
        clf
        subplot(2,2,1)
        hold on
        plot(T.datadate, V1)
        plot(T.datadate, V0)
        subplot(2,2,2)
        hold on
        plot(T.datadate, F_star(V1, T.DD, 1, sigma_V1, mu_V)./T.F)
        subplot(2,2,3)
        hold on
        plot(T.datadate, V1./V0-1)
        title(sprintf('%.10f', diff_V))
        subplot(2,2,4)
        hold on
        plot(T.datadate, T.F)
        yyaxis right
        plot(CDS.datadate, CDS.PoD)
        

        
        
        export_fig('./temp.png')
        %}
        

    
        V0 = .5*V0+.5*V1;
    
    counter = counter + 1;
    if counter > 100
        ME = 'algorithm did not converge';
        T = table;
        return
    end

    end
    
    %ln_V_F1 = T.DD.*sigma_V0-mu_V+.5*sigma_V0.^2;
    %F1 = V1 ./ exp(ln_V_F1);
    

    

    T.V = V1;
    T.F_star = F_star(V1, T.DD, 1, sigma_V1, mu_V);
    
    T.sigma_V = sigma_V1;
    T.mu_V = mu_V;
    

    vars_to_keep = {'datadate', 'mcap', 'F', 'F_star', 'V', 'mu_V', 'sigma_V', 'DD'};
    T = T(:, vars_to_keep);
catch ME
   rethrow(ME) 
   T = table; 
end
    
    

