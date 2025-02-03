function [T, ME] = compute_merton_vars(gvkey, glob, comp_dir, fundq_vars, F_formula, DGS1)
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

if ~exist('username', 'var') || isempty(username)
    username = getenv("USER");
end
if ~exist('lib_data', 'var') || strcmp(lib_data, '')
    lib_data = fullfile('/media', username, 'D', 'data');
end

if ~exist('comp_dir', 'var') || strcmp(comp_dir, '')
    comp_dir = fullfile(lib_data, 'comp');
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
if isa(DGS1, 'char')
    DGS1 = readtable(DGS1);
end

T = innerjoin(T, DGS1, 'Keys', 'datadate');
T = renamevars(T, 'DGS1', 'r');
T.r = log(1+T.r/100);

if isempty(T)
    ME = 'No overlap between fundq and DGS1';
    return
end


%
idx = ~isnan(T{:, {'F', 'r'}});
idx = idx & T{:, {'F', 'r'}} ~= 0;
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

% mcap=E, so the non-linear equation to solve is...
% note: need to think about frequency.


%}

% Value of equity: Black-Scholes-Merton
E = @(V_, F_, r_, T_, d_1_, d_2_) V_ .* normcdf(d_1_) - exp(-r_ .* T_) .* F_ .* normcdf(d_2_);

d_1 = @(V_, F_, r_, T_, sigma_V_) ( log(V_./F_) + (r_ + .5 * sigma_V_.^2) .* T_ ) ./ ( sigma_V_ .* sqrt(T_) );

d_2 = @(V_, F_, r_, T_, sigma_V_) d_1(V_, F_, r_, T_, sigma_V_) - sigma_V_ .* sqrt(T_);



% Black-Scholes-Merton Formula, for a given day
BSM = @(V_,sigma_V_, t) E(V_, T.F(t), T.r(t), 1, d_1(V_, T.F(t), T.r(t), 1, sigma_V_), d_2(V_, T.F(t), T.r(t), 1, sigma_V_)) - T.mcap(t);


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
diff_sigma_V = 1;

V0 = T.mcap + T.F; % initial condition

precision = 1e-6; % precision of iterative process

r=.5;

 
try
    counter = 1;
    %%
    while diff_sigma_V > precision

        V1 = NaN(size(T, 1), 1);
        fval = NaN(size(T, 1), 1);

        %{
        Note: extreme changes in V may affect the ability of the process to
        converge, even though those changes are both in V0 and V1. e.g. if
        V is very close to zero at a certain part of the time series, then
        even minute changes in v would dramatically shift sigma_V.
        
        Censoring these changes when computing sigma_V should take care of
        this.
        
        %}

        
        dV0 = log(V0(2:end)./V0(1:end-1));
        idx = dV0 < prctile(dV0,99) & dV0 > prctile(dV0,1);
        dV0(~idx) = NaN;

        % compute sigma_V_0 based on V0
        % conditional volatility
        sigma_V_0 = arrayfun(@(datadate) std(dV0(...
            T.datadate < datadate &...
            T.datadate >= datadate-365), 'omitnan'), T.datadate);
        
        % fill with previous value of sigma_V0
        sigma_V_0 = fillmissing(sigma_V_0, 'previous');

        % unconditional volatility
        %sigma_V_0 = nanstd(dV0(idx));
        sigma_V_0 = sigma_V_0*sqrt(252); % annualize

        % for each t, solve equation (2) in Vassalou and Xing (2004) for V_, given sigma_V_0
        for t = 252:size(T,1)
            %
            t
            V0(t)
            sigma_V_0(t)
            %}
            [V1(t), fval(t), ~, ~] = fsolve(@(V_) BSM(V_, sigma_V_0(t), t), V0(t), opts);
            %[V1(t), fval(t), exit_flag, output] = fsolve(@(V_) BSM(V_, sigma_V_0, t), V0(t), opts);
            %[counter, t, V1(t), V0(t), fval(t), exit_flag] % output

        end
        

        % compute new sigma_V_1 from V1
        dV1 = log(V1(2:end)./V1(1:end-1));
        idx = dV1 < prctile(dV1,99) & dV1 > prctile(dV1,1);
        dV1(~idx) = NaN;
        
        % compute sigma_V_1 based on V1
        % conditional volatility
        sigma_V_1 = arrayfun(@(datadate) std(dV1(...
            T.datadate < datadate &...
            T.datadate >= datadate-365), 'omitnan'), T.datadate);
        
        % unconditional volatility
        %sigma_V_1 = std(dV1, 'omitnan');

        sigma_V_1 = sigma_V_1*sqrt(252); % annalize
        diff_sigma_V = sum(abs(sigma_V_1-sigma_V_0)); % iterate until

        V0 = r*V0+(1-r)*V1; % set new initial condition
        %{
        clf
        subplot(2,2,1)
        hold on
        plot(T.datadate, V1)
        plot(T.datadate, V0)
        subplot(2,2,2)
        hold on
        plot(T.datadate, V1-V0)
        subplot(2,2,3)
        hold on
        plot(T.datadate(2:end), V0(2:end)./V0(1:end-1)-1)
        plot(T.datadate(2:end), V1(2:end)./V1(1:end-1)-1)
        subplot(2,2,4)
        hold on
        plot(T.datadate, T.F)

        
        
        export_fig('./temp.png')
        %}
    counter = counter + 1;
    if counter > 50
        ME = 'algorithm did not converge';
        T = table;
        return
    end
    end
    
    
    T.V = V1;
    
    T.dV = [NaN; V1(2:end)./V1(1:end-1)]-1;
    mu_V= rowfun(@(datadate) nanmean(T.dV(T.datadate < datadate)), T, 'InputVariables', {'datadate'}, 'OutputVariableNames', 'mu_V');
    sigma_V= rowfun(@(datadate) nanstd(T.dV(T.datadate < datadate)), T, 'InputVariables', {'datadate'}, 'OutputVariableNames', 'sigma_V');
    T.mu_V = mu_V.mu_V;
    T.sigma_V = sigma_V.sigma_V;
    T.DD = (log(T.V./T.F) + (T.mu_V-.5*T.sigma_V.^2)*252)./(T.sigma_V*sqrt(252));

    vars_to_keep = {'datadate', 'mcap', 'F', 'V', 'dV', 'mu_V', 'sigma_V', 'DD'};
    T = T(:, vars_to_keep);
catch ME
   rethrow(ME)
   T = table; 
end
    
    

