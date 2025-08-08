function [fit_function, lag, growth_rate, max_load] = fit_bacterial_growth_Gompertz_v3(gompertz_model, time_points, OD_values, tfDEBUG, lower_bounds, upper_bounds, yBlanked, A600)
% Wrapper function for fitting bacterial growth curves using a logistic model.
% The function returns the fit function, and the individual parameters:
% lag50, growth rate, max load, and lag05.
%
% Parameters:
%  - time_points: A vector of time points corresponding to OD measurements.
%  - OD_values: A vector of measured OD values corresponding to the time points.
%
% Returns:
%  - fit_function: A handle to the fitted logistic function.
%  - lag50: The time at which the bacterial OD reaches half of the max load.
%  - growth_rate: The growth rate parameter.
%  - max_load: The maximum bacterial load.
%  - lag05: The time at which the bacterial OD reaches 0.05 of the max load.

% p(1) is the maximum population
% p(2) is the growth rate
% p(3) is the lag: time between when a microbial population is transferred to a new habitat recovers and when a considerable cell division occurs

%% Estimate lag phase by fitting Gompertz model to blanked data
% Initial guess for fitting parameters
initial_paramsB = [max(yBlanked), log(2), mean(time_points)]; % [max_load, growth_rate, lag50]

% Bounds for the parameters
lower_boundsB = [0, 0, 0]; % No negative values allowed
upper_boundsB = [max(yBlanked), 2, max(time_points)]; % Change upperbound of max load to 1.2

% Set optimization options
options = optimoptions('lsqcurvefit', 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'off');

% Fit the Gompertz model to the data 
fitted_paramsB = lsqcurvefit(gompertz_model, initial_paramsB, time_points, yBlanked', lower_boundsB, upper_boundsB, options);

% Extract fitted parameters
max_loadB = fitted_paramsB(1);
growth_rateB = fitted_paramsB(2);
lag_guess = fitted_paramsB(3);

% Define the fitted Gompertz model
fit_functionB = @(t) max_loadB*exp(-exp(((exp(1)*growth_rateB)/max_loadB)*(lag_guess-t)+1));
%  if(tfDEBUG)
% figure; hold on;
% plot(time_points, yBlanked,'or');
% plot(time_points,fit_functionB(time_points),'-k');
% ylim([-0.05 0.8])
% xline(lag_guess)
%  end
%% Fit gompertz model to relative population (ln(N/N0))
% Initial guess for fitting parameters
initial_params = [max(OD_values), log(2), mean(time_points)]; % [max_load, growth_rate, lag50]

% Weighing exponential region 
weights = 1 ./ (1 + exp(-1 * (time_points - lag_guess))); % Weighting function
% figure; hold on;
% plot(time_points,weights.*OD_values)
% plot(time_points,OD_values, 'or')

options = optimoptions('lsqnonlin', 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'off');
gompertz_residuals = @(params, t, y, weights) weights .* (gompertz_model(params, t) - y);

fitted_params = lsqnonlin(@(params) gompertz_residuals(params, time_points, OD_values, weights), initial_params, lower_bounds, upper_bounds, options);

% Extract fitted parameters
max_load = fitted_params(1);
growth_rate = fitted_params(2);
lag = fitted_params(3);

% Define the fitted Gompertz model
fit_function = @(t) max_load*exp(-exp(((exp(1)*growth_rate)/max_load)*(lag-t)+1));

% Calculate lag05: the time at which the OD reaches lagParameter of the max loa
% target_OD = lagParameter * max_load;
% lag05 = fminbnd(@(t) abs(fit_function(t) - target_OD), min(time_points), max(time_points));

if(tfDEBUG)
    figure; 
    subplot(1,2,1)
    hold on;
    plot(time_points, OD_values.*weights,'or');
    plot(time_points,fit_function(time_points),'-k');
    ylim([-0.05 6])
    xline(lag_guess)
    yline(log(max_loadB/A600))
    axis square
    subplot(1,2,2)
    hold on;
    plot(time_points, exp(OD_values)*A600,'or');
    plot(time_points,exp(fit_function(time_points))*A600,'-k');
    ylim([-0.05 0.8])
    xline(lag_guess)
    yline(max_loadB)
    axis square
end
end