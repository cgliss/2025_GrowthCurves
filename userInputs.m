function [timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder] = userInputs(~)
% userInputs: contains all user inputs to be used across scripts 
%   Defines Gompertz model, timepoints to be used and parameters for
%   calculating Gompertz model 
%% User inputs
timeResolution = 1/3;
lastPossibleTimePoint = 73; % from longest time period of all experiments
lastTimePoint = 63; % last timepoint to use for determine halfAUC
cellConcentration = 200;
halfAUC = 0.5;
colors = [58 134 255;255 190 11; 255 0 110; ]./255; %lag, generation time, max load

% time
time = 0:timeResolution:(lastPossibleTimePoint-1)*timeResolution;
timeLimited =  0:timeResolution:(lastTimePoint-1)*timeResolution; % limited to max timepoint across all expirements 

% Bounds for gompertz model found based on looking at maximum, median, min of control (no drug) fit from
% function: gettingNoDrugFit.m
% [max_load, growth_rate, lag50]
% Cell concentration:  40,200,1000
boundOrder = [40, 200, 1000];
upper_bound = {[6 2 max(time)],[5.8, 1.7, max(time)],[6 1.5, max(time)]}; % max_load<control, growth_rate<control
lower_bound = {[0 0 0.7],[0 0 1], [0 0 2.4]}; % lag>control

% Gompertz Growth Model (Zwietering)
gompertz_model = @(p, t) p(1) * exp(-exp(((exp(1)*p(2))/p(1)) * (p(3)-t)+1));

% p(1) is the maximum population
% p(2) is the growth rate
% p(3) is the lag: time between when a microbial population is transferred to a new habitat recovers and when a considerable cell division occurs

%% Getting theoretical absorbance for starting population
load("fittedODcurve100uLM9.mat")
initialConc = 1/200; % 1:200 from OD 1
% Based on formula: ODdiluted =  ((A600 - blank) - lm.Coefficients.Estimate(1))./lm.Coefficients.Estimate(2);
% without blank since everything will be blank subtracted 
A600 = ((initialConc*lm.Coefficients.Estimate(2)) + lm.Coefficients.Estimate(1));

end