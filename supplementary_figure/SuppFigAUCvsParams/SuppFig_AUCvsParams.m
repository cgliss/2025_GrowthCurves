%% Supplementary Figure: AUC vs Params
%% User inputs
% Getting user inputs common across all scripts
cd ../../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("AllGrowthCurves_ln_20250320.mat")
load("ABXmechColorMap.mat")
cd supplementary_figure/SuppFigAUCvsParams
%% Getting no drug fit using median 
inx = find((allCCData.uM == 0) & (allCCData.CellDilutionFactor == cellConcentration) & strcmp(allCCData.Strain,'MiNoLi wt'));

r = size(inx,1); % number of drug conc
c = length(time); % number of time points

opticalDensity = nan(r,c);
for i=1:size(allCCData.MedianData(inx),1)
    opticalDensity(i,1:length(allCCData.MedianData{inx(i)})) = allCCData.MedianData{inx(i)};
end
opticalDensity = opticalDensity(:,1:length(timeLimited));
medianOD = median(opticalDensity,'omitnan');
yBlanked = exp(medianOD)*A600;
stdOD = std(exp(opticalDensity)*A600, 'omitnan');
[control.fit_function, control.lag , control.growth_rate, control.max_load] = fit_bacterial_growth_Gompertz_v3(gompertz_model,timeLimited, medianOD, false, lower_bound{2}, upper_bound{2}, yBlanked', A600);

control.doubling_time = log(2)/control.growth_rate;

figure; hold on;
plot(timeLimited,exp(opticalDensity)*A600,'k');
errorbar(timeLimited,yBlanked, stdOD,'or','LineWidth',2);
plot(timeLimited,exp(control.fit_function(timeLimited))*A600,'-b','LineWidth',2);
xline(control.lag)
yline(exp(control.max_load)*A600)
grid on
box on
title('Gompertz fit to no drug')
xlabel('hours')
ylabel('OD (blank substracted)')

% Threshold for AUC
t100 = linspace(0, (lastTimePoint-1)*timeResolution, 100); % Time vector (0 to 20.667hrs)
growthCurve_control = control.fit_function(t100);
initial_auc = trapz(t100, exp(growthCurve_control)*A600);
auc50 = initial_auc * 0.5;
errorRates = 0.005;

%% Plot AUC change with parameters
vecLength = 20;
yCNTRL = exp(gompertz_model([control.max_load, control.growth_rate,control.lag], timeLimited))*A600;
AUCcntrl = trapz(timeLimited, yCNTRL);

figure
tiledlayout(1,3)
lag = 0:1:max(timeLimited);
gr= 0:0.5:10;
growth_rate = log(2)./gr;
maxload = 0:0.05:0.7;
maxload = log(maxload./A600);
nexttile; hold on;
sAUC = [];
for i =lag
yFitted = exp(gompertz_model([control.max_load, control.growth_rate,i], timeLimited))*A600;
AUC = trapz(timeLimited, yFitted);
scatter(i, AUC/AUCcntrl, 'k','filled')
sAUC = [sAUC, AUC/AUCcntrl];
end 
plot(lag, sAUC, '-')
xlim([0, max(timeLimited)])
ylim([0 1])
grid on; box on; axis square;
xlabel('lag (h)')
ylabel('AUC')
title('lag component')

nexttile; hold on;sAUC = [];
for i =growth_rate
yFitted = exp(gompertz_model([control.max_load,i, control.lag], timeLimited))*A600;
AUC = trapz(timeLimited, yFitted);
scatter(log(2)/i, AUC/AUCcntrl, 'k','filled')
sAUC = [sAUC, AUC/AUCcntrl];
end 
plot(log(2)./growth_rate, sAUC, '-')
xlim([0, 10.5])
ylim([0 1])
grid on; box on; axis square;
xlabel('generation time (h)')
ylabel('AUC')
title('growth component')

nexttile; hold on;sAUC = [];
for i =maxload
yFitted = exp(gompertz_model([i,control.growth_rate, control.lag], timeLimited))*A600;
AUC = trapz(timeLimited, yFitted);
scatter(-exp(i)*A600, AUC/AUCcntrl, 'k','filled')
sAUC = [sAUC, AUC/AUCcntrl];
end 
plot(-exp(maxload).*A600, sAUC, '-')
xlim([-0.71 0])
ylim([0 1])
xlabel('-max load (Abs)')
ylabel('AUC')
grid on; box on; axis square;
title('load component')