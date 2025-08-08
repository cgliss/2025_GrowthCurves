function [timeFine,od,conc, drugHalfAUC, noDrugGompertz] = get_growthcurve_matching_halfauc_GompertzInterp_v7(allData,drugLabel,tfDEBUG,dilutionFactor,AUCfactor, Strain,timeResolution, timeLimited, gompertz_model, startColor, A600)

inx = find(strcmp(allData.DrugName,drugLabel) & (allData.CellDilutionFactor == dilutionFactor) & strcmp(allData.Strain,Strain));
concentration = allData.uM(inx);

% Making time match time period of the experiment 
time = 0:timeResolution:(length(allData.MedianData{inx(1)})-1)*timeResolution;

% Getting OD
r = size(concentration,1); % number of drug conc
c = length(time); % number of time points
opticalDensity = nan(r,c);
for i=1:size(allData.MedianData(inx),1)
    opticalDensity(i,:) = allData.MedianData{inx(i)}(1:length(time));
end
% removing timepoints that don't have OD readings
opticalDensity = opticalDensity(:,~isnan(opticalDensity(1,:)));

% getting halfAUC
noDrugIdx = inx(concentration==0);
noDrug = allData{noDrugIdx, 'MedianData'}{1};
noDrugFitFunction = allData.fit_function{noDrugIdx};
noDrugGompertz = {noDrugFitFunction, allData.lag(noDrugIdx), allData.growth_rate(noDrugIdx), allData.max_load(noDrugIdx)};

% NEED TO CALCULATE using timeLimited 
noDrugLimited = noDrug(1:length(timeLimited));
aucLowConcentration = trapz(timeLimited(~isnan(noDrugLimited)), exp(noDrugLimited(~isnan(noDrugLimited)))*A600); % AUC at lowest drug concentration using timeLimited!!! 

% remove no drug control 
tfZero = (concentration==0);
concentration = concentration(~tfZero);
opticalDensity = opticalDensity(~tfZero,:);
inx = inx(~tfZero,:);

% Getting fitted Gompertz Model
rConc = size(allData.MedianData(inx),1);
fitGompertz = cell(size(rConc,1));
lags = nan(rConc,1);
growth_rates = nan(rConc,1);
max_loads = nan(rConc,1);
for i=1:length(inx)
    % getting fit to Gompertz model
    fitGompertz{i} = allData.fit_function(inx(i));
    lags(i) = allData.lag(inx(i));
    growth_rates(i) = allData.growth_rate(inx(i));
    max_loads(i) = allData.max_load(inx(i));
end

% Removing no gompertz fit
nanIDX = isnan(lags); % removing concentrations that don't have gompetz model
concentration(nanIDX) = [];
lags(nanIDX) = [];
growth_rates(nanIDX) = [];
max_loads(nanIDX) =[];

% Getting finer concentrations 
timeFine = linspace(min(time), max(time), 100);
% if trapz(timeLimited, opticalDensity(min(concentration)==concentration,1:length(timeLimited)))<halfAUC
%      % If AUC(min(concentration)<halfAUC then need to extrapolate to
%      % concentration half step below minimum
%     a = mean([min(concentration)/(concentration(1)/concentration(2)), min(concentration)]); 
% else 
%     a = min(concentration);
% end 

a = min(concentration)/(concentration(1)/concentration(2));%mean([min(concentration)/(concentration(1)/concentration(2)), min(concentration)]); 
b = max(concentration)*(concentration(1)/concentration(2));
n = 1000;
t = linspace(0, 1, n); % Linear parameter from 0 to 1
concentrationFine = a + (b - a) * (exp(t * 3) - 1) / (exp(3) - 1); % Exponential mapping to get increase in gradual spacing as approach max concentration 


% Add ranges for fit to limit 
lagsFit = fit(concentration, lags, 'linear'); 
growth_ratesFit = fit(concentration, growth_rates, 'linear');
max_loadsFit = fit(concentration, max_loads, 'linear');
fitLags = feval(lagsFit,concentrationFine);
fitGrowth_rates = feval(growth_ratesFit,concentrationFine);
fitMax_loads = feval(max_loadsFit,concentrationFine);

% Plotting parameter vs concentration 
fig1 = figure;
tiledlayout(1,3)
nexttile; hold on; 
plot(lagsFit, concentration, lags, 'o')
xlabel('concentration (uM)')
ylabel('lag (h)')
grid on; box on;
nexttile; hold on; 
plot(growth_ratesFit, concentration, growth_rates, 'o')
xlabel('concentration (uM)')
ylabel('growth rate (h)')
grid on; box on;
nexttile; hold on; 
plot(max_loadsFit, concentration, max_loads, 'o')
xlabel('concentration (uM)')
ylabel('max load (abs)')
grid on; box on;
sgtitle(sprintf("%s %s", drugLabel, Strain));
saveas(fig1, sprintf("ParameterVConc_%s_%s.jpg",drugLabel, Strain))
close(fig1)

% Define the interpolated growth function
fitOD = nan(length(concentrationFine), length(timeFine));
p=nan(1,3);
%figure; hold on;
for j = 1:length(concentrationFine)
    % p(1) is the maximum population
    % p(2) is the growth rate
    % p(3) is the lag: time between when a microbial population is transferred to a new habitat recovers and when a considerable cell division occurs

    p(3) = fitLags(j);
    p(2) = fitGrowth_rates(j);
    p(1) = fitMax_loads(j);
    
    % Get OD 
    fitOD(j,:) = gompertz_model(p,timeFine);
   
    % Making those with negative paramters into 0 because can't have negative
    % parameters
     %plot(exp(fitOD(j,:)),'k')
    if sum(p<0)>=1
        fitOD(j,:) = zeros(1, length(timeFine));
       % plot(exp(fitOD(j,:)),'r')
    end 
end 


% Find concentration corresponding to half AUC
aucAtConcentration =nan(length(concentrationFine),1);
for j = 1:length(concentrationFine)
    [~, ind] = min(abs(timeFine - max(timeLimited)));
    aucAtConcentration(j) = trapz(timeFine(1:ind), exp(fitOD(j, 1:ind))*A600); % Iterate over concentration axis
end
% Normalizing AUC to no drug AUC
normAUCatConc = aucAtConcentration/aucLowConcentration;  
halfAUCidx = find(normAUCatConc<AUCfactor,1);
%% Send back full period of time 
od = fitOD(halfAUCidx,:);
conc = concentrationFine(halfAUCidx);
drugHalfAUC = normAUCatConc(halfAUCidx);

%% Plotting sample of curves 

% Generate the colormap
nColors = 12;
colorSpectrum = [linspace(startColor(1), 1, nColors)', ...
    linspace(startColor(2), 1, nColors)', ...
    linspace(startColor(3), 1, nColors)'];
colorSpectrum = flipud(colorSpectrum);
fh = figure; hold on;
timeLimitedFine = linspace(min(timeLimited), max(timeLimited), 100);


% Getting index of representative curves
aucFrac = linspace(0, max(aucAtConcentration),10);

predictedIdx = nan(length(aucFrac),1);
for j5 = 1:length(aucFrac)
    temp = find(normAUCatConc>aucFrac(j5),1, 'last'); 
    if ~isempty(temp)
        predictedIdx(j5) = temp;
    end
end

predictedIdx(isnan(predictedIdx)) =[];
% Plotting 10 representative curves
for j1 = 1:length(predictedIdx)
    plot(timeFine(1:ind), exp(fitOD(predictedIdx(j1), 1:ind))*A600, 'Color', colorSpectrum(j1,:),'LineWidth',2)
end
xline(fitLags(halfAUCidx),'r')
yline(fitMax_loads(halfAUCidx),'r')
plot(timeFine(1:ind), exp(fitOD(halfAUCidx, 1:ind))*A600, 'r','LineWidth',2)
xlabel('Time'); ylabel('Concentration'); 
ylim([-0.05 1]);
xlim([min(timeFine), timeFine(ind)])
grid on;
title(sprintf("%s %s", drugLabel, Strain));

if(tfDEBUG)
    saveas(fh, sprintf("CurveScape_%s_%s.jpg",drugLabel, Strain))
    close(fh);
elseif(~tfDEBUG)
    close(fh);
end

end