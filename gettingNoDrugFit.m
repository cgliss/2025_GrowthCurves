%% Getting no drug fit 
inx = find((allCCData.uM == 0) & (allCCData.CellDilutionFactor == 200) & strcmp(allCCData.Strain,'MiNoLi wt'));

r = size(inx,1); % number of drug conc
c = length(time); % number of time points

opticalDensity = nan(r,c);
for i=1:size(allCCData.MedianData(inx),1)
    opticalDensity(i,1:length(allCCData.MedianData{inx(i)})) = allCCData.MedianData{inx(i)};
end

nonRelativeOD = exp(opticalDensity)*A600;

medianOD = min(opticalDensity,[], 'omitnan');%median(opticalDensity, 'omitnan');
stdOD = std(nonRelativeOD, 'omitnan');
[control.fit_function, control.lag , control.growth_rate, control.max_load] = fit_bacterial_growth_Gompertz(gompertz_model,time, medianOD, false, [0 0 0], [10 10 max(time)]);
control.doubling_time = log(2)/control.growth_rate;

figure; hold on;
plot(time,median(nonRelativeOD,'omitnan'),'k');
errorbar(time,median(nonRelativeOD,'omitnan'),stdOD,'or','LineWidth',2);
plot(time,exp(control.fit_function(time))*A600,'-b','LineWidth',2);
xline(control.lag)
yline(control.max_load)
ylim([-0.05 0.8])
grid on
box on
title('Gompertz fit to no drug')
xlabel('hours')
ylabel('OD (blank substracted)')