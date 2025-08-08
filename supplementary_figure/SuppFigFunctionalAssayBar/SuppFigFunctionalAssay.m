%% This script is to plot the results of the functional assay
cd ../../
load("FunctionalAssayData_2025-03-20.mat")
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
cd supplementary_figure/SuppFigFunctionalAssayBar/
%% Graphing growth curves + LFC
c=["#A2142F",'#C0C0C0'];
drugNames = FunctionalAssayResults.DrugName;
IC = FunctionalAssayResults.IC;
cntrl_cond6 = []; 
cntrl_cond24 = [];
drug6 = [];
drug24 = [];
counter=1;
uniStrain = "MiNoLi wt";
for g3 = 1:length(drugNames)
    % Plotting functional assay
    for g = find(contains(lfctable.DrugName,drugNames(g3)) & contains(lfctable.Strain,uniStrain) & lfctable.IC==IC(g3))'
        lfcparams = lfctable{g,"LFC"}; %[noneTP bacAUC noneAUC log2(bacMedAUC./noneMedAUC) bacMedAUC noneMedAUC];
        timeG = lfctable{g,"CollectionTimepoint_hr"};
        bacAUC = lfcparams{2};
        noneAUC = lfcparams{3}; 
        lfcCB = log2(mean(noneAUC)/mean(bacAUC));
        if timeG == 6
            cntrl_cond6 = [cntrl_cond6, lfcCB];
            drug6 = [drug6, drugNames(g3)];
        else
            cntrl_cond24 = [cntrl_cond24, lfcCB];
            drug24 = [drug24, drugNames(g3)];
        end
    end
end
%% ordered
figure; 
tiledlayout(2,1)
nexttile; 
hold on;
[sCC, sIND] = sort(cntrl_cond6); 
b=bar(drug6(sIND),sCC);
b.FaceColor = 'flat';
temp = logical(FunctionalAssayResults.DrugDegradation6);
b.CData(temp(sIND),:) = repmat([1 1 1],sum(temp),1);
b.CData(~temp(sIND),:) = repmat([0 0 0],sum(~temp),1);
ylabel('log2(mean(control AUC)/ mean(bacteriaAUC))')
box on; grid on; 
ylim([-5 1])
title('6 h')

nexttile; 
hold on;
[sCC, sIND] = sort(cntrl_cond24); 
b=bar(drug24(sIND),sCC);
b.FaceColor = 'flat';
temp = logical(FunctionalAssayResults.DrugDegradation24);
b.CData(temp(sIND),:) = repmat([1 1 1],sum(temp),1);
b.CData(~temp(sIND),:) = repmat([0 0 0],sum(~temp),1);
ylabel('log2(mean(control AUC)/ mean(bacteriaAUC))')
ylim([-5 1])
box on; grid on; 
title('24 h')
