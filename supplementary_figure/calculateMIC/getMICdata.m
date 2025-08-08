%% Getting MIC of each drug based on which concentration is above the inhibitory concentration
cd ../../
load('AllGrowthCurves_Part2.mat.mat','allAUCdata')
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
cd supplementary_figure/calculateMIC
%%
drugs = unique(allAUCdata.DrugName);
fh = figure('color','w','Position', get(0, 'Screensize'));
ic=1;MICconcS = []; tis ={};firstConcs=[];
for i = 1:length(drugs)
    if i== 21
        fh = figure('color','w','Position', get(0, 'Screensize'));
        ic =1; 
    end 
    RIs = find(contains(allAUCdata.DrugName, drugs{i}) & allAUCdata.CellDilutionFactor==200);
    uniStrain = unique (allAUCdata{RIs, 'Strain'});
    for k1 = 1:length(uniStrain)
    subplot(5,5,ic)
     hold on;
    ic = ic+1; 
    RI = find(contains(allAUCdata.DrugName, drugs{i}) & allAUCdata.CellDilutionFactor==200 & contains(allAUCdata.Strain, uniStrain{k1}));
    AUC = allAUCdata{RI, "AUC"}; % getting AUC values: rows with numbers are ones that last OD reading<0.2 
    noGrowthInd = find((abs(AUC)-0)==min(abs(AUC)-0));
    dilutions = allAUCdata{RI, "uM"};
    %MICconc = max(dilutions(noGrowthInd)); 
    ncInd = RI(dilutions == 0); %no drug control index 
    ncAUC = allAUCdata{ncInd, "AUC"};
    [sDil,sind] = sort(dilutions,'ascend'); 
    legs = {};MICconc = []; firstConc = []; 
    fn=0; 
    for i1 = RI(sind)'
        aucV = allAUCdata{i1, "AUC"};
        ratio = aucV/ncAUC; 
        if ratio<0.05
            if isempty(MICconc)
                MICconc = allAUCdata{i1, "uM"};
                AUCmic = allAUCdata{i1, "AUC"};
                ratioMIC = ratio;
            end
            plot(timeLimited, allAUCdata{i1,"BlankedMedianData"}{1}(1:lastTimePoint), 'r', 'LineWidth',1.5)
        elseif aucV<ncAUC 
            if fn == 0
                firstConc = [firstConc; allAUCdata{i1, "uM"}]; 
                fn=1; 
                 plot(timeLimited, allAUCdata{i1,"BlankedMedianData"}{1}(1:lastTimePoint), 'b', 'LineWidth',1.5)
            else
            plot(timeLimited, allAUCdata{i1,"BlankedMedianData"}{1}(1:lastTimePoint), 'g', 'LineWidth',1.5)
            end
        else
             plot(timeLimited, allAUCdata{i1,"BlankedMedianData"}{1}(1:lastTimePoint), 'k', 'LineWidth',1.5)
        end 
        
        legs = [legs, sprintf('%0.3fuM; %0.2fAUC; ratio: %0.2f',  allAUCdata{i1, "uM"},  aucV, aucV/ncAUC)];
    end 
    grid on; box on; 
    if isempty(MICconc)
        title(sprintf('%s %s\n NO MIC', drugs{i}, uniStrain{k1}), 'Color','r')
    else
        title(sprintf('%s %s\nMIC= %.3f uM; Ratio to Cntrl=%0.2f', drugs{i},uniStrain{k1}, MICconc,ratioMIC))
    end 
    if isempty(MICconc)
    MICconcS = [MICconcS; 0];
    else
    MICconcS = [MICconcS; MICconc];
    end
    firstConcs = [firstConcs; firstConc];
    legend(legs,'Location', 'eastoutside', 'FontSize',8)
    xlim([0 max(timeLimited)])
    ylim([-0.2 1.2])
    xlabel('time(h)')
    ylabel('absorbance')
    tis = [tis, [drugs{i} ' ' uniStrain{k1}]];
    end 
end 
