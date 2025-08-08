%% Figure 3: Heatmap + Triangle
%% User inputs
% Getting user inputs common across all scripts
cd ../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("AllGrowthCurves_ln_20250320.mat")
load("ABXmechColorMap.mat")
load("FunctionalAssayData_2025-03-20.mat")
cd figure3
errorRates = 0.005;
%% Fitting Gompertz fit for the cell concentration of interest for WT
% Getting drug labels
drugLabels = unique(allCCData.DrugName);

myMechs =[];
myType =[];
for dl1 = 1:length(drugLabels)
    RI= find(contains(allCCData.DrugName,drugLabels{dl1}));
    typeA = allCCData{RI(1), "Type"};
    if strcmp(typeA,'nonABX')
        mech = "nonABX";
    else
        mech = allCCData{RI(1), "Mechanism"};
    end 
    myMechs = [myMechs, mech];
    myType = [myType, allCCData{RI(1), "Type"}];
end
myMechs(contains(myType, 'nonABX')) = 'nonABX';

fh = figure('color','w','Position', get(0, 'Screensize'), "Visible", 'on'); hold on;
tiledlayout(5,10)
lags = []; growth_rates = []; max_loads = []; 
startColor = [0 0 0]; storeCG = [];
for i= 1:length(drugLabels)
    % Find halfAUC using 20.33 hrs: by using fitted Gompertz Model for each
    % concentration and interpolating Gompertz parameters for
    % concentrations to find one corresponding to halfAUC
    % fitted Gompertz Model from "AnalyzeDiffCellConc_Gompertz.m" in "AllGrowthCurves_0.05lag_Set1to35.mat"
    [timeVector,od,conc, aucAtConcentration, controlGompertz] = get_growthcurve_matching_halfauc_GompertzInterp_v7(allCCData,drugLabels{i},true, cellConcentration, halfAUC,"MiNoLi wt", timeResolution, timeLimited, gompertz_model, startColor, A600);
    storeCG = [storeCG; controlGompertz{2:end} ];
    if abs(aucAtConcentration-halfAUC)>errorRates
        fprintf('ERROR- AUC50 is not close enough')
    end
    % Getting Gompertz model of halfAUC
    boundIdx = boundOrder == cellConcentration;
    yBlanked = exp(od)*A600;
    [fit_function, lags(i), growth_rates(i), max_loads(i)] = fit_bacterial_growth_Gompertz_v3(gompertz_model,timeVector, od, false, lower_bound{boundIdx}, upper_bound{boundIdx}, yBlanked', A600);

    % Comparing to controlGompertz and making it the no drug control Gompertz parameter if out of bounds
    % controlGompertz = {fit_function, lag, growth_rate, max_load}
    tColor = 'k';
    nexttile
    hold on;
    
    if lags(i)<controlGompertz{2} % less then no drug lag is reset to no drug lag
        lags(i) = controlGompertz{2};
        tColor = 'r';
    end
    if growth_rates(i)>controlGompertz{3} % more than no drug GR is reset to no drug GR
        growth_rates(i) = controlGompertz{3};
         tColor = 'g';
    end
    if max_loads(i)>controlGompertz{4} % more than no drug load is reset to no drug load
        max_loads(i)=controlGompertz{4};
        tColor = 'b';
    end
     if strcmp(myType{i}, 'nonABX')
        startColor =[0 0 0];
    else
        startColor = ABXmech.Colors{strcmp(ABXmech.Mechanism,myMechs{i})};
    end 
    idx = find(contains(allCCData.DrugName, drugLabels{i}) & contains(allCCData.Strain, "MiNoLi wt") & allCCData.CellDilutionFactor ==cellConcentration);
    %yFitted = exp(fit_function(timeLimited))*A600;
    yFitted = exp(gompertz_model([max_loads(i), growth_rates(i),lags(i)], timeLimited))*A600;
    plot(timeLimited, yFitted, 'Color',startColor, 'LineWidth',2)
    fill([timeLimited fliplr(timeLimited)], [yFitted zeros(1,length(timeLimited))] ,startColor, 'LineWidth',2, 'FaceAlpha',0.5, 'EdgeColor','none');
    %title(sprintf("%s\nnormAUC= %.2f", drugLabels{i},aucAtConcentration), 'Color',tColor)
    title(sprintf("%s %.4f", drugLabels{i},aucAtConcentration-0.5), 'Color',tColor)
    grid on; box on;
    ylim([0 0.8])
    xlim([0 max(timeLimited)])
end


% Fitting Gompertz fit for the cell concentration of interest for Cassette
lags_CC = []; growth_rates_CC = []; max_loads_CC = [];ccLabels = {};
ccDrugLabels = unique(cassetteAllData.DrugName);
myind = 1;
startColor =[0 0 0];
for i= 1:length(ccDrugLabels)
    Strains = unique(cassetteAllData{contains(cassetteAllData.DrugName,ccDrugLabels{i}),"Strain"});
    for i3 = 1:length(Strains)
        % Find halfAUC using 20.33 hrs: by using fitted Gompertz Model for each
        % concentration and interpolating Gompertz parameters for
        % concentrations to find one corresponding to halfAUC
        % fitted Gompertz Model from "AnalyzeDiffCellConc_Gompertz.m" in "AllGrowthCurves_0.05lag_Set1to35.mat"
         [timeVector,od,conc, aucAtConcentration, controlGompertz] = get_growthcurve_matching_halfauc_GompertzInterp_v7(cassetteAllData,ccDrugLabels{i},true, cellConcentration, halfAUC,Strains{i3}, timeResolution, timeLimited, gompertz_model, startColor,A600);
         if abs(aucAtConcentration-halfAUC)>errorRates
             fprintf('ERROR- AUC50 is not close enough')
         end
        if ~isempty(od)
            % Getting Gompertz model of halfAUC
            boundIdx = boundOrder == cellConcentration;
            yBlanked = exp(od)*A600;
            [fit_function, lags_CC(myind), growth_rates_CC(myind), max_loads_CC(myind)] = fit_bacterial_growth_Gompertz_v3(gompertz_model,timeVector, od, false, lower_bound{boundIdx}, upper_bound{boundIdx}, yBlanked', A600);

            ccLabels = [ccLabels; {ccDrugLabels{i} Strains{i3}}];

            % Comparing to controlGompertz and making it the no drug control Gompertz parameter if out of bounds
            % controlGompertz = {fit_function, lag, growth_rate, max_load}
            tColor = 'k';
            if lags_CC(myind)<controlGompertz{2} % less then no drug lag is reset to no drug lag
                lags_CC(myind) = controlGompertz{2};
                tColor = 'r';
            end
            if growth_rates_CC(myind)>controlGompertz{3} % more than no drug GR is reset to no drug GR
                growth_rates_CC(myind) = controlGompertz{3};
                tColor = 'g';
            end
            if max_loads_CC(myind)>controlGompertz{4} % more than no drug load is reset to no drug load
                max_loads_CC(myind)=controlGompertz{4};
                tColor = 'b';
            end

            nexttile
            hold on;
            idx = find(contains(cassetteAllData.DrugName, ccDrugLabels{i}) & contains(cassetteAllData.Strain, Strains{i3}) & cassetteAllData.CellDilutionFactor ==cellConcentration);
            yFitted = exp(fit_function(timeLimited))*A600;
            plot(timeLimited, yFitted, 'Color',startColor, 'LineWidth',2)
            fill([timeLimited fliplr(timeLimited)], [yFitted zeros(1,length(timeLimited))] ,startColor, 'LineWidth',2, 'FaceAlpha',0.5, 'EdgeColor','none');
            %title(sprintf("%s: %s\nnormAUC = %.2f",ccDrugLabels{i}, Strains{i3}, aucAtConcentration), 'Color', tColor)
            title(sprintf("%s: %s %.4f",ccDrugLabels{i}, Strains{i3}, aucAtConcentration-0.5), 'Color', tColor)
            grid on; box on;
            ylim([0 0.8])
            xlim([0 max(timeLimited)])

            myind = myind +1;
        end
    end
end
sgtitle('Red <control lag; Green >control growth rate; Blue > control max load')
gen_times_CC = log(2)./growth_rates_CC;
gen_times = log(2)./growth_rates;

% Assuming fig is the handle to the figure you are trying to save...
%set the figure's Renderer property to 'painters' before saving as EPS/SVG:
%set(fh, 'Renderer', 'painters');
% % % If saving as an EPS file, in R2020a and later, you can use the exportgraphics command and specify the 'ContentType' as 'vector':
%exportgraphics(fh, 'SupplementaryFig2.eps', 'ContentType', 'vector');
%clearvars fh
%saveas(fh,'GompertzFit.jpg')
%% No Drug average gompertz
aCG = mean(storeCG);
fprintf('\nlag = %0.3f +/- std %0.3f ', aCG(1), std(storeCG(:,1)))
fprintf('\ngeneration time = %0.3f +/- std %0.3f', log(2)/aCG(2), std(log(2)/storeCG(:,2)))
fprintf('\nmax load = %0.3f +/- std %0.3f', exp(aCG(3))*A600, std(exp(storeCG(:,3))*A600))

aCG = min(storeCG);
fprintf('\nmin lag = %0.2f', aCG(1))
fprintf('\nmin generation time = %0.2f', log(2)/aCG(2))
fprintf('\nmin max load = %0.2f', exp(aCG(3))*A600)

aCG = max(storeCG);
fprintf('\nmax lag = %0.2f', aCG(1))
fprintf('\nmax generation time = %0.2f', log(2)/aCG(2))
fprintf('\nmax max load = %0.2f', exp(aCG(3))*A600)
%% Fig 3a: 
subLabels = {"Sulfamethoxazole","Fosfomycin", "Azacitidine", "Pentamidine isethionate"};

fh = figure('color','w','Position', get(0, 'Screensize'), "Visible", 'on'); hold on;
tiledlayout(2,2)

for i= 1:length(subLabels)
    if i == 4
        startColor =[122, 50, 255]./255;
    else
        startColor =colors(i,:);
    end 
    % Find halfAUC using 20.33 hrs: by using fitted Gompertz Model for each
    % concentration and interpolating Gompertz parameters for
    % concentrations to find one corresponding to halfAUC
    % fitted Gompertz Model from "AnalyzeDiffCellConc_Gompertz.m" in "AllGrowthCurves_0.05lag_Set1to35.mat"
    [timeVector,od,conc, aucAtConcentration, controlGompertz] = get_growthcurve_matching_halfauc_GompertzInterp_v6(allCCData,subLabels{i},true, cellConcentration, halfAUC,"MiNoLi wt", timeResolution, timeLimited, gompertz_model, startColor, A600);

    yBlanked = exp(od)*A600;
    % Getting Gompertz model of halfAUC
    [fit_function, lag, growth_rate, max_load] = fit_bacterial_growth_Gompertz_v3(gompertz_model,timeVector, od, false, lower_bound{boundIdx}, upper_bound{boundIdx}, yBlanked', A600);
    

    % Comparing to controlGompertz and making it the no drug control Gompertz parameter if out of bounds
    % controlGompertz = {fit_function, lag, growth_rate, max_load}
    nexttile
    hold on;
    yFitted = exp(fit_function(timeLimited))*A600;
    plot(timeLimited, yFitted, 'Color',startColor, 'LineWidth',2)
    fill([timeLimited fliplr(timeLimited)], [yFitted zeros(1,length(timeLimited))] ,startColor, 'LineWidth',2, 'FaceAlpha',0.5, 'EdgeColor','none');

    title(subLabels{i})
    grid on; box on; 
    ylim([0 0.8])
    xticks([0 5 10 15 20])
    xlim([0 max(timeLimited)])
end
set(gcf,'Position', [1 1 1294 976])

%% Fig 3c: Heatmap Clustering by Gompertz paramters from half AUC

% fosInd = find(strcmp(drugLabels, 'Fosfomycin'));
% idxs = 1:length(drugLabels);
% idxs(fosInd) = [];
% drugLabelsNfos = drugLabels(idxs);
tempData = {lags', gen_times', max_loads'};
data = [];
for i2 = 1:3
    temp = tempData{i2};
    norm = (temp-min(temp))/(max(temp)-min(temp));
    data = [data, norm];
end

colLabels = ["lag", "doubling time", "load"];

% Create the figure
figure('color','w');
% making max_load negative
data(:,3) = -1*data(:,3)+1;
% Impute missing values using knnimpute
imputedData = knnimpute(data);

% Perform hierarchical clustering for rows
rowDist = pdist(imputedData, 'euclidean'); % Pairwise distance for rows
rowLink = linkage(rowDist, 'ward'); % Hierarchical clustering linkage
%rowDist = pdist(imputedData, 'correlation'); % Pairwise distance for rows
%rowLink = linkage(rowDist, 'centroid'); % Hierarchical clustering linkage
subplot(1, 3, 1);
rowOrder1 = optimalleaforder(rowLink, rowDist);
cutoff = 0.75;
clusters= cluster(rowLink, 'cutoff',cutoff, 'criterion', 'distance'); %same order as drugLabels
[H1,T,outperm]= dendrogram(rowLink, 0, 'Orientation', 'left','Reorder', rowOrder1,'Labels', drugLabels, 'Colorthreshold', cutoff); % Get the row order
% Extract colors used in the dendrogram
dendroColors = arrayfun(@(h) get(h, 'Color'), H1, 'UniformOutput', false);

rowOrder = fliplr(outperm);
title('Row Dendrogram');
%axis off;
% Reorder the data matrix
clusteredData = data(rowOrder,:);
% Plot the heatmap
subplot(1, 3, 2:3);
temp = drugLabels(rowOrder);
rowLabels = cell(length(temp),3);
for p1 = 1:length(temp)
    rowLabels{p1,1} = sprintf("%i: %s",p1, temp{p1});
    rowLabels{p1,2} = p1;
    rowLabels{p1,3} = temp{p1};
end
h = heatmap(colLabels, rowLabels(:,2), clusteredData, 'CellLabelColor','none', 'GridVisible','off');
h.Colormap = gray; % Set colormap
title('Heatmap:half AUC');

%% Fig 3b: Triangle no color
%tempDataCC = {[gen_times'; gen_times_CC'], [max_loads'; max_loads_CC'], [lags'; lags_CC']};
tempDataCC = {[gen_times'], [max_loads'], [lags']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    norm = (temp-min(temp))/(max(temp)-min(temp));
    dataCC = [dataCC, norm];
end
% making max_load negative
dataCC(:,2) = -1*dataCC(:,2)+1;

[x, y, A, B, C] = plotDrugEffects(dataCC);

% Plot the triangle
figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
    'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
% Plot the data points
radius = 0.012;
for k = 1:length(rowLabels)
    if k<=length(drugLabels)
        label = rowLabels{strcmp(rowLabels(:,3), drugLabels{k}), 2};
        rectangle('Position', [x(k)-radius, y(k)-radius, 2*radius, 2*radius],'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1);
        text(x(k), y(k), string(label) ,'HorizontalAlignment', 'center', 'FontSize',12)
        %text(x(k), y(k), drugLabelsNfos{k} ,'HorizontalAlignment', 'center')
    end
end
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'r', 'filled');
p(3) = scatter(nan,nan,'kd', 'filled');
p(4) = scatter(nan,nan,'rd', 'filled');
%legend(p, {'wt: no drug degradation','wt: degrades drug at 24hr','resistance cassette: no drug degradation','resistance cassette: degrades drug at 6hr'})
% Annotate the plot
axis equal;
axis off;
clear fh
%save('halfAUCInterpolation_GompertzInterpolation.mat')
%% Fig 3d: Triangle Cluster color
% Plot the triangle
figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
    'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
% Plot the data points
radius = 0.012;
unicmap =[];
unicmap = unique(cell2mat(dendroColors), 'rows');
unicmap = unicmap([2 6 4 5 7 3],:); % reordering to match clusters

for k = 1:length(rowLabels)
    if k<=length(drugLabels)
        label = rowLabels{strcmp(rowLabels(:,3), drugLabels{k}), 2};
        rectangle('Position', [x(k)-radius, y(k)-radius, 2*radius, 2*radius],'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor',unicmap(clusters(k),:));
        text(x(k), y(k), string(label) ,'HorizontalAlignment', 'center')
        %text(x(k), y(k), drugLabelsNfos{k} ,'HorizontalAlignment', 'center')
    end
end
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');

% Annotate the plot
axis equal;
axis off;

figure; hold on;
for j= 1:6
scatter(1,j,'filled', 'MarkerFaceColor',unicmap(j,:))
end
%% Fig 3e: Triangle colored by Mechanism 
% Plot the triangle
figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
    'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
% Plot the data points
radius = 0.012;

for k = 1:length(rowLabels)
    if k<=length(drugLabels)
        idxs = strcmp(rowLabels(:,3), drugLabels{k});
        label = rowLabels{idxs, 2};
        mech = allCCData{strcmp(rowLabels{idxs,3}, allCCData.DrugName), 'Mechanism'}{1};
        fColor = ABXmech.Colors(contains(ABXmech.Mechanism, mech));
        if isempty(fColor)
            fColor ='none';
        else
            fColor = fColor{1};
        end
        rectangle('Position', [x(k)-radius, y(k)-radius, 2*radius, 2*radius],'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor',fColor);
        text(x(k), y(k), string(label) ,'HorizontalAlignment', 'center')
        %text(x(k), y(k), drugLabels{k} ,'HorizontalAlignment', 'center')
    end
end
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'r', 'filled');
p(3) = scatter(nan,nan,'kd', 'filled');
p(4) = scatter(nan,nan,'rd', 'filled');
%legend(p, {'wt: no drug degradation','wt: degrades drug at 24hr','resistance cassette: no drug degradation','resistance cassette: degrades drug at 6hr'})
% Annotate the plot
axis equal;
axis off;

%% Saving data
save("halfAUCInterpolation_LN_2025_03_20.mat")
fprintf('saved!')