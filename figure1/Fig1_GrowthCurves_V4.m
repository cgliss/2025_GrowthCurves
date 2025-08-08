%% Figure 1: Growth curves of 1:200
%% User inputs
% Getting user inputs common across all scripts
cd ../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("AllGrowthCurves_ln_20250320.mat")
load("ABXmechColorMap.mat")
cd figure1

%% 1b. Plot Normalized Growth Curves 
% smooth the line, do all the drugs 
x1=cellfun(@length, allCCData.NormalizedData);
maxTime = min(x1); % getting maximum timepoint for all strains
drugLabels = unique(allCCData{:, "DrugName"});
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
[sortedMechs,I] = sort(myMechs);
uniMechs = unique(sortedMechs);
drugLabels = drugLabels(I);
myType = myType(I);
cellConc = 200;
fig = figure('color','w','Position', get(0, 'Screensize'), "Visible", 'on');
tiledlayout(6,7)
i4 =1;
x1=cellfun(@length, allCCData.NormalizedData);
maxTime = min(x1); % getting maximum timepoint for all strains
maxTimeCC = [40 9; 200 10; 1000 12];
ind= []; ind1=1; rows = 6; cols = 8;

for i = 1:length(drugLabels)
    if strcmp(myType{i}, 'nonABX')
        startColor =[0 0 0];
    else
        startColor = ABXmech.Colors{strcmp(ABXmech.Mechanism,sortedMechs{i})};
    end 
    
    RI= find(contains(allCCData.DrugName,drugLabels{i}));
    if ~isempty(allCCData{RI(1),"LogPhaseFileName"}{1}) % skipping no data
         uniDil = unique(allCCData{RI, "uM"});
        nColors = length(uniDil)+1;
        % Generate the colormap
        cmap= [linspace(startColor(1), 1, nColors)', ...
            linspace(startColor(2), 1, nColors)', ...
            linspace(startColor(3), 1, nColors)'];
        ind = [ind ,i];
       
        colMap = [uniDil, cmap(1:(length(uniDil)),:)];
        uniStrain = "MiNoLi wt";
        for i2 = 1:height(uniStrain)
            nexttile
            hold on;
            i4 = i4 + 1;
            rowInds=RI(contains(allCCData{RI,"Strain"}, uniStrain(i2)) & allCCData{RI,"CellDilutionFactor"}==cellConc);
            curConc = allCCData{rowInds, "uM"};
            [ScurConc,inxSort] = sort(curConc, 'descend');
            inx = rowInds(inxSort);
            numNoCurve = 0; % number of curves that are not fit based on criteria
            lagParameter  = 0.05; % used to calculate time to reach lagParameter* maxOD
            for i1 = 1:length(inx)
                iRow = inx(i1);
                dilRow = allCCData{iRow,"uM"};
                y = allCCData{iRow,"BlankedMedianData"}{1,1};
                timeInterval = 20/60;
                % plotting data
                x = 0:timeInterval:((length(y)-1)*timeInterval);
                yNorm = allCCData{iRow,"BlankedMedianData"}{1};
                x(isnan(yNorm)) = [];
                yNorm(isnan(yNorm)) = [];
                for i8=1:width(yNorm)
                    ydata = smoothdata(yNorm(:,i8), 'movmean', 3); % smoothing by getting the mean of each hour
                    if dilRow == 0
                    plot(x,ydata,'-','Color',[176, 157, 136]/255, 'LineWidth',2)
                    else
                        plot(x,ydata,'-','Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)
                    end 
                end
            end
            title(drugLabels{i})
            grid on; box on; %axis square;
            xlim([0,lastTimePoint*timeResolution])
            ylim([-0.02,0.8])

        end
    end
end

% Assuming fig is the handle to the figure you are trying to save...
%set the figure's Renderer property to 'painters' before saving as EPS/SVG:
 set(fig, 'Renderer', 'painters');
% % % If saving as an EPS file, in R2020a and later, you can use the exportgraphics command and specify the 'ContentType' as 'vector':
 exportgraphics(fig, 'figure1b_medianData_NoDrugDotted.eps', 'ContentType', 'vector');

%% Figure 1c: Piechart of mechanisms
cats = categorical(sortedMechs);
counts = countcats(cats);
labels = categories(cats);
labels = labels([1 2 5 3 6 7 4 8]); % reordering 
counts = counts([1 2 5 3 6 7 4 8]);
% Plot the pie chart
figure;
h = pie(counts);

% Add count labels
for i = 1:2:length(h) % h contains text objects at odd indices
    h(i+1).String = sprintf("%s %i", labels{(i+1)/2},counts((i+1)/2)); 
end


%% Figure 1a
x1=cellfun(@length, allCCData.NormalizedData);
maxTime = min(x1); % getting maximum timepoint for all strains
drugLabels = unique(allCCData{:, "DrugName"});
myMechs =[];
mycolormap = zeros(length(allCCData.Type), 3);
mycolormap(:,2) = 0.5;
mycolormap(contains(allCCData.Type,"nonABX"),1) = 1;
mycolormap(contains(allCCData.Type,"Cassette"),3) = 1;
for dl1 = 1:length(drugLabels)
    RI= find(contains(allCCData.DrugName,drugLabels{dl1}));
    myMechs = [myMechs, allCCData{RI(1), "Type"}];
end
[~,I] = sort(myMechs);
drugLabels = drugLabels(I);
cellConc = 200;
fig = figure('color','w','Position', get(0, 'Screensize'), "Visible", 'on');
%tiledlayout(4,10)
i4 =1;
x1=cellfun(@length, allCCData.NormalizedData);
maxTime = min(x1); % getting maximum timepoint for all strains
maxTimeCC = [40 9; 200 10; 1000 12];
ind= []; ind1=1; rows = 6; cols = 8;

for i = 10%1:length(drugLabels)
    RI= find(contains(allCCData.DrugName,drugLabels{i}));
    if ~isempty(allCCData{RI(1),"LogPhaseFileName"}{1}) % skipping no data
        ind = [ind ,i];
        uniDil = unique(allCCData{RI, "uM"});
        
        uniStrain = unique(allCCData{RI, "Strain"});
        nColors = length(uniDil)+1;
        startColor= [96, 57, 19]./255;
        % Generate the colormap
        colMap= [linspace(startColor(1), 1, nColors)', ...
            linspace(startColor(2), 1, nColors)', ...
            linspace(startColor(3), 1, nColors)'];
        colMap = [uniDil, colMap(1:length(uniDil),:)]; 
        for i2 = 1:height(uniStrain)
                nexttile 
                hold on;
                i4 = i4 + 1;
                rowInds=RI(contains(allCCData{RI,"Strain"}, uniStrain(i2)) & allCCData{RI,"CellDilutionFactor"}==cellConc);
                curConc = allCCData{rowInds, "uM"};
                [ScurConc,inxSort] = sort(curConc);
                inx = rowInds(inxSort);
                numNoCurve = 0; % number of curves that are not fit based on criteria
                lagParameter  = 0.05; % used to calculate time to reach lagParameter* maxOD
                subplot(1,2,1); hold on;
                for i1 = 1:4
                    iRow = inx(i1);
                    dilRow = allCCData{iRow,"uM"};
                    yAll = allCCData{iRow,"BlankedData"}{1,1};
                    yStdev = std(yAll,[],2)/sqrt(3);
                    y = allCCData{iRow,"BlankedMedianData"}{1,1};
                    timeInterval = 20/60;
                    x = 0:timeInterval:((length(y)-1)*timeInterval);
                    % plotting median data
                    x(isnan(y)) = [];
                    y(isnan(y)) = [];
                    for i8=1:width(yAll)
                        plot(x, yAll(:,i8),'Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)
                    end
                end
                %title(drugLabels{i})
                xlabel('time')
                ylabel('Abs(OD600)')
                box on; axis square;
                xlabel('time')
                ylabel('Abs(OD600)')
                box on; axis square; grid on;
                xlim([0,21])
                ylim([-0.02,0.8])

                subplot(1,2,2); hold on;
                for i1 = 1:length(inx)-2
                    iRow = inx(i1);
                    dilRow = allCCData{iRow,"uM"};
                    timeInterval = 20/60;
                    yAll = allCCData{iRow,"BlankedData"}{1,1};
                    yStdev = std(yAll,[],2)/sqrt(3);
                    yB = allCCData{iRow,"BlankedMedianData"}{1,1};
                    y = allCCData{iRow,"MedianData"}{1,1};
                    x = 0:timeInterval:((length(y)-1)*timeInterval);
                    x(isnan(y)) = [];
                    y(isnan(y)) = [];
                    fit_function = allCCData{iRow,"fit_function"}{1,1};
                    if ~isempty(fit_function)
                        plot(x,exp(fit_function(x))*A600,'Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)
                    end 
                end
                %title(drugLabels{i})
                xlabel('time')
                ylabel('Abs(OD600)')
                box on; axis square; grid on;
                xlim([0,20])
                ylim([-0.02,0.8])
        end
    end
end

