%% Figure 4: Cassettes + Functional assay
%% User inputs
% Getting user inputs common across all scripts
cd ../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("halfAUCInterpolation_LN_2025_03_20.mat")
cd figure4
%% 4b Cassette Triangles 
% Removing NDM1,rmtB, tet(A)   & fosfomycin 
ccIdx = 1:length(gen_times_CC);%[1 3 4 5 7];
tempDataCC = {[gen_times'; gen_times_CC(ccIdx)'], [max_loads'; max_loads_CC(ccIdx)'], [lags'; lags_CC(ccIdx)']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    norm = (temp-min(temp))/(max(temp)-min(temp));
    dataCC = [dataCC, norm];
end
% making max_load negative
dataCC(:,2) = -1*dataCC(:,2)+1;

[x, y, A, B, C] = plotDrugEffects(dataCC);
 
allLabels = [drugLabels; ccLabels(ccIdx,2)];
wtL = length(drugLabels);

% Plot the triangles
figure;
lineInds = cell(length(ccLabels),1);
for i5 = 1:length(ccLabels)
    lineInds{i5} = [find(strcmp(allLabels, ccLabels{i5,1})) find(strcmp(allLabels, ccLabels{i5,2}))];
end 


radius = 0.018;
for o2 = 1:length(lineInds)
    subplot(1,5,o2)
    fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
    colors = [255 190 11; 255 0 110; 58 134 255]./255;
    patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
          'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
    hold on;
    % Plot the data points
    for k = lineInds{o2}
        if k>wtL
            %cassette
            cassXY = [x(k), y(k)];
           scatter(x(k), y(k), 50,'+k');
           text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Rotation', 0)
        else
            wtXY = [x(k) y(k)];
            rectangle('Position', [x(k)-radius, y(k)-radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1);
            label = rowLabels{strcmp(rowLabels(:,3), drugLabels), 2};
            text(x(k), y(k), string(label),'HorizontalAlignment', 'center')
            text(x(k), y(k), allLabels{k},'HorizontalAlignment', 'center')
        end
        
    end
    plot(x(lineInds{o2}), y(lineInds{o2}), 'k');
    % getting angle between cassette and x axis
    dx = cassXY(1) - wtXY(1);
    dy = cassXY(2) - wtXY(2);
    distance = sqrt(dx^2 + dy^2);
    angle = rad2deg(atan2(dy, dx));
    title(sprintf('Angle from x axis: %0.2f\n Magnitude: %0.2f \n YMagnitude/YTriangle = %0.2fpercent', angle, distance, dy/(C(2)-A(2)*100)))
    text(A(1), A(2), 'generation time(h)', 'VerticalAlignment', 'top');
    text(B(1), B(2), 'max load (Abs(OD600))', 'VerticalAlignment', 'top');
    text(C(1), C(2), 'lag(h)', 'VerticalAlignment', 'bottom');
    
    % Annotate the plot
    axis equal;
    axis off;
end
%% Making table of parameter information
paramTable = table(allLabels, [gen_times'; gen_times_CC(ccIdx)'],  [lags'; lags_CC(ccIdx)'], exp([max_loads'; max_loads_CC(ccIdx)'])*A600, 'VariableNames', {'drugNames','generation time(h)',' lag(h)', 'max load(Abs)'});
writetable(paramTable,'ParamTable.xlsx')
%% 4c Functional Assay Triangle
% Plot the triangle
figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
funcTriangle = [x,y,nan(length(x),1)];
% Plot the data points
for k = 1:length(allLabels)
    if k>length(drugLabels)
        % cassette
        ccInd = strcmp(ccLabels{ccIdx(k-length(drugLabels)),2},ccFunctionalAssayResults.Strain);
        if max(ccFunctionalAssayResults{ccInd,[4 5]}, [],'all')==1 % 6hr
            scatter(x(k), y(k), 60,'w+', 'LineWidth',3);
            funcTriangle(k,3) = 1;
        else
            scatter(x(k), y(k), 60,'k+', 'LineWidth',3);
            funcTriangle(k,3) = 0;
        end
    elseif sum(FunctionalAssayResults{contains(FunctionalAssayResults.DrugName,allLabels{k}),[3 4]} == 1)>0
        scatter(x(k), y(k), 60,'w', 'filled', 'LineWidth',10);
        funcTriangle(k,3) = 1;
    else
        scatter(x(k), y(k), 60,'k', 'filled', 'LineWidth',10);
        funcTriangle(k,3) = 0;
    end 
    %text(x(k), y(k),allLabels{k})
end 
text(A(1), A(2), 'generation time(h)', 'VerticalAlignment', 'top');
text(B(1), B(2), 'max load (Abs(OD600))', 'VerticalAlignment', 'top');
text(C(1), C(2), 'lag(h)', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'w', 'filled');
p(3) = scatter(nan,nan,'kd', 'filled');
p(4) = scatter(nan,nan,'wd', 'filled');
legend(p, {'wt: no drug degradation','wt: degrades','resistance cassette: no drug degradation','resistance cassette'})
% Annotate the plot
axis equal;
axis off;
%sgtitle('Drug Degradation')


%% 4 Supplement. Plot Normalized Growth Curves of wt + cassette
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
tiledlayout(3,4)
uniCABX = unique(AllcassetteAllData.DrugName);
for i = 1:length(uniCABX)
    idx= find(strcmp(drugLabels,uniCABX(i)));
    startColor = ABXmech.Colors{strcmp(ABXmech.Mechanism,sortedMechs{idx})};

    RI= find(contains(allCCData.DrugName,drugLabels{idx}));
    uniDil = unique(allCCData{RI, "uM"});
    
    RIc= find(contains(AllcassetteAllData.DrugName,drugLabels{idx}));
    uniDilC = unique(AllcassetteAllData{RIc, "uM"});
    allDil = [uniDil; uniDilC(uniDilC~=0)];
    allDil = sort(allDil);
    % Generate the colormap
    colMap = [allDil, turbo(length(allDil))];
    
    % Plotting wt
    uniStrain = "MiNoLi wt";
    for i2 = 1:height(uniStrain)
        nexttile
        hold on;
        rowInds=RI(contains(allCCData{RI,"Strain"}, uniStrain(i2)) & allCCData{RI,"CellDilutionFactor"}==cellConc);
        curConc = allCCData{rowInds, "uM"};
        [ScurConc,inxSort] = sort(curConc);
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
                plot(x,ydata,'-','Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)
            end
        end
        title(drugLabels{idx})
        grid on; box on; %axis square;
        xlim([0,lastTimePoint*timeResolution])
        ylim([-0.02,1])
    end

    % Plotting cassette 
    uniStrain = unique(AllcassetteAllData{RIc, "Strain"});
    for i2 = 1:height(uniStrain)
        nexttile
        hold on;
        rowInds=RIc(contains(AllcassetteAllData{RIc,"Strain"}, uniStrain(i2)) & AllcassetteAllData{RIc,"CellDilutionFactor"}==cellConc);
        curConc = AllcassetteAllData{rowInds, "uM"};
        [ScurConc,inxSort] = sort(curConc);
        inx = rowInds(inxSort);
        numNoCurve = 0; % number of curves that are not fit based on criteria
        lagParameter  = 0.05; % used to calculate time to reach lagParameter* maxOD
        for i1 = 1:length(inx)
            iRow = inx(i1);
            dilRow = AllcassetteAllData{iRow,"uM"};
            y = AllcassetteAllData{iRow,"BlankedMedianData"}{1,1};
            timeInterval = 20/60;
            % plotting data
            x = 0:timeInterval:((length(y)-1)*timeInterval);
            yNorm = AllcassetteAllData{iRow,"BlankedMedianData"}{1};
            x(isnan(yNorm)) = [];
            yNorm(isnan(yNorm)) = [];
            for i8=1:width(yNorm)
                ydata = smoothdata(yNorm(:,i8), 'movmean', 3); % smoothing by getting the mean of each hour
                plot(x,ydata,'-','Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)
            end
        end
        title(uniStrain{i2})
        grid on; box on; %axis square;
        xlim([0,lastTimePoint*timeResolution])
        ylim([-0.02,1])
    end

end
colormap turbo
