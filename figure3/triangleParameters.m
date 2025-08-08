%% Plots 3 growth cuve paramters in a triangle
% Requires lags', gen_times', max_loads' from half AUC interpolation
load('auc_simulation_0.015_errorRate_4.7halfAUC.mat')
load('FunctionalAssayData_0.05lag_2025-02-06.mat')
%% Plotting 3D as triangle with cassette
% Removing NDM1 & fosfomycin 
ccIdx = [1 3:8];
noFosIdx = ~contains(drugLabels, 'Fosfomycin');
drugLabelsNfos = drugLabels(~contains(drugLabels, 'Fosfomycin'));
tempDataCC = {[gen_times(noFosIdx)'; gen_times_CC(ccIdx)'], [max_loads(noFosIdx)'; max_loads_CC(ccIdx)'], [lags(noFosIdx)'; lags_CC(ccIdx)']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    if i2 == 3 % less then no drug lag is reset to no drug lag
        temp(temp<control.lag) = control.lag;
    elseif i2 == 1 % more than no drug GR is reset to no drug GR
        temp(temp<control.doubling_time) = control.doubling_time;
    elseif i2 == 2 % more than no drug load is reset to no drug load
        temp(temp>control.max_load) = control.max_load;
    end
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
allLabels = [drugLabelsNfos; ccLabels(ccIdx,2)];
tInd = [3 9 12 22 25 34 39:45];
% Plot the data points
for k = 1:length(allLabels)
    if k>length(drugLabelsNfos)
        % cassette
        ccInd = strcmp(ccLabels{ccIdx(k-length(drugLabelsNfos)),2},ccFunctionalAssayResults.Strain);
        if max(ccFunctionalAssayResults{ccInd,[4 5]}, [],'all')==1 % 6hr
            scatter(x(k), y(k), 50,'rd', 'filled');
        else
            scatter(x(k), y(k), 50,'kd', 'filled');
        end
    elseif sum(FunctionalAssayResults{contains(FunctionalAssayResults.DrugName,allLabels{k}),[3 4]} == 1)>0
        scatter(x(k), y(k), 50,'r', 'filled');
    else
        scatter(x(k), y(k), 50,'k', 'filled');
    end 
    % if sum(k == tInd)>0
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 0)
    %end
end 
text(A(1), A(2), 'doubling time', 'VerticalAlignment', 'top');
text(B(1), B(2), 'Load', 'VerticalAlignment', 'top');
text(C(1), C(2), 'Lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'r', 'filled');
p(3) = scatter(nan,nan,'kd', 'filled');
p(4) = scatter(nan,nan,'rd', 'filled');
legend(p, {'wt: no drug degradation','wt: degrades drug at 24hr','resistance cassette: no drug degradation','resistance cassette: degrades drug at 6hr'})
% Annotate the plot
axis equal;
axis off;
%sgtitle('Drug Degradation')

%% Plotting 3D as triangle with cassette: colored by functional assay
% Removing NDM1 & fosfomycin 
ccIdx = [1 3:8];
noFosIdx = ~contains(drugLabels, 'Fosfomycin');
drugLabelsNfos = drugLabels(~contains(drugLabels, 'Fosfomycin'));
tempDataCC = {[gen_times(noFosIdx)'; gen_times_CC(ccIdx)'], [max_loads(noFosIdx)'; max_loads_CC(ccIdx)'], [lags(noFosIdx)'; lags_CC(ccIdx)']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    if i2 == 3 % less then no drug lag is reset to no drug lag
        temp(temp<control.lag) = control.lag;
    elseif i2 == 1 % more than no drug GR is reset to no drug GR
        temp(temp<control.doubling_time) = control.doubling_time;
    elseif i2 == 2 % more than no drug load is reset to no drug load
        temp(temp>control.max_load) = control.max_load;
    end
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
allLabels = [drugLabelsNfos; ccLabels(ccIdx,2)];
tInd = [3 9 12 22 25 34 39:45];
% Plot the data points
for k = 1:length(allLabels)
    if k>length(drugLabelsNfos)
        % cassette
        ccInd = strcmp(ccLabels{ccIdx(k-length(drugLabelsNfos)),2},ccFunctionalAssayResults.Strain);
        if max(ccFunctionalAssayResults{ccInd,[4 5]}, [],'all')==1 % 6hr
            scatter(x(k), y(k), 50,'rd', 'filled');
        else
            scatter(x(k), y(k), 50,'kd', 'filled');
        end
    elseif sum(FunctionalAssayResults{contains(FunctionalAssayResults.DrugName,allLabels{k}),[3 4]} == 1)>0
        scatter(x(k), y(k), 50,'r', 'filled');
    else
        scatter(x(k), y(k), 50,'k', 'filled');
    end 
    % if sum(k == tInd)>0
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 0)
    %end
end 
text(A(1), A(2), 'doubling time', 'VerticalAlignment', 'top');
text(B(1), B(2), 'Load', 'VerticalAlignment', 'top');
text(C(1), C(2), 'Lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'r', 'filled');
p(3) = scatter(nan,nan,'kd', 'filled');
p(4) = scatter(nan,nan,'rd', 'filled');
legend(p, {'wt: no drug degradation','wt: degrades drug at 24hr','resistance cassette: no drug degradation','resistance cassette: degrades drug at 6hr'})
% Annotate the plot
axis equal;
axis off;
%sgtitle('Drug Degradation')

%% Plotting 3D as triangle with individual cassette: 
% Removing NDM1
% ccIdx = [1 3:8];
% tempDataCC = {[gen_times'; gen_times_CC(ccIdx)'], [max_loads'; max_loads_CC(ccIdx)'], [lags'; lags_CC(ccIdx)']};
% dataCC = [];
% for i2 = 1:3
%     temp = tempDataCC{i2};
%     if i2 == 3 % less then no drug lag is reset to no drug lag
%         temp(temp<control.lag) = control.lag;
%     elseif i2 == 1 % more than no drug GR is reset to no drug GR
%         temp(temp<control.doubling_time) = control.doubling_time;
%     elseif i2 == 2 % more than no drug load is reset to no drug load
%         temp(temp>control.max_load) = control.max_load;
%     end
%     norm = (temp-min(temp))/(max(temp)-min(temp));
%     dataCC = [dataCC, norm];
% end
% % making max_load negative
% dataCC(:,2) = -1*dataCC(:,2)+1;

[x, y, A, B, C] = plotDrugEffects(dataCC);
 
% Plot the triangle
figure;
lineInds = {[9 39]; [24 43]; [21 42];[33 44];[3 38]; [9 40]; [12 41]};
for o2 =1:7
subplot(3,3,o2)
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;

allLabels = [drugLabelsNfos; ccLabels(ccIdx,2)];
wtL = length(drugLabelsNfos);
unicmap = unique(cell2mat(dendroColors), 'rows');
%unicmap = unicmap([1 3 2],:);
unicmap = unicmap([1 3 5 4 2],:); % reordering to match clusters
% Plot the data points
p=[];
tInd = [3 38 39 9 40 12 41 21 42 24 43 33 44];
for k = lineInds{o2}
    if k>wtL
        %cassette
        p(1) = scatter(x(k), y(k), 50,'b', 'filled');
    else
        p(2)=scatter(x(k), y(k), 50,'k', 'filled');
        %scatter(x(k), y(k), 50,unicmap(clusters(k),:), 'filled');
    end 
    % if sum(k == tInd)>0
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Rotation', 0)
    %end
end 

% for o1 = 1:6
    plot(x(lineInds{o2}), y(lineInds{o2}), 'k');
%end 
text(A(1), A(2), 'generation time(h)', 'VerticalAlignment', 'top');
text(B(1), B(2), 'max load (Abs(OD600))', 'VerticalAlignment', 'top');
text(C(1), C(2), 'lag(h)', 'VerticalAlignment', 'bottom');
% o1=0;
% for o = [1 4 3 5 2]
%     o1 = o1+1;
%     p(o1+1) = scatter(nan, nan, 50,unicmap(o,:), 'filled');
% end 
legend(p, {'resistance cassette','wt'}, 'Location', 'northeast')
%legend(p, {'resistance cassette','wt: cluster 1', 'wt: cluster 2', 'wt: cluster 3', 'wt: cluster 4', 'wt: cluster 5'})
% Annotate the plot
axis equal;
axis off;
%sgtitle('Cluster')
end 
%% Plotting 3D as triangle with cassette: colored by cluster  
% Removing NDM1
% ccIdx = [1 3:8];
% tempDataCC = {[gen_times'; gen_times_CC(ccIdx)'], [max_loads'; max_loads_CC(ccIdx)'], [lags'; lags_CC(ccIdx)']};
% dataCC = [];
% for i2 = 1:3
%     temp = tempDataCC{i2};
%     if i2 == 3 % less then no drug lag is reset to no drug lag
%         temp(temp<control.lag) = control.lag;
%     elseif i2 == 1 % more than no drug GR is reset to no drug GR
%         temp(temp<control.doubling_time) = control.doubling_time;
%     elseif i2 == 2 % more than no drug load is reset to no drug load
%         temp(temp>control.max_load) = control.max_load;
%     end
%     norm = (temp-min(temp))/(max(temp)-min(temp));
%     dataCC = [dataCC, norm];
% end
% % making max_load negative
% dataCC(:,2) = -1*dataCC(:,2)+1;

[x, y, A, B, C] = plotDrugEffects(dataCC);
 
% Plot the triangle
figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;

allLabels = [drugLabelsNfos; ccLabels(ccIdx,2)];
wtL = length(drugLabelsNfos);
unicmap = unique(cell2mat(dendroColors), 'rows');
unicmap = unicmap([1 3 5 4 2],:); % reordering to match clusters
% Plot the data points
p=[];
tInd = [3 9 12 22 25 34 39:45];
for k = 1:length(allLabels)
    if k>wtL
        % cassette
        p(1) = scatter(x(k), y(k), 50,'kd', 'filled');
    else
        scatter(x(k), y(k), 50,unicmap(clusters(k),:), 'filled');
    end 
    %if sum(k == tInd)>0
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 0)
    %end
end 
% lineInds = {[3 39];[40 9 41]; [12 42]; [22 43]; [25 44]; [34 45]};
% for o1 = 1:6
%     plot(x(lineInds{o1}), y(lineInds{o1}), 'k');
% end 
text(A(1), A(2), 'doubling time', 'VerticalAlignment', 'top');
text(B(1), B(2), 'Load', 'VerticalAlignment', 'top');
text(C(1), C(2), 'Lag', 'VerticalAlignment', 'bottom');

for o = [1 4 3 5 2]
    p(o+1) = scatter(nan, nan, 50,unicmap(o,:), 'filled');
end 
legend(p, {'resistance cassette','wt: cluster 1', 'wt: cluster 2', 'wt: cluster 3', 'wt: cluster 4', 'wt: cluster 5'})
% Annotate the plot
axis equal;
axis off;
sgtitle('Cluster')
%% Plotting 3D as triangle with cassette: colored by slope of IC50 
% Removing NDM1
ccIdx = [1 3:8];
tempDataCC = {[gen_times'; gen_times_CC(ccIdx)'], [max_loads'; max_loads_CC(ccIdx)'], [lags'; lags_CC(ccIdx)']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    if i2 == 3 % less then no drug lag is reset to no drug lag
        temp(temp<control.lag) = control.lag;
    elseif i2 == 1 % more than no drug GR is reset to no drug GR
        temp(temp<control.doubling_time) = control.doubling_time;
    elseif i2 == 2 % more than no drug load is reset to no drug load
        temp(temp>control.max_load) = control.max_load;
    end
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

% check if ICtable drugNames is same order as drugLabels
if sum(~strcmp(ICtable.DrugName, drugLabels))==0
    fprintf('ICtable drugNames is same order as drugLabel\n')
end
allLabels = [drugLabels; ccLabels(ccIdx,3)];% drugLabels same order as wt in dataCC
wtL = length(drugLabels);
cmap = ICtable.HillCoefs(:,2);

% Plot the data points
p=[];
for k = 1:length(allLabels)
    if k>wtL
        % cassette
        scatter(x(k), y(k), 50,'kd', 'filled');
    else
        scatter(x(k), y(k), 50,cmap(k), 'filled', 'MarkerEdgeColor', 'k');
    end 
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 0)
end 
colormap pink
c=colorbar;
c.Label.String = 'Hill Coefficient (slope) of IC50(1:200)';
c.Label.FontSize = 12;
text(A(1), A(2), 'doubling time', 'VerticalAlignment', 'top');
text(B(1), B(2), 'Load', 'VerticalAlignment', 'top');
text(C(1), C(2), 'Lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'kd', 'filled');
legend(p, {'wt','resistance cassette'})
sgtitle('Slope IC50(1:200)')
% Annotate the plot
axis equal;
axis off;
%% Plotting 3D as triangle with cassette: colored by toxicity of compound log10(IC50)
% Removing NDM1
ccIdx = [1 3:8];
tempDataCC = {[gen_times'; gen_times_CC(ccIdx)'], [max_loads'; max_loads_CC(ccIdx)'], [lags'; lags_CC(ccIdx)']};
dataCC = [];
for i2 = 1:3
    temp = tempDataCC{i2};
    if i2 == 3 % less then no drug lag is reset to no drug lag
        temp(temp<control.lag) = control.lag;
    elseif i2 == 1 % more than no drug GR is reset to no drug GR
        temp(temp<control.doubling_time) = control.doubling_time;
    elseif i2 == 2 % more than no drug load is reset to no drug load
        temp(temp>control.max_load) = control.max_load;
    end
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

% check if ICtable drugNames is same order as drugLabels
if sum(~strcmp(ICtable.DrugName, drugLabels))==0
    fprintf('ICtable drugNames is same order as drugLabel\n')
end
allLabels = [drugLabels; ccLabels(ccIdx,3)];% drugLabels same order as wt in dataCC
wtL = length(drugLabels);
cmap = log10(ICtable.IC_1to200);

% Plot the data points
p=[];
for k = 1:length(allLabels)
    if k>wtL
        % cassette
        p(1) = scatter(x(k), y(k), 50,'kd', 'filled');
    else
        p(3) = scatter(x(k), y(k), 50,cmap(k), 'filled','MarkerEdgeColor','k');
    end 
    text(x(k), y(k), allLabels{k}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Rotation', 0)
end 
colormap pink
c=colorbar;
c.Label.String = 'log10(IC50 of 1:200)';
c.Label.FontSize = 12;
text(A(1), A(2), 'doubling time', 'VerticalAlignment', 'top');
text(B(1), B(2), 'Load', 'VerticalAlignment', 'top');
text(C(1), C(2), 'Lag', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'kd', 'filled');
legend(p, {'wt','resistance cassette'})
sgtitle('log10(IC50(1:200))')
% Annotate the plot
axis equal;
axis off;