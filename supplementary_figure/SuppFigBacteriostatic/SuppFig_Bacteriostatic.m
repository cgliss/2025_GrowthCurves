%% Supplementary Figure: Bacteriostatic 
%% User inputs
% Getting user inputs common across all scripts
cd ../../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("halfAUCInterpolation_LN_2025_03_20.mat")
cd supplementary_figure/SuppFigBacteriostatic/
bacteriostatic = {'Azithromycin';'Chloramphenicol'; 'Clarithromycin'; 'Sulfamethoxazole'; 'Sulfamonomethoxine'; 'Trimethoprim'; 'Kasugamycin'; 'Tetracycline'};
both = {'Nourseothricin', 'Furazolidone', 'Nitrofurantoin'};
%% Bacteriostatic vs cidal triangle 
% Plot the triangle
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

figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
% Plot the data points
for k = 1:length(drugLabels)
    if sum(contains(allCCData{contains(allCCData.DrugName, drugLabels{k}), 'Type'}, 'nonABX'))==0 
    if sum(contains(bacteriostatic, drugLabels{k}))>0
        scatter(x(k), y(k), 60,'w', 'filled', 'LineWidth',10);
    elseif sum(contains(both, drugLabels{k}))>0
        scatter(x(k), y(k), 60,[128 128 128]/255, 'filled', 'LineWidth',10);
    else
        scatter(x(k), y(k), 60,'k', 'filled', 'LineWidth',10);
    end 
    %text(x(k), y(k),allLabels{k})
    end
end 
text(A(1), A(2), 'generation time(h)', 'VerticalAlignment', 'top');
text(B(1), B(2), 'max load (Abs(OD600))', 'VerticalAlignment', 'top');
text(C(1), C(2), 'lag(h)', 'VerticalAlignment', 'bottom');
p=[];
p(1) = scatter(nan,nan,'k', 'filled');
p(2) = scatter(nan,nan,'w', 'filled');
legend(p, {'wt: bacteriocidal','wt: bacteriostatic'})
% Annotate the plot
axis equal;
axis off;
%sgtitle('Drug Degradation')

