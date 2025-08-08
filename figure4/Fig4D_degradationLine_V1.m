%% To determine if there is the line along y axis that can best predict inactivation vs no
% Run Fig4_Cassettes_V2.m first
%% Determine false positive/ negative based on line location 
% funcTriangle = [x, y, inactivate[0=no 1=yes]
% True Positives (TP): Correctly classified positive samples.
% False Positives (FP): Incorrectly classified negative samples as positive.
% False Negatives (FN): Incorrectly classified positive samples as negative.
% True Negatives (TN): Correctly classified negative samples.

figure;
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
%yline(lineY(max(accuracy)==accuracy))
axis equal;

%%
lineY = linspace(A(2), C(2),1000);
posNeg = nan(length(lineY),4);
accuracy = nan(length(lineY),1);
specificity = nan(length(lineY),1);
FPR = nan(length(lineY),1);
FNR = nan(length(lineY),1);
PR = nan(length(lineY),1);
TPR = nan(length(lineY),1);
euclideanDistance = nan(length(lineY),1);
YoudenJStats = nan(length(lineY),1);
sensitivity = nan(length(lineY),1);
for j = 1:length(lineY)
    % Divide into groups [yes no] invactivation
    idx=funcTriangle(:,2)>=lineY(j);
    groupKnown = funcTriangle(:,3);
    groupPredict = zeros(length(funcTriangle),1)+idx; % assume above line is inactivate, below line is degrade 
    CM= confusionmat(groupKnown, groupPredict);
    TP = CM(2,2);
    FP = CM(1,2);
    FN = CM(2,1);
    TN = CM(1,1);
    posNeg(j,:)= [TP FP FN TN];

    accuracy(j) = (TP + TN) / sum(CM(:));
    sensitivity(j) = TP/(TP+FN);
    specificity(j) = TN / (TN + FP); % true negative rate 
    PR(j) = TP/(TP+FP);
    FPR(j) = FP / (FP + TN);
    TPR(j) = TP/(TP + FN);
    FNR(j) = FN / (FN + TP);
    % find best cutoff line using eucliden distance to [0,1]
    euclideanDistance(j) = pdist([FPR(j), TPR(j); 0,1], 'euclidean');

    % find Youden's J Statistic 
    YoudenJStats(j) =  specificity(j) + sensitivity(j)-1; 
end 

figure 
tiledlayout(4,1)
nexttile
plot(lineY,posNeg(:,1))
title('true pos')
nexttile
plot(lineY,posNeg(:,2))
title('false pos')
nexttile
plot(lineY,posNeg(:,3))
title('false neg')
nexttile
plot(lineY,posNeg(:,4))
title('true neg')

figure
tiledlayout(5,1)
nexttile 
hold on;
plot(lineY, accuracy, '-', 'LineWidth',2)
plot(lineY(max(accuracy)==accuracy), accuracy(max(accuracy)==accuracy), 'r.', 'LineWidth',2)
grid on;
ylim([0 1])
xlim([0, C(2)])
title('accuracy')
nexttile 
plot(lineY, specificity, '-', 'LineWidth',2)
grid on;xlim([0, C(2)])
title('specificity')
nexttile 
plot(lineY, FPR, '-', 'LineWidth',2)
grid on;xlim([0, C(2)])
title('false positive rate')
nexttile 
plot(lineY, FNR, '-', 'LineWidth',2)
grid on;xlim([0, C(2)])
title('false negative rate')
ylabel('line y')
nexttile 
plot(lineY, PR, '-', 'LineWidth',2)
grid on;xlim([0, C(2)])
title('precision')
ylabel('line y')

%% Plotting ROC curve 
figure; 
subplot(1,2,1)
hold on;
plot(FPR, TPR)
plot(0:0.1:1, 0:0.1:1,'k--')
AUC = trapz(flipud(FPR),flipud(TPR));
title(sprintf('AUC = %.2f', AUC))
maxIdx = find(euclideanDistance==min(euclideanDistance));
scatter(FPR(maxIdx(1)), TPR(maxIdx(1)), 'red','filled')
text(FPR(maxIdx(1)), TPR(maxIdx(1)), sprintf('optimal cutoff (minimum distance = %.2f; Youden J Stats = %.2f)',euclideanDistance(maxIdx(1)), YoudenJStats(maxIdx(1))))
maxIdxJ = find(YoudenJStats==max(YoudenJStats));
scatter(FPR(maxIdxJ), TPR(maxIdxJ), 'green','filled')
text(FPR(maxIdxJ(1)), TPR(maxIdxJ(1)), sprintf('optimal cutoff (minimum distance = %.2f; Youden J Stats = %.2f)\n sensitivity=%.2f specificity=%.2f',euclideanDistance(maxIdxJ(1)), YoudenJStats(maxIdxJ(1)), sensitivity(maxIdxJ(1)), specificity(maxIdxJ(1))))

grid on; box on; axis square;
ylabel('true positive rate (sensitivity)')
xlabel('false positive rate(1-specificity: false alarm)')

subplot(1,2,2)
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
yline(lineY(maxIdx(1)))
yline(lineY(maxIdxJ(1)))
axis equal;
axis off;

%% Bar plot of black vs white 
line10 = linspace(A(2), C(2),21);

degrade =nan(length(line10),1); 
noDegrade = nan(length(line10),1);  
for k1 = 1:length(line10)
    idx=funcTriangle(funcTriangle(:,2)>=line10(k1), 3);
    degrade(k1) = sum(idx)/length(idx)*100;
    noDegrade(k1) = sum(~idx)/length(idx)*100;
end 
figure
tiledlayout(1,2)
nexttile
bh=barh([degrade, noDegrade], 1,'stacked');
set(bh, 'FaceColor', 'Flat')
xticks([0 25 50 75 100])
xlim([0 100])
bh(1).CData = [1 1 1];  % Change color to first level
bh(2).CData = [0 0 0];  % Change color to second level, etc...


nexttile
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
for k2 = line10(1:end)
    yline(k2)
end 
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
%yline(lineY(max(accuracy)==accuracy))
axis equal;
axis off; 
%% Ratio of black to white in total triangle
figure
black = sum(funcTriangle(:,3)==0)/length(funcTriangle);
white = sum(funcTriangle(:,3)==1)/length(funcTriangle);
bh = bar(log2(black/white));
%set(bh, 'FaceColor', 'Flat')
%bh(1).CData = [1 1 1];  % Change color to first level
%bh(2).CData = [0 0 0];  % Change color to second level, etc...
%ylim([0.9 1])
grid on;
xlabel('ratio')

%% Bar plot of black vs white 
line10 = linspace(A(2), C(2),21);

degrade =nan(length(line10),1); 
noDegrade = nan(length(line10),1);  
for k1 = 1:length(line10)
    idx=funcTriangle(funcTriangle(:,2)<=line10(k1), 3);
    degrade(k1) = sum(idx)/length(idx)*100;
    noDegrade(k1) = sum(~idx)/length(idx)*100;
end 
figure
tiledlayout(1,2)
nexttile
bh=barh(log2(degrade./noDegrade));
%set(bh, 'FaceColor', 'Flat')
%xticks([0 25 50 75 100])
%xlim([0 100])
%bh(1).CData = [1 1 1];  % Change color to first level
%bh(2).CData = [0 0 0];  % Change color to second level, etc...


nexttile
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
for k2 = line10(1:end)
    yline(k2)
end 
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
%yline(lineY(max(accuracy)==accuracy))
axis equal;
axis off; 
%% Bar plot of black vs white 
line10 = linspace(A(2), C(2),21);

degrade =nan(length(line10),1); 
noDegrade = nan(length(line10),1);  
for k1 = 1:length(line10)
    idx=funcTriangle(funcTriangle(:,2)<=line10(k1), 3);
    degrade(k1) = sum(idx);
    noDegrade(k1) = sum(idx ==0);
end 
figure
tiledlayout(1,2)
nexttile; hold on
plot(degrade/sum(funcTriangle(:,3)==1),'k');
plot(noDegrade/sum(funcTriangle(:,3)==0),'b');
grid on;
%set(bh, 'FaceColor', 'Flat')
%xticks([0 25 50 75 100])
%xlim([0 100])
%bh(1).CData = [1 1 1];  % Change color to first level
%bh(2).CData = [0 0 0];  % Change color to second level, etc...


nexttile
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
for k2 = line10(1:end)
    yline(k2)
end 
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
%yline(lineY(max(accuracy)==accuracy))
axis equal;
axis off; 
%% ratio / ratio black vs white 
line10 = linspace(A(2), C(2),21);

degrade =nan(length(line10),1); 
noDegrade = nan(length(line10),1);  
ratioT = nan(length(line10),1); 
ratioU = nan(length(line10),1); 
ratioD = nan(length(line10),1); 
for k1 = 1:length(line10)
    uidx=funcTriangle(funcTriangle(:,2)>=line10(k1), 3);
    didx=funcTriangle(funcTriangle(:,2)<line10(k1), 3);
    ratioU(k1) = sum(uidx==0)/sum(uidx);
    ratioD(k1) = sum(didx==0)/sum(didx);
    ratioT(k1) = ratioU(k1)/ratioD(k1);
end 
figure
tiledlayout(1,3)
nexttile; hold on;
plot(line10,ratioT,'o-k');
grid on;

nexttile; hold on;
plot(line10,ratioU,'o-b');
plot(line10,ratioD,'o-g');
legend({'above', 'below'})
grid on;
%set(bh, 'FaceColor', 'Flat')
%xticks([0 25 50 75 100])
%xlim([0 100])
%bh(1).CData = [1 1 1];  % Change color to first level
%bh(2).CData = [0 0 0];  % Change color to second level, etc...


nexttile
fill([A(1), B(1), C(1)], [A(2), B(2), C(2)], 'w', 'EdgeColor', 'k'); % Triangle outline
hold on;
colors = [255 190 11; 255 0 110; 58 134 255]./255;
patch('Faces', [1 2 3], 'Vertices', [A; B; C], 'FaceVertexCData', colors, ...
      'FaceColor', 'interp', 'EdgeColor', 'none'); % Interpolated color
bw = [0 0 0; 1 1 1];
for i= 1:length(funcTriangle)
scatter(funcTriangle(i,1), funcTriangle(i,2), 'filled', 'MarkerFaceColor',bw(funcTriangle(i,3)+1,:))
end
for k2 = line10(1:end)
    yline(k2)
end 
text(A(1), A(2), 'slow growth', 'VerticalAlignment', 'top');
text(B(1), B(2), 'low capacity', 'VerticalAlignment', 'top');
text(C(1), C(2), 'long lag', 'VerticalAlignment', 'bottom');
%yline(lineY(max(accuracy)==accuracy))
axis equal;
axis off; 
%% Enrichment plot 
load('FuncAssayResults.mat')
y = funcTriangle(:,2);
tf = funcTriangle(:,3);
N = length(y);
% Convert tf to +1/-1
tf_adjusted = double(tf);
tf_adjusted(tf_adjusted == 0) = -1;
% ==== Sort by position ====
[y_sorted, sortIdx_pos] = sort(y);
tf_sorted_pos = tf(sortIdx_pos);
tf_adjusted_sorted = tf_adjusted(sortIdx_pos);
% ==== 1. Moving Average by Position ====
window_half_width = 0.1;  % units in y
movAvg_by_pos = zeros(N,1);
for i = 1:N
    low = y_sorted(i) - window_half_width;
    high = y_sorted(i) + window_half_width;
    in_window = (y_sorted >= low) & (y_sorted <= high);
    movAvg_by_pos(i) = mean(tf_adjusted_sorted(in_window));
end
% ==== 2. Manhattan-like Plot (effect vs. rank) ====
[effect_sorted, sortIdx_effect] = sort(y, 'descend');
tf_sorted_effect = tf(sortIdx_effect);
barcodePos = find(tf_sorted_effect == 1);
% ==== 3. Barcode Plot ====
% (uses barcodePos from Manhattan data)
barcodePos_by_y = y(find(tf_sorted_effect == 1));
% ==== 4. Enrichment Score by Position ====
Nh = sum(tf_sorted_pos);
Nm = N - Nh;
hit_step = 1 / Nh;
miss_step = -1 / Nm;
ES = zeros(N,1);
for i = 2:N
    if tf_sorted_pos(i)
        ES(i) = ES(i-1) + hit_step;
    else
        ES(i) = ES(i-1) + miss_step;
    end
end
% Force ES to return to 0
ES = ES - ES(end);
% ==== PLOT ====
figure;
% 1. Moving Average by Position
subplot(2,2,1);
plot(y_sorted, movAvg_by_pos, 'LineWidth', 1.5);
title('Running Average by Position');
xlabel('Position (y)');
ylabel('Moving Avg (1/-1)');
xlim([min(y_sorted) max(y_sorted)]);
set(gca,'xlim',[0 1]);
% 2. Manhattan-like Plot
subplot(2,2,2);
stem(effect_sorted, '.', 'Marker', 'none');
hold on;
stem(barcodePos, effect_sorted(barcodePos), 'r.');
title('Manhattan-like Plot');
ylabel('y value');
xlim([1 N]);
% 3. Barcode Plot
subplot(2,2,3);
hold on;
for i = 1:length(barcodePos)
    bar(barcodePos(i),1, 'FaceColor', 'k','EdgeColor', 'w','BarWidth',1);
end
title('Barcode Plot');
ylim([0 1]);
xlim([0 N]);
set(gca, 'ytick', []);
% 4. Enrichment Score by Position
subplot(2,2,4);
stairs(y_sorted, ES, 'LineWidth', 1.5);
xlabel('Position (y)');
ylabel('Enrichment Score');
title('Enrichment Score Plot (by Position)');
yline(0, '--k');
set(gca,'xlim',[0 1],'ylim',[-0.65 0.1]);
sgtitle('Skewness Visualizations for Ordered Data');