%% Figure 2: Presenting 3 extreme cases 
%% User inputs
% Getting user inputs common across all scripts
cd ../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("AllGrowthCurves_ln_20250320.mat")
load("ABXmechColorMap.mat")
cd figure2
drugLabels = {"Azacitidine", "Sulfamethoxazole", "Fosfomycin"};
%% Plotting normalized growth curves and growth curve parameters vs AUC
% allCCData.AUC: from timeperiod of 1:maximum time reach capactity for each
% cell concentration
ICpoint = 50;
timeInterval = 20/60; 

gen_time = log(2)./allCCData.growth_rate;
lag = allCCData.lag;
AUC =  allCCData.AUC;
gompertzFits =allCCData.fit_function;
maxLoad = -1*exp(allCCData.max_load)*A600;
fig = figure('color','w');
rows = 3; cols = 3;counter = 1;
plotConfigs = {'lag(h)', [0, 12], lag, 7, 2;...
    'generation time(h)', [0, 2.5], gen_time, 8, 3;...
    '-max load (Abs)', [-1, 0 ], maxLoad, 9, 4};

uniStrain = "MiNoLi wt";
cellConc1 = 200;
legPlots = [];axisS =[];
for i =1:length(drugLabels)
    RI= find(contains(allCCData.DrugName,drugLabels{i}));
    colorsUsed =[];
    for i2 = 1:height(uniStrain)
        subplot(rows, cols, counter); hold on;
        rowInds=RI(contains(allCCData{RI,"Strain"}, uniStrain(i2)) & allCCData{RI,"CellDilutionFactor"}==cellConc1);
        p=[];dils=[];
        for i7=1:length(rowInds)
            % A: Plotting normalized growth curves
            y=log(allCCData{rowInds(i7),"MedianData"}{1});
            x = 0:timeInterval:((length(y)-1)*timeInterval);
            dilRow = allCCData{rowInds(i7),"uM"};
            if dilRow == 0
                continue % skipping no drug
            end 

            % Making color map
            gfIdx = length(rowInds);%~isnan(allCCData{rowInds,"lag"});
            nColors = gfIdx+1; % concentrations that have gompertz fit 
            startColor = colors(i,:);
            % Generate the colormap
            colorSpectrum = [linspace(startColor(1), 1, nColors)', ...
                linspace(startColor(2), 1, nColors)', ...
                linspace(startColor(3), 1, nColors)'];
            cmap = colorSpectrum;
            uniDil = allCCData{rowInds, "uM"};
            colMap = [flipud(uniDil), cmap(1:gfIdx,:)];
            % plot gompertz fit
            gompertzFit = allCCData{rowInds(i7),"fit_function"}{1};
            if ~isempty(gompertzFit) % skipping lines without gompertzfit
                
                %plot normalized growth
                yData=allCCData{rowInds(i7),"NormalizedData"}{1};
                yData = exp(yData(1:length(x),:))*A600;
                yData = smoothdata(yData, 1, 'movmean', 6); % smoothing by getting the mean of 2 hours
                % plot(x, yData)

                stdY = std(yData, [],2);
                yMData = mean(yData, 2);
                % Shade area between y1 and y2
                fill([x, fliplr(x)], [yMData-stdY; flipud(yMData+stdY)]', colMap(colMap(:,1)==dilRow,2:end), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                x = 0:timeInterval:22;
                p=[p, plot(x,exp(gompertzFit(x))*A600,'color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',1)];
                %p=[p, plot(x,exp(gompertzFit(x))*A600, 'LineWidth',1)];
                colorsUsed = [colorsUsed; colMap(colMap(:,1)==dilRow,2:end)];
                dils = [dils, allCCData{rowInds(i7),"uM"}];
            end 
        end
        xlabel('time(h)'); ylabel('normalized absorbance')
        grid on; box on; 
        xlim([0,max(timeLimited)])
        ylim([-0.05,0.8])
        title(drugLabels{i});
        legend(p, string(dils))
        % B: Interpolation 
        subplot(rows,cols,counter+3); hold on;
        % Find halfAUC using 20.33 hrs: by using fitted Gompertz Model for each
        % concentration and interpolating Gompertz parameters for
        % concentrations to find one corresponding to halfAUC
        % fitted Gompertz Model from "AnalyzeDiffCellConc_Gompertz.m" in "AllGrowthCurves_0.05lag_Set1to35.mat"
       [timeVector,od,conc, aucAtConcentration,  controlGompertz, predictedParams] = GC_matching_halfauc_Gompertz_v6_fig2(allCCData,drugLabels{i},true, cellConcentration, halfAUC,"MiNoLi wt", timeResolution, timeLimited, gompertz_model, startColor, A600);

        controlGompertz{4} = -1*exp(controlGompertz{4})*A600; % making max load negative and back to Abs
        controlGompertz{3} = log(2)/controlGompertz{3}; % making growth rate into generation time
        controlFitFunc = controlGompertz{1}; 
        controlAUC = trapz(timeLimited, exp(controlFitFunc(timeLimited))*A600);

        % C:Plotting growth curve parameters vs AUC
        % Getting indexes & concentrations
        curConc = allCCData{rowInds, "uM"};
        % removing no drug
        rowInds(curConc ==0) = [];
        curConc(curConc==0) =[];
        [sortCurConc,inxSort] = sort(curConc);
        inx = rowInds(inxSort);
            
        for plotIdx = 1:height(plotConfigs)
            % Get plot configuration
            yLabel = plotConfigs{plotIdx, 1};
            yLimits = plotConfigs{plotIdx, 2};
            yDataFunc = plotConfigs{plotIdx, 3};
            Range = linspace(0, 12, 100);

            % Get data
            x= NaN(length(inx),1);
            for if1 = 1:length(inx)
                % Using AUC from gompertz model for 1:timeLimited
                fitFunction = gompertzFits{inx(if1)};
                if ~isempty(fitFunction)
                AUCx = trapz(timeLimited,exp(fitFunction(timeLimited))*A600);
                x(if1) = AUCx./controlAUC; 
                end 
            end 
            
            y = yDataFunc(inx);
            nI = isnan(y);
            x(nI) = []; y(nI) =[];
            x = [controlAUC/controlAUC; x];
            y = [controlGompertz{plotConfigs{plotIdx, 5}}; y];
            % sorting 
            [x, sYind] = sort(x);
            y = y(sYind);

            % Plotting Gompertz parameters 
            subplot(rows,cols,plotConfigs{plotIdx, 4})
            hold on;
            
            %plot(x,y,'LineWidth', 2,'Color',cmap(1,:));
            mk  =0;
            for k=find(x~=1)'
                mk = mk+1;
                plot(x(k),y(k),'o','MarkerSize', 10,'LineStyle', 'none','MarkerFaceColor',colorsUsed(mk,:),'MarkerEdgeColor',cmap(1,:));
            end

            % Plot Predicted paramters
            % predictedParams(j1,:) = [aucAtConcentration(j1), fitLags(j1), log(2)/fitGrowth_rates(j1) , fitMax_loads(j1)];
            % aucAtConcentration for limited time
            [~, uIds] = unique(predictedParams(:,1)); % getting only unique parameters
            predictedParams = predictedParams(uIds,:); 
            xP = []; yP = [];
            for k1 = 1:height(predictedParams)
                plot(predictedParams(k1,1),predictedParams(k1,plotConfigs{plotIdx, 5}),'MarkerSize', 6,'LineStyle','none','Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',cmap(1,:));
                xP = [xP, predictedParams(k1,1)];
                yP = [yP, predictedParams(k1,plotConfigs{plotIdx, 5})];
            end 
            xP =[xP, 1];
            yP = [yP, controlGompertz{plotConfigs{plotIdx, 5}}];
            xFit = linspace(min(xP),1, 50);
            p = polyfit(xP,yP,2);
            yFit = polyval(p, xFit);
            plot(xFit,yFit, 'LineWidth', 2,'Color',cmap(1,:))

             % plot control 
            plot(controlAUC/controlAUC,controlGompertz{plotConfigs{plotIdx, 5}}, 'o','MarkerSize', 10,'LineStyle', 'none','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])


            %xlim([0.5 1])
            xlim([0.2 1])
            ylim(yLimits)
            
            ylabel(yLabel);
            xlabel('AUC');
            box on; grid on;
            axisS = [axisS, gca];
        end
        counter = counter+1;
    end
end

set(gcf,'Position', [1 1 1294 976])

% Assuming fig is the handle to the figure you are trying to save...
% % set the figure's Renderer property to 'painters' before saving as EPS/SVG:
set(fig, 'Renderer', 'painters');
% If saving as an EPS file, in R2020a and later, you can use the exportgraphics command and specify the 'ContentType' as 'vector': 
%exportgraphics(fig, 'figure2_v7_ln.eps', 'ContentType', 'vector');