%% Supplementary Figure: Fig2B/C for all drugs 
%% User inputs
% Getting user inputs common across all scripts
cd ../../
[timeResolution,lastPossibleTimePoint, lastTimePoint, cellConcentration, halfAUC,upper_bound, lower_bound, time, timeLimited, gompertz_model, colors, A600, boundOrder]= userInputs;
load("AllGrowthCurves_ln_20250320.mat")
load("ABXmechColorMap.mat")
cd supplementary_figure/SuppFig2BC
%% Plotting normalized growth curves and growth curve parameters vs AUC
% allCCData.AUC: from timeperiod of 1:maximum time reach capactity for each
% cell concentration
drugLabels = unique(allCCData.DrugName);
ICpoint = 50;
timeInterval = 20/60; 

gen_time = log(2)./allCCData.growth_rate;
lag = allCCData.lag;
AUC =  allCCData.AUC;
gompertzFits =allCCData.fit_function;
maxLoad = -1*exp(allCCData.max_load)*A600;
fig = figure('color','w','Position', [1 1 1384 850]);
rows = 5; cols = 5;counter = 1;
plotConfigs = {'lag(h)', [0, 12], lag, 3, 2;...
    'generation time(h)', [0, 2.5], gen_time, 4, 3;...
    '-max load (Abs)', [-1, 0 ], maxLoad, 5, 4};

uniStrain = "MiNoLi wt";
cellConc1 = 200;
legPlots = [];axisS =[];figCount = 1;cot =1;
for i =[22 23 25 28 34 3 14 9 12 17 27 30 4 6 8 20 32 33 37 21 26 1 2 5 7 10 11 13 24 29 31 15 16 18 19 35 36 38]
    if rem(counter,25)==1 && counter ~=1
        saveas(fig,['SuppFig2_' num2str(figCount) '.jpg'])
        % % % set the figure's Renderer property to 'painters' before saving as EPS/SVG:
        set(fig, 'Renderer', 'painters');
        % % If saving as an EPS file, in R2020a and later, you can use the exportgraphics command and specify the 'ContentType' as 'vector': 
        exportgraphics(fig, ['SuppFig2_' num2str(figCount) '.eps'], 'ContentType', 'vector');
        fig = figure('color','w','Position',  [1 1 1384 850]);
        rows = 5; cols = 5;counter = 1; figCount = figCount+1; 
    end 
    RI= find(contains(allCCData.DrugName,drugLabels{i}));
    colorsUsed =[];
    for i2 = 1:height(uniStrain)
        subplot(rows, cols, counter); hold on;
        rowInds=RI(contains(allCCData{RI,"Strain"}, uniStrain(i2)) & allCCData{RI,"CellDilutionFactor"}==cellConc1);
        
        % Making color map
        gfIdx = length(rowInds);%~isnan(allCCData{rowInds,"lag"});
        nColors = gfIdx+1; % concentrations that have gompertz fit
        myType = allCCData{RI(1), "Type"};
        typeA = allCCData{RI(1), "Type"};
        if strcmp(myType, 'nonABX')
            startColor =[0 0 0];
        else
            mech = allCCData{RI(1), "Mechanism"};
            startColor = ABXmech.Colors{strcmp(ABXmech.Mechanism,mech)};
        end
        % Generate the colormap
        colorSpectrum = [linspace(startColor(1), 1, nColors)', ...
            linspace(startColor(2), 1, nColors)', ...
            linspace(startColor(3), 1, nColors)'];
        cmap = colorSpectrum;
        uniDil = sort(allCCData{rowInds, "uM"}, 'descend');
        colMap = [flipud(uniDil), cmap(1:gfIdx,:)];
        
        p=[];dils=[];
        for i7=1:length(rowInds)
            % A: Plotting normalized growth curves
            y=log(allCCData{rowInds(i7),"MedianData"}{1});
            x = 0:timeInterval:((length(y)-1)*timeInterval);
            dilRow = allCCData{rowInds(i7),"uM"};

            
            %plot normalized growth
            yData=allCCData{rowInds(i7),"BlankedMedianData"}{1};
            ydata = smoothdata(yData, 'movmean', 3); % smoothing by getting the mean of each hour
             x = 0:timeInterval:((length(yData)-1)*timeInterval);
             x(isnan(ydata)) =[];
            ydata(isnan(ydata)) =[];
           
            if dilRow == 0
                p=[p,plot(x,ydata,'-','Color',[158, 140, 120]/255, 'LineWidth',2)];
            else
                p=[p,plot(x,ydata,'-','Color',colMap(colMap(:,1)==dilRow,2:end), 'LineWidth',2)];
            end
            colorsUsed = [colorsUsed; colMap(colMap(:,1)==dilRow,2:end)];
            dils = [dils, round(allCCData{rowInds(i7),"uM"},2)];
        end
        dilFac = max(uniDil)/ uniDil(2);
        text(0.25, 0.7, sprintf('max conc. = %g uM\ndil. factor = %g', max(uniDil), dilFac ))
        xlabel('time(h)'); ylabel('normalized absorbance')
        grid on; box on; 
        xlim([0,max(timeLimited)])
        ylim([-0.05,0.8])
        title(drugLabels{i});
        counter = counter+1;

        % B: Interpolation 
        subplot(rows,cols,counter); hold on;
        % Find halfAUC using 20.33 hrs: by using fitted Gompertz Model for each
        % concentration and interpolating Gompertz parameters for
        % concentrations to find one corresponding to halfAUC
        % fitted Gompertz Model from "AnalyzeDiffCellConc_Gompertz.m" in "AllGrowthCurves_0.05lag_Set1to35.mat"
       [timeVector,od,conc, aucAtConcentration,  controlGompertz, predictedParams] = GC_matching_halfauc_Gompertz_sf(allCCData,drugLabels{i},true, cellConcentration, halfAUC,"MiNoLi wt", timeResolution, timeLimited, gompertz_model, startColor, A600);
  
        %  % Getting Gompertz model of halfAUC
        % boundIdx = boundOrder == cellConcentration;
        % yBlanked = exp(od)*A600;
        % [fit_function, lags(i), growth_rates(i), max_loads(i)] = fit_bacterial_growth_Gompertz_v3(gompertz_model,timeVector, od, false, lower_bound{boundIdx}, upper_bound{boundIdx}, yBlanked', A600);


        controlGompertz{4} = -1*exp(controlGompertz{4})*A600; % making max load negative and back to Abs
        controlGompertz{3} = log(2)/controlGompertz{3}; % making growth rate into generation time
        controlFitFunc = controlGompertz{1}; 
        controlAUC = trapz(timeLimited, exp(controlFitFunc(timeLimited))*A600);
        counter = counter +1;

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
            subplot(rows,cols,counter)
            counter = counter +1;
            hold on;
            
            %plot(x,y,'LineWidth', 2,'Color',cmap(1,:));
            mk  =0;xt =[]; yt =[];
            for k=find(x~=1)'
                mk = mk+1;
                xt = [xt, x(k)]; 
                yt = [yt , y(k)];
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
            l = plot(xFit,yFit, 'LineWidth', 2,'Color',cmap(1,:));

             % plot control 
            plot(controlAUC/controlAUC,controlGompertz{plotConfigs{plotIdx, 5}}, 'o','MarkerSize', 10,'LineStyle', 'none','MarkerFaceColor',[158, 140, 120]/255,'MarkerEdgeColor',[158, 140, 120]/255)
            
            y_pred = polyval(p, xt);

            % get goodness of fit: R^2 
            SS_res = sum((yt - y_pred).^2);         % Residual sum of squares
            SS_tot = sum((yt - mean(yt)).^2);        % Total sum of squares
            RMSE= sqrt(mean((yt - y_pred).^2)); 
            goodnessOfFit(cot, plotIdx) = RMSE; 
            if plotIdx ==1 
                text(0.67,11,  sprintf('RMSE = %0.3f', RMSE))
            elseif plotIdx == 2
                 text(0.67,2.25, sprintf('RMSE = %0.3f', RMSE))
            else
                 text(0.67,-0.1, sprintf('RMSE = %0.3f', RMSE))
            end 
            %xlim([0.5 1])
            xlim([0.2 1])
            ylim(yLimits)
            
            ylabel(yLabel);
            xlabel('AUC');
            box on; grid on;
            axisS = [axisS, gca];
        end
        cot = cot +1;
    end
end
saveas(fig,['SuppFig2_' num2str(figCount) '.jpg'])
        % % % set the figure's Renderer property to 'painters' before saving as EPS/SVG:
        set(fig, 'Renderer', 'painters');
        % % If saving as an EPS file, in R2020a and later, you can use the exportgraphics command and specify the 'ContentType' as 'vector': 
        exportgraphics(fig, ['SuppFig2_' num2str(figCount) '.eps'], 'ContentType', 'vector');