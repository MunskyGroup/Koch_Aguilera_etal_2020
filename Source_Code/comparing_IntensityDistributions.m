function [fit_LL_Distribution,fit_LL_Distribution_CAP,fit_LL_Distribution_IRES] = comparing_IntensityDistributions(Intensity_1G,Intensity_2G,strExperimentalData,folderName,plottingCondition)
sizeInt = size(Intensity_1G,1); % Number of spots
Sampling_numberOfSpots_target =1000; % representative number of spots per cell
Sampling_numberOfSpots = min ([Sampling_numberOfSpots_target, sizeInt]);
n_Botstrapping=50;
edg= 1:2:80;

% deffining a thereshold to remove elements below the number of repetions of the tag in the construct.
threshold_1G = 1; % select a threshold to remove from data.
threshold_2G = 1; % select a threshold to remove from data

% removing elements spots with intensity in both signals.
temp_Intensity_1G = Intensity_1G(:,end); % selecting last element or tp
temp_Intensity_1G(temp_Intensity_1G<threshold_1G)=0; % making zero intensities below detection threshold
temp_Intensity_1G(temp_Intensity_1G>0)=1; % making one intensities above the detection threshold

temp_Intensity_2G = Intensity_2G(:,end); % selecting last element or tp
temp_Intensity_2G(temp_Intensity_2G<threshold_2G)=0; % making zero intensities below detection threshold
temp_Intensity_2G(temp_Intensity_2G>0)=1; % making one intensities above the detection threshold

spots_with_Both_signals =  temp_Intensity_1G + temp_Intensity_2G; % summing both vectores. 0 = no translation, 1= translation in one of the genes. 2= translation in both genes.
spots_with_Both_signals(spots_with_Both_signals<2)=0; % keeping only those elements with values equal to 2.
IndexWith_Both_Intensities = find (spots_with_Both_signals); % finding index with both intensities.

% making zero elements with both intensities. This will pass only spots
% with intensity in CAP-only and IRES-only. this will be used for fitting
% the experimental distributions.

Intensity_1G_lastTimePoint= Intensity_1G(:,end);
Intensity_2G_lastTimePoint= Intensity_2G(:,end);

Intensity_1G_lastTimePoint(IndexWith_Both_Intensities) =0; % making zeros elements with intensities for both genes.
Intensity_2G_lastTimePoint(IndexWith_Both_Intensities) =0; % making zeros elements with intensities for both genes.

intensity_Simulation_1G {n_Botstrapping}= [];
intensity_Simulation_2G {n_Botstrapping}= [];

for bt=1:n_Botstrapping
    %   ii = randi([1, sizeInt], 1, sizeInt);
    counter = 1;
    for ii =1:Sampling_numberOfSpots
        ii= randi(sizeInt);
        vec_intensity_Simulation_1G(counter) =  Intensity_1G_lastTimePoint(ii);
        vec_intensity_Simulation_2G(counter)  =  Intensity_2G_lastTimePoint(ii);
        counter = counter +1;
    end
    vec_intensity_Simulation_1G (vec_intensity_Simulation_1G<threshold_1G)= [];
    vec_intensity_Simulation_2G (vec_intensity_Simulation_2G<threshold_2G)= [];
    
    if isempty (vec_intensity_Simulation_1G) ==1
        vec_intensity_Simulation_1G = zeros(1,sizeInt);
    end
    if isempty (vec_intensity_Simulation_2G) ==1
        vec_intensity_Simulation_2G = zeros(1,sizeInt);
    end
    intensity_Simulation_1G {bt}=vec_intensity_Simulation_1G ;
    intensity_Simulation_2G {bt}=vec_intensity_Simulation_2G ;
    
    %% Compare both intensities.
    [probability_Sim_1G,~]= histcounts (vec_intensity_Simulation_1G,edg,'Normalization','probability');
    [probability_Sim_2G,~] = histcounts (vec_intensity_Simulation_2G,edg,'Normalization','probability');
    probability_Sim_1G(probability_Sim_1G==0)=1e-5; probability_Sim_2G(probability_Sim_2G==0)=1e-5; % removing zeros to avoid conflicts with log function.
    [hist_exp_1G,~] = histcounts(strExperimentalData.Int_cap,edg);
    [hist_exp_2G,~] = histcounts(strExperimentalData.Int_ires,edg);
    LL_CAP(1,bt) =  (-  dot(hist_exp_1G,log(probability_Sim_1G))) ;
    LL_IRES(1,bt) =  (-  dot(hist_exp_2G,log(probability_Sim_2G))) ;
end

fit_LL_Distribution_CAP = mean (LL_CAP);
fit_LL_Distribution_IRES = mean (LL_IRES);
hisData_1G = zeros(n_Botstrapping,length(edg));
hisData_2G = zeros(n_Botstrapping,length(edg));

%% Converting data to distributions
for i=1: n_Botstrapping
    [hisData_1G(i,1:end-1),edg1 ]= histcounts (intensity_Simulation_1G {i},edg,'Normalization','probability');
    [xb_1G,Data_stairs_1G(:,i)] = stairs(edg1,hisData_1G(i,:));
end
Data_stairs_1G= Data_stairs_1G';
mean_stairs_1G = mean(Data_stairs_1G);
std_stairs_1G =  std(Data_stairs_1G);

for i=1: n_Botstrapping
    [hisData_2G(i,1:end-1),edg2] = histcounts (intensity_Simulation_2G {i},edg,'Normalization','probability');
    [xb_2G,Data_stairs_2G(:,i)] = stairs(edg2,hisData_2G(i,:));
end
Data_stairs_2G= Data_stairs_2G';
mean_stairs_2G = mean(Data_stairs_2G);
std_stairs_2G =  std(Data_stairs_2G);

fit_LL_Distribution = (fit_LL_Distribution_CAP + fit_LL_Distribution_IRES);


if plottingCondition==1
    width = 1.8;
    %height = 1.5;
    height = 0.75;
    scaleFig =2;
    font_gca = 14;
    font_labels =16;
    font_legend =10;
    BoxLineWidth =2;
    maxY =0.5;
    color_Ires = [0 0 1];
    color_Cap = [0, 0.6,0];
    exp_colorPlot = [0 0 0];
    
    %% Plotting CAP
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, width, height].*scaleFig;
    hold on
    % ploting experimental data
    lineProps.col= {color_Cap};
    h_Sim = mseb(xb_1G', mean_stairs_1G,std_stairs_1G,lineProps,0);
    [hisData_norm_1G] =histcounts (strExperimentalData.Int_cap,edg,'Normalization','probability');
    hisData_norm_1G(end+1)=0;
    h_data= stairs(edg,hisData_norm_1G, 'Color' ,exp_colorPlot,'LineWidth',2);
    xlim([1, 50])
    ylim([0, maxY])
    lgd= legend([h_data, h_Sim.mainLine], {'CAP Data','CAP Model'});
    set(lgd,'FontSize',font_legend);
    box on
    grid off
    set (gca ,'FontSize',font_gca, 'FontName', 'Arial','linewidth',BoxLineWidth);
    xlabel('Intensity (ump)','FontName', 'Arial','FontSize',font_labels);
    ylabel('Probability','FontName', 'Arial','FontSize',font_labels);
    nameplot = horzcat('Dist_CAP');
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    % Saving all data used for the plot
    save('data_Figure_5_C_cap.mat','hist_exp_1G','xb_1G','Data_stairs_1G','mean_stairs_1G','std_stairs_1G','hisData_norm_1G','width','height','scaleFig','font_labels','font_legend','font_gca','BoxLineWidth','folderName')
    movefile('data_Figure_5_C_cap.mat',folderName,'f');
    
    %% Plotting IRES
    figure('visible', 'off');
    fig2= gcf;
    fig2.PaperUnits = 'inches';
    fig2.PaperPosition = [0, 0, width, height].*scaleFig;
    hold on
    exp_colorPlot = [0 0 0];
    % ploting experimental data
    lineProps.col= {color_Ires};
    h_Sim = mseb(xb_2G', mean_stairs_2G,std_stairs_2G,lineProps,0);
    [hisData_norm_2G] =histcounts (strExperimentalData.Int_ires,edg,'Normalization','probability');
    hisData_norm_2G(end+1)=0;
    h_data= stairs(edg,hisData_norm_2G, 'Color' ,exp_colorPlot,'LineWidth',2);
    xlim([1, 50])
    ylim([0, maxY])
    lgd= legend([h_data, h_Sim.mainLine], {'IRES Data','IRES Model'});
    set(lgd,'FontSize',font_legend);
    box on
    grid off
    set (gca ,'FontSize',font_gca, 'FontName', 'Arial','linewidth',BoxLineWidth);
    xlabel('Intensity (ump)','FontName', 'Arial','FontSize',font_labels);
    ylabel('Probability','FontName', 'Arial','FontSize',font_labels);
    nameplot = horzcat('Dist_IRES');
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    % Saving all data used for the plot
    save('data_Figure_5_C_ires.mat','hist_exp_2G','xb_2G','Data_stairs_2G','mean_stairs_2G','std_stairs_2G','hisData_norm_2G','width','height','scaleFig','font_labels','font_legend','font_gca','BoxLineWidth','folderName')
    movefile('data_Figure_5_C_ires.mat',folderName,'f');
        
end

end