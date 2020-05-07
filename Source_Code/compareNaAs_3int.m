function [nfit_val_Total,nfit_val_1G, nfit_val_2G, nfit_val_BG] = compareNaAs_3int(param, folderName, strExperimentalData,sim_CAP_ONLY, sim_IRES_ONLY, sim_CI,error_CAP_ONLY, error_IRES_ONLY,error_CI,plottingCondition)
pre_namePlot ='NaAs';
%% Intensity Deffinition
data_1G = strExperimentalData.mean_exp_NaAs_CAP;
data_2G = strExperimentalData.mean_exp_NaAs_IRES;
data_BG = strExperimentalData.mean_exp_NaAs_CAP_IRES;
err_data_1G = strExperimentalData.std_exp_NaAs_CAP;
err_data_2G = strExperimentalData.std_exp_NaAs_IRES;
err_data_BG = strExperimentalData.std_exp_NaAs_CAP_IRES;
ds_expTime = strExperimentalData.timeExp_NaAs;

%% Plotting
if plottingCondition ==1
    
    width = 1.8;
    height = 1.5;
    scaleFig =2;
    font_gca = 14;
    font_labels =20;
    font_legend =8;
    BoxLineWidth =2;
    %% PLOTTING. 2 frames square.  USED IN PAPER.
    gray =[0.7,0.7,0.7];
    dark_yellow = [1,0.8,0.2];
    black = [ 0,0,0];
    magenta = [0.6,0,0.8];
    color_cap = [0, 0.6,0];
    color_ires = [0 0 1];
    color_CI  = [0.7,0.7,0.7];
    clf
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0, width,height]*scaleFig;
    hold on
    % plotting simulations
    for i =1: size(sim_CAP_ONLY,1)
        lineProps.col= {color_cap}; lineProps.width = 0.5;
        A1 = mseb(ds_expTime,sim_CAP_ONLY(i,:),error_CAP_ONLY(i,:),lineProps,1);
        lineProps.col = {color_ires}; lineProps.width = 0.5;
        A2 = mseb(ds_expTime,sim_IRES_ONLY(i,:),error_IRES_ONLY(i,:),lineProps,1);
        lineProps.col = {color_CI}; lineProps.width = 0.5;
        A3 = mseb(ds_expTime,sim_CI(i,:),error_CI(i,:),lineProps,1);
    end
    % plotting data
    h1 = errorbar(ds_expTime, data_1G, err_data_1G, '^', 'MarkerEdgeColor', color_cap,'Color', color_cap,'MarkerFaceColor', color_cap, 'MarkerSize',4, 'LineStyle','none', 'LineWidth',1);
    h2 = errorbar(ds_expTime, data_2G, err_data_2G, 's', 'MarkerEdgeColor',color_ires, 'Color',color_ires, 'MarkerFaceColor',color_ires,'MarkerSize',4, 'LineStyle','none', 'LineWidth',1);
    h3 = errorbar(ds_expTime, data_BG, err_data_BG, 's', 'MarkerEdgeColor',color_CI, 'Color',color_CI, 'MarkerFaceColor',color_CI,'MarkerSize',4, 'LineStyle','none', 'LineWidth',1);
    plot ([0,0],[-0.05,2.5],'-','Color', [0.5,0.5,0.5],'LineWidth',0.5)
    box on
    set(gca,'linewidth',1)
    % LABLES
    xlabel('Time (min)','FontSize',font_labels);
    ylabel('Norm. Total Intensity (a.u.)','FontSize',font_labels);
    ylim([0 2.5])
    xlim([-900/60 3420/60]);
    grid on
    set(gca,'XGrid','off')
    set(gca,'YGrid','off')
    yticks([0 0.50 1.00 1.5 2 2.5 3 ]);
    yticklabels({'0', '0.5', '1.0','1.5' , '2' ,'2.5', '3' });
    set (gca ,'FontSize',font_gca, 'FontName', 'Arial','linewidth',BoxLineWidth);
    if param.InhibitorType ==1
        nameplot = horzcat(pre_namePlot,'_TypeI');
    else
        nameplot = horzcat(pre_namePlot,'_TypeII');
    end
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    
end
%% Comparing simulations and experiments.
for i=1:size(sim_CAP_ONLY,1)
    dataPoints = 25;
    nsim_CAP_ONLY = sim_CAP_ONLY(i,:);
    nsim_IRES_ONLY = sim_IRES_ONLY(i,:);
    nsim_CI = sim_CI(i,:);
    fit_val_1G(i)  = sum( ((data_1G(1:dataPoints) - nsim_CAP_ONLY(1:dataPoints)).^2)  ./ err_data_1G(1:dataPoints).^2);
    fit_val_2G(i) =  sum( ((data_2G(1:dataPoints) - nsim_IRES_ONLY(1:dataPoints)).^2)  ./ err_data_2G(1:dataPoints).^2);
    fit_val_BG(i) =  sum( ((data_BG(1:dataPoints) - nsim_CI(1:dataPoints)).^2)  ./ err_data_BG(1:dataPoints).^2);
end
nfit_val_1G = mean(fit_val_1G);
nfit_val_2G = mean(fit_val_2G);
nfit_val_BG = mean(fit_val_BG);
nfit_val_Total = nfit_val_1G+nfit_val_2G+nfit_val_BG;

end
