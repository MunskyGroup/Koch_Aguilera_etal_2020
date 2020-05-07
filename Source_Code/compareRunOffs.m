function [nfit_val_total, nfit_val_1G, nfit_val_2G] = compareRunOffs(folderName,strExperimentalData,mean_sim_1G, mean_sim_2G, err_sim_1G, err_sim_2G,plottingCondition)
%% Intensity Deffinition
data_1G = strExperimentalData.mean_exp_HT_CAP;
data_2G = strExperimentalData.mean_exp_HT_IRES;
err_data_1G = strExperimentalData.std_exp_HT_CAP;
err_data_2G = strExperimentalData.std_exp_HT_IRES;
ds_expTime = strExperimentalData.timeExp_Harringtonine;
dataTrajectories_CAP = strExperimentalData.trajectories_exp_CAP;
dataTrajectories_IRES =strExperimentalData.trajectories_exp_IRES;

%% Plotting
if plottingCondition ==1
    width = 1.8;
    height = 1.5;
    scaleFig =2;
    font_gca = 14;
    font_labels =20;
    font_legend =10;
    BoxLineWidth =2;
    %% PLOTTING. 2 frames square.  USED IN PAPER.
    gray =[0.7,0.7,0.7];
    color_cap = [0, 0.6,0];
    color_ires = [0 0 1];
    black = [ 0,0,0];
    %CAP plot
    pre_namePlot = 'HT_CAP';
    figureHT (width,height,scaleFig,1,ds_expTime,dataTrajectories_CAP,color_cap,mean_sim_1G,err_sim_1G,font_labels,font_legend,font_gca,BoxLineWidth,pre_namePlot,folderName)
    % Saving all data used for the plot
    save('data_Figure_5_D.mat','width','height','scaleFig','ds_expTime','dataTrajectories_CAP','color_cap','mean_sim_1G','err_sim_1G','font_labels','font_legend','font_gca','BoxLineWidth','pre_namePlot','folderName')
    movefile('data_Figure_5_D.mat',folderName,'f');
    % IRES plot
    pre_namePlot = 'HT_IRES';
    figureHT (width,height,scaleFig,2,ds_expTime,dataTrajectories_IRES,color_ires,mean_sim_2G,err_sim_2G,font_labels,font_legend,font_gca,BoxLineWidth,pre_namePlot,folderName)
    % Saving all data used for the plot
    save('data_Figure_5_E.mat','width','height','scaleFig','ds_expTime','dataTrajectories_IRES','color_ires','mean_sim_2G','err_sim_2G','font_labels','font_legend','font_gca','BoxLineWidth','pre_namePlot','folderName')
    movefile('data_Figure_5_E.mat',folderName,'f');
end
%% Comparing simulations and experiments.
dataPoints = 35;
n_rep= size(mean_sim_1G,1);
fit_val_1G = zeros(1,n_rep);
fit_val_2G = zeros(1,n_rep);
for i=1:n_rep % this loop iterates for the number of independent repetitions
    nsim_CAP_ONLY = mean_sim_1G(i,:);
    nsim_IRES_ONLY = mean_sim_2G(i,:);
    fit_val_1G(i) =   sum( ((data_1G(1:dataPoints) - nsim_CAP_ONLY(1:dataPoints)).^2)  ./ err_data_1G(1:dataPoints).^2);
    fit_val_2G(i) =   sum( ((data_2G(1:dataPoints) - nsim_IRES_ONLY(1:dataPoints)).^2)  ./ err_data_2G(1:dataPoints).^2);
end
nfit_val_1G = mean(fit_val_1G);
nfit_val_2G = mean(fit_val_2G);
nfit_val_total = nfit_val_1G+nfit_val_2G;

    function figureHT (width,height,scaleFig,geneType,ds_expTime,dataTrajectories,color_plot,mean_sim,err_sim,font_labels,font_legend,font_gca,BoxLineWidth,pre_namePlot,folderName)
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        fig1.PaperPosition = [0, 0,width,height]*scaleFig;
        hold on
        % plotting data
        for i =1: size(dataTrajectories,1)
            h1 =   plot(ds_expTime,dataTrajectories(i,:),'-','Color', [0.6,0.6,0.6],'LineWidth',0.5);
        end
        for i =1: size(mean_sim,1)
            lineProps.col= {color_plot}; lineProps.width = 2;
            A1 = mseb(ds_expTime,mean_sim(i,:),err_sim(i,:),lineProps,1);
        end
        plot ([0,0],[-0.05,1.5],'-','Color', [0,0,0],'LineWidth',1)
        box on
        set(gca,'linewidth',1)
        % LABLES
        xlabel('Time (min)','FontSize',font_labels);
        ylabel('Norm. Total Intensity (a.u.)','FontSize',font_labels);
        ylim([-0.05 1.4])
        xlim([-5 31]);
        if geneType ==1
            lgd= legend([h1, A1.mainLine  ], {'CAP Data','CAP Model'},'FontSize',font_legend);
        else
            lgd= legend([h1, A1.mainLine  ], {'IRES Data','IRES Model'},'FontSize',font_legend);
        end
        grid on
        set(gca,'XGrid','off')
        set(gca,'YGrid','off')
        yticks([0 0.25 0.50 0.75 1.00 1.25]);
        yticklabels({'0.0','0.25','0.5', '0.75', '1.0','1.25'});
        set (gca ,'FontSize',font_gca, 'FontName', 'Arial','linewidth',BoxLineWidth);
        nameplot = horzcat(pre_namePlot);
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end

end
