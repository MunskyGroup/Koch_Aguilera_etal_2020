function  compareNaAs_3int_sweep(param,value_sweep, folderName, strExperimentalData,sim_CAP_ONLY, sim_IRES_ONLY, sim_CI,error_CAP_ONLY, error_IRES_ONLY,error_CI,plottingCondition,SupplementaryFigure)
pre_namePlot ='NaAs_sweep_';

if SupplementaryFigure ==1
    pre_namePlot = ['Sup__',pre_namePlot];
end

%% Intensity Deffinition
data_1G = strExperimentalData.mean_exp_NaAs_CAP;
data_2G = strExperimentalData.mean_exp_NaAs_IRES;
data_BG = strExperimentalData.mean_exp_NaAs_CAP_IRES;
err_data_1G = strExperimentalData.sem_exp_NaAs_CAP;
err_data_2G = strExperimentalData.sem_exp_NaAs_IRES;
err_data_BG = strExperimentalData.sem_exp_NaAs_CAP_IRES;
ds_expTime = strExperimentalData.timeExp_NaAs;
n_sweep = size(sim_CAP_ONLY,1);

%% Plotting
if plottingCondition ==1
    width = 2;
    height = 1.5;
    font_gca = 14; %14
    font_labels =20; %20
    y_lable = 'Norm. Total Intensity (a.u.)';
    
    if SupplementaryFigure ==1
        width = 1;
        height = 1.2;
        font_gca = 12; %14
        font_labels =12; %20
        y_lable = 'Norm. Total Int. (a.u.)';
        pre_fname = 'sim_data_NaAS_';
    else
        pre_fname = 'sim_data_NaAS_';
    end
    
    
    scaleFig =2;
    font_legend =8;
    textSize =8;
    BoxLineWidth =2;
    color_cap = [0, 0.6,0];
    color_ires = [0 0 1];
    color_CI  = [0.3,0.3,0.3];
    color_vector = colormap(flip(autumn(n_sweep+2)));
    color_vector = color_vector(3:end,:);
    color_plot =color_cap;
    makeplot (1,value_sweep,ds_expTime,sim_CAP_ONLY,error_CAP_ONLY,data_1G,err_data_1G,color_vector,color_plot, width,height,scaleFig,font_gca,font_labels,font_legend,BoxLineWidth,y_lable,textSize)
    if param.InhibitorType ==1
        fname=[pre_fname,'CAP_typeI','.mat'];
    elseif  param.InhibitorType ==2
        fname=[pre_fname,'CAP_typeII','.mat'];
    end
    save(fname,'width','height','scaleFig','value_sweep','ds_expTime','sim_CAP_ONLY','error_CAP_ONLY','data_1G','err_data_1G','font_labels','font_legend','font_gca','BoxLineWidth','fname','folderName')
    movefile(fname,folderName,'f');
    
    color_vector = colormap(flip(autumn(n_sweep+2)));
    color_vector = color_vector(3:end,:);
    color_plot =color_ires;
    makeplot (2,value_sweep,ds_expTime,sim_IRES_ONLY,error_IRES_ONLY,data_2G,err_data_2G,color_vector,color_plot, width,height,scaleFig,font_gca,font_labels,font_legend,BoxLineWidth,y_lable,textSize)
    if param.InhibitorType ==1
        fname=[pre_fname,'IRES_typeI','.mat'];
    elseif param.InhibitorType ==2
        fname=[pre_fname,'IRES_typeII','.mat'];
    end
    save(fname,'width','height','scaleFig','value_sweep','ds_expTime','sim_IRES_ONLY','error_IRES_ONLY','data_2G','err_data_2G','font_labels','font_legend','font_gca','BoxLineWidth','fname','folderName')
    movefile(fname,folderName,'f');
    
    color_vector = colormap(flip(autumn(n_sweep+2)));
    color_vector = color_vector(3:end,:);
    color_plot =color_CI;
    makeplot (3,value_sweep,ds_expTime,sim_CI,error_CI,data_BG,err_data_BG,color_vector,color_plot, width,height,scaleFig,font_gca,font_labels,font_legend,BoxLineWidth,y_lable,textSize)
    if param.InhibitorType ==1
        fname=[pre_fname,'CAP_IRES_typeI','.mat'];
    elseif param.InhibitorType == 2
        fname=[pre_fname,'CAP_IRES_typeII','.mat'];
    end
    save(fname,'width','height','scaleFig','value_sweep','ds_expTime','sim_CI','error_CI','data_BG','err_data_BG','font_labels','font_legend','font_gca','BoxLineWidth','fname','folderName')
    movefile(fname,folderName,'f');
%     
end

    function makeplot (geneType,value_sweep,time,sim,error_sim,data,err_data,color_vector,color_plot, width,height,scaleFig,font_gca,font_labels,font_legend,BoxLineWidth,y_lable,textSize)
        value_sweep = round(value_sweep.*100,0);
        geneName = {'CAP','IRES','CAP-IRES'};
        figure('visible', 'off');
        fig1= gcf;
        fig1.PaperUnits = 'inches';
        fig1.PaperPosition = [0, 0, width,height]*scaleFig;
        hold on
        % plotting simulations
        for i =1: size(sim,1)
            lineProps.col= {color_vector(i,:)}; lineProps.width = 0.5;
            A1 = mseb(time,sim(i,:),error_sim(i,:),lineProps,0);
        end
        % plotting data
        h1 = errorbar(time, data, err_data, 's', 'MarkerEdgeColor', color_plot,'Color', color_plot,'MarkerFaceColor', color_plot, 'MarkerSize',4, 'LineStyle','none', 'LineWidth',1);
        plot ([0,0],[-0.05,2.5],'-','Color', [0.5,0.5,0.5],'LineWidth',1)
        box on
        set(gca,'linewidth',1)
        for i =1: size(sim,1)
            text(time(end-5),max([0.07,sim(i,end-3)]),[num2str(value_sweep(i)),'%'],'FontSize',textSize)
        end
        % LABLES
        xlabel('Time (min)','FontSize',font_labels);
        ylabel(y_lable,'FontSize',font_labels);
        xlim([-15 57]);
        ylim([-0.05 2]);
        %lgd = legend([h1,A1.mainLine, h2,A2.mainLine,h3,A3.mainLine  ], {'CAP Data','CAP Model','IRES Data','IRES Model','CI Data','CI Model' },'Location','FontSize',font_legend);
        grid on
        set(gca,'XGrid','off')
        set(gca,'YGrid','off')
        yticks([0 0.50 1.00 1.5 2 2.5 3 ]);
        yticklabels({'0', '0.5', '1.0','1.5' , '2' ,'2.5', '3' });
        set (gca ,'FontSize',font_gca, 'FontName', 'Arial','linewidth',BoxLineWidth);
        if param.InhibitorType ==1
            nameplot = horzcat(pre_namePlot,geneName{geneType},'_TypeI');
        elseif param.InhibitorType ==2
            nameplot = horzcat(pre_namePlot,geneName{geneType},'_TypeII');
        elseif param.InhibitorType ==3
            nameplot = horzcat(pre_namePlot,geneName{geneType},'_TypeIII');
        end
        print('-dpng','-r300',nameplot)
        movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    end
end
