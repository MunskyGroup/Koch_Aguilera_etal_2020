function  [fit_Percent,conditionPercentage,fit_Percent_1G,fit_Percent_2G,fit_Percent_BG] = percentagePerFrame(numberOfSpots,folderName,intensityVector_0,intensityVector_1,plottingCondition,strExperimentalData )
for i =1:numberOfSpots
    StrData(1,i).trajectory_1 = intensityVector_0(i,:);
    StrData(1,i).trajectory_2 = intensityVector_1(i,:);
    % calculating size of trajectories
    StrData(1,i).Length = length(StrData(1,i).trajectory_1);
end
for i =1:numberOfSpots
    if max( StrData(1,i).trajectory_1) >0
        StrData(1,i).trajectory_1 =  StrData(1,i).trajectory_1./max( StrData(1,i).trajectory_1);
    end
    if max( StrData(1,i).trajectory_2) >0
        StrData(1,i).trajectory_2 = StrData(1,i).trajectory_2 ./ max(StrData(1,i).trajectory_2);
    end
end
observationTime = 0; % in seconds
[sim_1G,sim_2G,sim_BG,sim_NT, err_sim_1G,err_sim_2G,err_sim_BG,err_sim_NT]  = simulationsToPercentage (numberOfSpots,StrData,observationTime,plottingCondition);

% Converting percentages to probabilities
sim_prob_1G = sim_1G/100;
sim_prob_2G = sim_2G/100;
sim_prob_BG = sim_BG/100;
sim_prob_NT = sim_NT/100;

% removing zeros from the simulation to avoid conflicts with log
sim_prob_1G(sim_prob_1G==0)= 1e-5;sim_prob_2G(sim_prob_2G==0)= 1e-5; sim_prob_BG(sim_prob_BG==0)= 1e-5; sim_prob_NT(sim_prob_NT==0)= 1e-5;

%% Experimental Data in percentage
exp_1G = strExperimentalData.Percentage_Cap;
exp_2G = strExperimentalData.Percentage_Ires;
exp_BG = strExperimentalData.Percentage_Cap_Ires;
exp_NT = 100-(exp_1G+exp_2G+exp_BG);
% exp_NT = strExperimentalData.Percentage_NT;
expSpots=strExperimentalData.NumberForPercentage;

err_exp_1G = strExperimentalData.Percentage_err_Cap;
err_exp_2G = strExperimentalData.Percentage_err_Ires;
err_exp_BG = strExperimentalData.Percentage_err_Cap_Ires;
err_exp_NT = strExperimentalData.Percentage_err_NT;
% error propagation
% err_exp_NT = sqrt(err_exp_1G^2+err_exp_2G^2+err_exp_BG^2);

%% Evaluating objective function
fit_Percent_1G = (exp_1G/100)*expSpots * -log(sim_prob_1G); % CAP
fit_Percent_2G = (exp_2G/100)*expSpots * -log(sim_prob_2G); % IRES
fit_Percent_BG = (exp_BG/100)*expSpots * -log(sim_prob_BG); % CAP-IRES
fit_Percent_NT = (100-(exp_1G+exp_2G+exp_BG))/100* expSpots * -log(1-sim_prob_1G-sim_prob_2G-sim_prob_BG); % non Translating
% fit_Percent_NT = (exp_NT/100)* expSpots * -log(sim_prob_NT); % non Translating
fit_Percent = (fit_Percent_1G + fit_Percent_2G + fit_Percent_BG+ fit_Percent_NT); % Sum of all
fit_Percent(isreal(fit_Percent)==0) =inf;
conditionPercentage =1;

if plottingCondition ==1
    width = 1.8;
    height = 1.5;
    scaleFig =2;
    font_gca = 12;
    font_labels =16;
    font_legend =10;
    BoxLineWidth =2;
    y = [ exp_1G, sim_1G ;exp_2G, sim_2G; exp_BG, sim_BG; exp_NT, sim_NT];
    std_dev = [ err_exp_1G,err_sim_1G; err_exp_2G,err_sim_2G;  err_exp_BG,err_sim_BG; err_exp_NT, err_sim_NT ];
    figure('visible', 'off');
    fig1= gcf;
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0, 0,width,height]*scaleFig;
    num = 4; %number of different subcategories
    c = 1:num;
    axes1 = axes;
    hold on
    box on
    % Bar(s)
    barColor1 = [0.1 0.1 0.1];     %black
    barColor2 = [0.8 0.8 0.8];     %gray
    for i = 1:num
        bar(c(i)-0.2,y(i,1),0.3,'FaceColor',barColor1);
        bar(c(i)+0.2,y(i,2),0.3,'FaceColor',barColor2);
    end
    errH1 = errorbar(c-0.2,y(:,1),std_dev(:,1),'.','Color','k');
    errH2 = errorbar(c+0.2,y(:,2),std_dev(:,2),'.','Color','k');
    errH1.LineWidth = 0.7;
    errH2.LineWidth = 0.7;
    errH1.Color = [0.5 0.5 0.5];
    errH2.Color = [0.5 0.5 0.5];
    box on
    [~, hobj, ~, ~] =    legend({'Data', 'Model'},'Fontsize',font_legend);
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',0.1);
    set(axes1,'Xlim',[0.5 4.5]);
    set(axes1,'XTick',[1 1.5 2 2.5 3 3.5 4 4.5],'XTickLabel',...
        {'CAP',' ','IRES',' ','Cap+IRES',' ', 'NT', ' '});
    ylim([0, 80])
    set (gca ,'FontSize',font_gca,'FontName', 'Arial','linewidth',BoxLineWidth);
    ylabel('Translation %','FontSize',font_labels);
    nameplot = horzcat('Percent');
    print('-dpng','-r300',nameplot)
    movefile(horzcat(nameplot, '.png'),horzcat(folderName),'f');
    % Saving all data used for the plot
    save('data_Figure_5_B.mat','exp_1G', 'sim_1G' ,'exp_2G', 'sim_2G', 'exp_BG', 'sim_BG', 'exp_NT', 'sim_NT', ...
        'err_exp_1G','err_sim_1G', 'err_exp_2G','err_sim_2G',  'err_exp_BG','err_sim_BG', 'err_exp_NT','err_sim_NT',...
        'width','height','scaleFig','font_labels','font_legend','font_gca','BoxLineWidth','folderName')
    movefile('data_Figure_5_B.mat',folderName,'f');
    
end

end