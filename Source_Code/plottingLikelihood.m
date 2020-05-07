%% plotting Function
function plottingLikelihood (typeStress,matirxWithLikelihoods,matirxWithLikelihoods_SD,mean_prediction_1,sd_prediction_1,mean_prediction_2,sd_prediction_2,folderName,preNamePlot,labels_Model,n_freeParameters)
%% Plotting Normalized
close all;
[val,labelOrder] = sort(n_freeParameters);
N_Models = size(matirxWithLikelihoods,2);
x_values  = linspace (1,N_Models,N_Models);
n_sweep = size(mean_prediction_1,1);

%% locating change in number of free parameters.
diff_change_noParam = diff(val);  diff_change_noParam= [1,diff_change_noParam];
x_change_freeParam = x_values(diff_change_noParam~=0); %x_change_freeParam = [x_change_freeParam, length(n_freeParameters)]% detecting the location of the changes to create rectangles indicating the number of free parameters
lab_text= val(x_change_freeParam);

% labelOrder reorder the matrix to consider the desire order during the plotting.
matirxWithLikelihoods = matirxWithLikelihoods(:,labelOrder);
matirxWithLikelihoods_SD = matirxWithLikelihoods_SD(:,labelOrder);
mean_prediction_1 = mean_prediction_1(:,labelOrder);
sd_prediction_1 = sd_prediction_1(:,labelOrder);
mean_prediction_2 = mean_prediction_2(:,labelOrder);
sd_prediction_2 = sd_prediction_2(:,labelOrder);
minVal = min([min(min(mean_prediction_1)),min(min(mean_prediction_2))]);
mean_prediction_1 = mean_prediction_1-minVal;
mean_prediction_2 = mean_prediction_2-minVal;

% maxYval =max([max(max(mean_prediction_1)),max(max(mean_prediction_2))]);
% minYval =  min([min(min(mean_prediction_1)),min(min(mean_prediction_2))]);
maxYval = 500;
minYval =-10;

markerS = 6;
sizeWide=14;
sizeHeight=2.4;
fontSize =12;
fontLabel =16;
font_legend=10;
textSize =8;

% color_1 = [0, 0.4470, 0.7410];
% color_2 = [0.8500, 0.3250, 0.0980];
% color_3 = [0.9290, 0.6940, 0.1250];
% color_4 = [0.4940, 0.1840, 0.5560];
% color_5 = [0, 0, 0];

color_1 = [0, 0.4470, 0.7410;...
0.8500, 0.3250, 0.0980;...
1, 0.2, 0.3980;...
0.9290, 0.6940, 0.1250;...
0.4940, 0.1840, 0.5560;...
0.5 0 0.5;...
0, 0, 0
1 0 1];


LinWidth =3;

figure('visible', 'off');
fig1= gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0, 0, sizeWide, sizeHeight];
hold on
color_vector_1 = colormap(flipud(lines(n_sweep+2)));
colors = get(gca,'colororder');
colors(8,:) =[0.9290, 0.6940, 0.1250];
%[x y w h]
colorRec = [0.9, 0.9, 0.9];
rectangle('Position',[0 minYval 5.5 maxYval+10],'FaceColor',colorRec,'LineWidth',0.5,'EdgeColor',[0.5,0.5,0.5]);
rectangle('Position',[7.5 minYval 2 maxYval+10],'FaceColor',colorRec,'LineWidth',0.5,'EdgeColor',[0.5,0.5,0.5]);
colorRec = [0, 1, 0] ;
rectangle('Position',[5.5 minYval 1 maxYval+10],'FaceColor',colorRec,'LineWidth',0.5);
colorRec = [1, 1,1] ;

rectangle('Position',[0 400 14 100],'FaceColor',colorRec,'LineWidth',0.5,'EdgeColor',[0.5,0.5,0.5]);


% text for the number of parameters
positionText_label =[1.4, 4.2, 8.2,   11.2, 12.6 ,13.5];
for kk=1: 1:length(x_change_freeParam)
    txt = [ num2str(lab_text(kk)) ' Param'];
    text(positionText_label(kk),maxYval-50,txt,'FontSize',textSize)
end

% plotitng the fit and error lines
Sel_Ts = plot ([1,14],[100,100],'--r', 'LineWidth',2);
LinWidth =3.5;
c=1;

%% Cross-validation for chemical stresses
for k =1: 3%n_sweep-1
     mean_prediction_1(:,matirxWithLikelihoods(4,:)>100)=inf;
     mean_prediction_2(:,matirxWithLikelihoods(4,:)>100)=inf;
%     mean_prediction_3(:,matirxWithLikelihoods(4,:)>100)=inf;
 %   color1 = {color_vector_1(c,:)}; %lineProps.width = LinWidth;
    L_p1{k}= errorbar (x_values+0.05*k+0.01,mean_prediction_1(k,:),sd_prediction_1(k,:),'o','MarkerSize',markerS,'MarkerFaceColor',color_1(c,:), 'Color',color_1(c,:),'LineWidth',1); % [0.2,0.2,0.2]
    c=c+1;
 %   color2 = {color_vector_1(c,:)}; %lineProps.width = LinWidth;
    L_p2{k}=errorbar (x_values+0.05*k+0.01,mean_prediction_2(k,:),sd_prediction_2(k,:),'^','MarkerSize',markerS,'MarkerFaceColor',color_1(c,:), 'Color',color_1(c,:),'LineWidth',1);
    c=c+1;
    %lineProps.col = {color_vector_1(c,:)}; lineProps.width = LinWidth;
    %L_p3{k}=errorbar (x_values+0.1*k+0.066,mean_prediction_3(k,:),sd_prediction_3(k,:),'*','MarkerSize',markerS,'MarkerFaceColor',[0.2,0.2,0.2], 'LineWidth',1);
    c=c+1;
end
LinWidth =4;
markerS = 8;

% color3 = {color_vector_1(c,:)};
L_obj = errorbar(x_values-0.1,matirxWithLikelihoods(4,:),matirxWithLikelihoods_SD(4,:),'s','MarkerSize',markerS,'MarkerEdgeColor','k','MarkerFaceColor',[1,0,0],'Color',[0 0 0], 'LineWidth',1);
% label for the selected model
txt = ['Selected Model'];
h=text(5.6,120,txt,'FontSize',textSize);
set(h,'Rotation',90);
if typeStress ==1
    lgd= legend([L_p1{1},L_p1{2},L_p1{3},L_p2{1},L_p2{2},L_p2{3}, L_obj,Sel_Ts], {'L_{NaAs-k_{ON-C}-0%}','L_{NaAs-k_{ON-C}-33%}','L_{NaAs-k_{ON-C}-67%}','L_{NaAs-k_{INIT-C}-0%}','L_{NaAs-k_{INIT-C}-33%}','L_{NaAs-k_{INIT-C}-67%}','L_{Optimized}','Selection threshold'});
else
    lgd= legend([L_p1{1},L_p1{2},L_p1{3},L_p2{1},L_p2{2},L_p2{3}, L_obj,Sel_Ts], {'L_{DTT-k_{ON-C}-0%}','L_{DTT-k_{ON-C}-33%}','L_{DTT-k_{ON-C}-67%}','L_{DTT-k_{INIT-C}-0%}','L_{DTT-k_{INIT-C}-33%}','L_{DTT-k_{INIT-C}-67%}','L_{Optimized}','Selection threshold'});
end
set(lgd,'FontSize',font_legend, 'Location', 'eastoutside');


% plot the lines for different models
for kk=3: 1:length(x_change_freeParam)
    if mod(kk,2) ==0; colorRec = [1,1,1] ;else; colorRec = [0.95, 0.95, 0.95];end
    plot ([x_change_freeParam(kk-1)-0.5,x_change_freeParam(kk-1)-0.5], [minYval,maxYval],'k','LineWidth',2);
end
plot ([x_change_freeParam(kk)-0.5,x_change_freeParam(kk)-0.5], [minYval,maxYval],'k','LineWidth',2);


xlim([1-0.2, N_Models+0.2])
ylim([minYval, maxYval])
box on
set(gca,'linewidth',1)
xticks([1:1:14]);
xticklabels(labels_Model(labelOrder))
set (gca ,'FontSize',fontSize,'FontName', 'Arial');
xlabel('Models','FontSize',fontLabel)
ylabel('-log(Likelihood)','FontSize',fontLabel)
%set(gca, 'YScale', 'log')
nameplot = [preNamePlot,'_Lin'];
print('-dtiff','-r300',nameplot)
movefile(horzcat(nameplot, '.tif'),horzcat(folderName),'f');

end