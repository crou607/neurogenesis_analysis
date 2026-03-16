% Script that analyses and plots the boxplots and correlations for
% NB delamination durations and transcription durations shown in 
% Figures 4 and S3

close all
clear all

%% NB Delamination Duration
load('CystsData.mat')
load('ControlData.mat')
load('CanoeData.mat')
load('BayData.mat')
load('BayControlData.mat')

%% get Cysts data
for i = 1:size(CombinedMatrixCysts,2)
    DelamDurationCyst(i) = min(find(isnan(CombinedMatrixCysts(:,i))))-1;
end
DelamDurationCyst = DelamDurationCyst * 0.5;

TranscriptionDurationCysts = TranscriptionDurationCysts * 0.5;
%% get Control data
for i = 1:size(ControlNBData,2)
    DelamDurationControl(i) = min(find(isnan(ControlNBData(:,i))))-1;
end
DelamDurationControl = DelamDurationControl * 0.5;

ControlNCData(find(ControlNCData<=0))=NaN;
for i = 1:size(ControlNBData,2)
    idNCs = find(ControlClustering==i);
    for j = 1:length(idNCs)
        start_point_individual(j) = min(find(~isnan(ControlNCData(:,idNCs(j)))));
        end_point_individual(j) = max(find(~isnan(ControlNCData(:,idNCs(j)))));
    end
    cluster_start_point = min(start_point_individual);
    cluster_end_point = max(end_point_individual);
    TranscriptionDurationControl(i) = cluster_end_point - cluster_start_point;
    clear start_point_individual end_point_individual cluster_start_point cluster_end_point
end
TranscriptionDurationControl= TranscriptionDurationControl * 0.5;

%% get Canoe data
for i = 1:size(CanoeNBData,2)
    DelamDurationCanoe(i) = min(find(isnan(CanoeNBData(:,i))))-1;
end
DelamDurationCanoe = DelamDurationCanoe * 0.5;

CanoeNCData(find(CanoeNCData<=0))=NaN;
for i = 1:size(CanoeNBData,2)
    idNCs = find(CanoeClustering==i);
    for j = 1:length(idNCs)
        start_point_individual(j) = min(find(~isnan(CanoeNCData(:,idNCs(j)))));
        end_point_individual(j) = max(find(~isnan(CanoeNCData(:,idNCs(j)))));
    end
    cluster_start_point = min(start_point_individual);
    cluster_end_point = max(end_point_individual);
    TranscriptionDurationCanoe(i) = cluster_end_point - cluster_start_point;
    clear start_point_individual end_point_individual cluster_start_point cluster_end_point
end
TranscriptionDurationCanoe = TranscriptionDurationCanoe * 0.5;
%% get Bay data
for i = 1:size(BayNBData,2)
    DelamDurationBay(i) = min(find(isnan(BayNBData(:,i))))-1;
end
% DelamDurationBay = DelamDurationBay * 0.5;

BayNCData(find(BayNCData<=0))=NaN;

for i = 1:size(BayNBData,2)
    idNCs = find(BayClustering==i);
    for j = 1:length(idNCs)
        start_point_individual(j) = min(find(~isnan(BayNCData(:,idNCs(j)))));
        end_point_individual(j) = max(find(~isnan(BayNCData(:,idNCs(j)))));
    end
    cluster_start_point = min(start_point_individual);
    cluster_end_point = max(end_point_individual);
    TranscriptionDurationBay(i) = cluster_end_point - cluster_start_point;
    clear start_point_individual end_point_individual cluster_start_point cluster_end_point
end
TranscriptionDurationBay = TranscriptionDurationBay * 0.5;

%% get Bay Control data
for i = 1:size(BayControlNBData,2)
    DelamDurationBayControl(i) = min(find(isnan(BayControlNBData(:,i))))-1;
end
% DelamDurationBayControl = DelamDurationBayControl * 0.5;

BayControlNCData(find(BayControlNCData<=0))=NaN;

for i = 1:size(BayControlNBData,2)
    idNCs = find(BayControlClustering==i);
    for j = 1:length(idNCs)
        start_point_individual(j) = min(find(~isnan(BayControlNCData(:,idNCs(j)))));
        end_point_individual(j) = max(find(~isnan(BayControlNCData(:,idNCs(j)))));
    end
    cluster_start_point = min(start_point_individual);
    cluster_end_point = max(end_point_individual);
    TranscriptionDurationBayControl(i) = cluster_end_point - cluster_start_point;
    clear start_point_individual end_point_individual cluster_start_point cluster_end_point
end
TranscriptionDurationBayControl = TranscriptionDurationBayControl * 0.5;



%% Delam duration plots

DelamDurationCyst(3) = [];
TranscriptionDurationCysts(3) =[];

xpos = [1 1.5 2];
Y = [DelamDurationControl DelamDurationCyst DelamDurationCanoe];
  
X = [xpos(1)*repmat(1,size(DelamDurationControl,2),1)' xpos(2)*repmat(1,size(DelamDurationCyst,2),1)' ...
    xpos(3)*repmat(1,size(DelamDurationCanoe,2),1)'];


fig1 = figure(1)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
% scatter(ones(size(DelamDurationWT)).*(xpos(1)+(rand(size(DelamDurationWT))-0.5)/10),DelamDurationWT,15,'k','filled')
scatter(ones(size(DelamDurationControl)).*(xpos(1)+(rand(size(DelamDurationControl))-0.5)/10),DelamDurationControl,15,'k','filled')
scatter(ones(size(DelamDurationCyst)).*(xpos(2)+(rand(size(DelamDurationCyst))-0.5)/10),DelamDurationCyst,15,'k','filled')
scatter(ones(size(DelamDurationCanoe)).*(xpos(3)+(rand(size(DelamDurationCanoe))-0.5)/10),DelamDurationCanoe,15,'k','filled')

xticks([xpos(1) xpos(2) xpos(3)])
xticklabels({'Control','Cysts','Canoe'})
title('Delamination duration from transcription onset')
ylabel('Duration (mins)')
ylim([0 70])
saveas(fig1,'NBDelaminatonDurationBoxplotsAll.pdf')
%ttests
[t,p] = ttest2(DelamDurationControl,DelamDurationCyst)
[t,p] = ttest2(DelamDurationControl,DelamDurationCanoe)

mean(DelamDurationControl)
std(DelamDurationControl)

mean(DelamDurationCyst)
std(DelamDurationCyst)

mean(DelamDurationCanoe)
std(DelamDurationCanoe)
%% Transcription Duration
xpos = [1 1.5 2];
Y = [TranscriptionDurationControl TranscriptionDurationCysts TranscriptionDurationCanoe];

X = [xpos(1)*repmat(1,size(TranscriptionDurationControl,2),1)' xpos(2)*repmat(1,size(TranscriptionDurationCysts,2),1)'...
    xpos(3)*repmat(1,size(TranscriptionDurationCanoe,2),1)'];


fig2 = figure(2)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
scatter(ones(size(TranscriptionDurationControl)).*(xpos(1)+(rand(size(TranscriptionDurationControl))-0.5)/10),TranscriptionDurationControl,15,'k','filled')
scatter(ones(size(TranscriptionDurationCysts)).*(xpos(2)+(rand(size(TranscriptionDurationCysts))-0.5)/10),TranscriptionDurationCysts,15,'k','filled')
scatter(ones(size(TranscriptionDurationCanoe)).*(xpos(3)+(rand(size(TranscriptionDurationCanoe))-0.5)/10),TranscriptionDurationCanoe,15,'k','filled')

xticks([xpos(1) xpos(2) xpos(3)])
xticklabels({'Control','Cysts','Canoe'})
title('Transcription Onset (from first MS2 spot appearing to last spot disappearing')
ylabel('Duration (mins)')
ylim([0 70])
saveas(fig2,'TranscriptionDurationBoxplotsAll.pdf')
%ttests
[t,p] = ttest2(TranscriptionDurationControl,TranscriptionDurationCysts)
[t,p] = ttest2(TranscriptionDurationControl,TranscriptionDurationCanoe)

mean(TranscriptionDurationControl)
std(TranscriptionDurationControl)

mean(TranscriptionDurationCysts)
std(TranscriptionDurationCysts)

mean(TranscriptionDurationCanoe)
std(TranscriptionDurationCanoe)


% %% Transcription levels
% xpos = [1 1.5];
% Y = [TranscriptionTracesMaxControl TranscriptionTracesMaxCanoe];
% 
% X = [xpos(1)*repmat(1,size(TranscriptionTracesMaxControl,2),1)' xpos(2)*repmat(1,size(TranscriptionTracesMaxCanoe,2),1)'];
% 
% 
% fig5 = figure(5)
% boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
% hold on
% scatter(ones(size(TranscriptionTracesMaxControl)).*(xpos(1)+(rand(size(TranscriptionTracesMaxControl))-0.5)/10),TranscriptionTracesMaxControl,15,'k','filled')
% scatter(ones(size(TranscriptionTracesMaxCanoe)).*(xpos(2)+(rand(size(TranscriptionTracesMaxCanoe))-0.5)/10),TranscriptionTracesMaxCanoe,15,'k','filled')

%% Correlate NB Delamination duration with Transcription duration

f3 = figure(3)
hold on
scatter(DelamDurationControl,TranscriptionDurationControl,'k','filled')
scatter(DelamDurationCyst,TranscriptionDurationCysts,'r','filled')
scatter(DelamDurationCanoe,TranscriptionDurationCanoe,'b','filled')

xlabel('Delamination Duration (mins)')
ylabel('Transcription Duration (mins)')
legend
axis tight

DelamDurationsCombined = [DelamDurationControl DelamDurationCyst DelamDurationCanoe];
TranscriptionDurationsCombined = [TranscriptionDurationControl TranscriptionDurationCysts TranscriptionDurationCanoe];

% --- Correlation ---
[R,P] = corrcoef(DelamDurationsCombined,TranscriptionDurationsCombined);
r = R(1,2);
p = P(1,2);

% --- Fit regression line ---
coeffs = polyfit(DelamDurationsCombined,TranscriptionDurationsCombined,1);
xFit = linspace(4,60,100);
yFit = polyval(coeffs,xFit);

% Plot regression line
plot(xFit,yFit,'k-','LineWidth',2)

% Add correlation text to plot
text(1,45,sprintf('r = %.2f\np = %.3f',r,p),...
    'FontSize',12,'BackgroundColor','w')

legend({'Control','Cyst','Canoe','Linear fit'},'Location','best')
ylim([0 70])
xlim([0 65])

saveas(f3,'CorrelationsAll.pdf')

%% Delam duration plots with Bay
xpos = [1 1.5 2 2.5 3];
Y = [DelamDurationControl DelamDurationCyst DelamDurationCanoe DelamDurationBayControl DelamDurationBay];
  
X = [xpos(1)*repmat(1,size(DelamDurationControl,2),1)' xpos(2)*repmat(1,size(DelamDurationCyst,2),1)' ...
    xpos(3)*repmat(1,size(DelamDurationCanoe,2),1)' xpos(4)*repmat(1,size(DelamDurationBayControl,2),1)' xpos(5)*repmat(1,size(DelamDurationBay,2),1)'];


fig11 = figure(11)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
% scatter(ones(size(DelamDurationWT)).*(xpos(1)+(rand(size(DelamDurationWT))-0.5)/10),DelamDurationWT,15,'k','filled')
scatter(ones(size(DelamDurationControl)).*(xpos(1)+(rand(size(DelamDurationControl))-0.5)/10),DelamDurationControl,15,'k','filled')
scatter(ones(size(DelamDurationCyst)).*(xpos(2)+(rand(size(DelamDurationCyst))-0.5)/10),DelamDurationCyst,15,'k','filled')
scatter(ones(size(DelamDurationCanoe)).*(xpos(3)+(rand(size(DelamDurationCanoe))-0.5)/10),DelamDurationCanoe,15,'k','filled')
scatter(ones(size(DelamDurationBayControl)).*(xpos(4)+(rand(size(DelamDurationBayControl))-0.5)/10),DelamDurationBayControl,15,'k','filled')
scatter(ones(size(DelamDurationBay)).*(xpos(5)+(rand(size(DelamDurationBay))-0.5)/10),DelamDurationBay,15,'k','filled')

xticks([xpos(1) xpos(2) xpos(3) xpos(4) xpos(5)])
xticklabels({'Control','Cysts','Canoe','BayControl','Bay'})
title('Delamination duration from transcription onset')
ylabel('Duration (mins)')
saveas(fig11,'NBDelaminatonDurationBoxplotsAllWithBay.pdf')
%ttests
[t,p] = ttest2(DelamDurationControl,DelamDurationCyst)
[t,p] = ttest2(DelamDurationControl,DelamDurationCanoe)


%% Transcription Duration with Bay
xpos = [1 1.5 2 2.5 3];
Y = [TranscriptionDurationControl TranscriptionDurationCysts TranscriptionDurationCanoe TranscriptionDurationBayControl TranscriptionDurationBay];

X = [xpos(1)*repmat(1,size(TranscriptionDurationControl,2),1)' xpos(2)*repmat(1,size(TranscriptionDurationCysts,2),1)'...
    xpos(3)*repmat(1,size(TranscriptionDurationCanoe,2),1)' xpos(4)*repmat(1,size(TranscriptionDurationBayControl,2),1)'...
    xpos(5)*repmat(1,size(TranscriptionDurationBay,2),1)'];


fig12 = figure(12)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
scatter(ones(size(TranscriptionDurationControl)).*(xpos(1)+(rand(size(TranscriptionDurationControl))-0.5)/10),TranscriptionDurationControl,15,'k','filled')
scatter(ones(size(TranscriptionDurationCysts)).*(xpos(2)+(rand(size(TranscriptionDurationCysts))-0.5)/10),TranscriptionDurationCysts,15,'k','filled')
scatter(ones(size(TranscriptionDurationCanoe)).*(xpos(3)+(rand(size(TranscriptionDurationCanoe))-0.5)/10),TranscriptionDurationCanoe,15,'k','filled')
scatter(ones(size(TranscriptionDurationBayControl)).*(xpos(4)+(rand(size(TranscriptionDurationBayControl))-0.5)/10),TranscriptionDurationBayControl,15,'k','filled')
scatter(ones(size(TranscriptionDurationBay)).*(xpos(5)+(rand(size(TranscriptionDurationBay))-0.5)/10),TranscriptionDurationBay,15,'k','filled')

xticks([xpos(1) xpos(2) xpos(3) xpos(4) xpos(5)])
xticklabels({'Control','Cysts','Canoe','BayControl','Bay'})
title('Transcription Onset (from first MS2 spot appearing to last spot disappearing')
ylabel('Duration (mins)')
ylim([0 70])
saveas(fig12,'TranscriptionDurationBoxplotsAllWithBay.pdf')
%ttests
[t,p] = ttest2(TranscriptionDurationControl,TranscriptionDurationCysts)
[t,p] = ttest2(TranscriptionDurationControl,TranscriptionDurationCanoe)


%% Correlate NB Delamination duration with Transcription duration with Bay

f13 = figure(13)
hold on
scatter(DelamDurationControl,TranscriptionDurationControl,'k','filled')
scatter(DelamDurationCyst,TranscriptionDurationCysts,'r','filled')
scatter(DelamDurationCanoe,TranscriptionDurationCanoe,'b','filled')
scatter(DelamDurationBayControl,TranscriptionDurationBayControl,'y','filled')
scatter(DelamDurationBay,TranscriptionDurationBay,'m','filled')

xlabel('Delamination Duration (mins)')
ylabel('Transcription Duration (mins)')
legend
axis tight
ylim([0 70])
xlim([0 70])
DelamDurationsCombined = [DelamDurationControl DelamDurationCyst DelamDurationCanoe DelamDurationBayControl DelamDurationBay];
TranscriptionDurationsCombined = [TranscriptionDurationControl TranscriptionDurationCysts TranscriptionDurationCanoe TranscriptionDurationBayControl TranscriptionDurationBay];
saveas(f13,'CorrelationsAllWithBay.pdf')

corrcoef(DelamDurationsCombined,TranscriptionDurationsCombined)


%% Delam duration plots just Bay
xpos = [1 1.5];
Y = [DelamDurationBayControl DelamDurationBay];
  
X = [xpos(1)*repmat(1,size(DelamDurationBayControl,2),1)' xpos(2)*repmat(1,size(DelamDurationBay,2),1)'];


fig111 = figure(111)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
scatter(ones(size(DelamDurationBayControl)).*(xpos(1)+(rand(size(DelamDurationBayControl))-0.5)/10),DelamDurationBayControl,15,'k','filled')
scatter(ones(size(DelamDurationBay)).*(xpos(2)+(rand(size(DelamDurationBay))-0.5)/10),DelamDurationBay,15,'k','filled')

xticks([xpos(1) xpos(2)])
xticklabels({'BayControl','Bay'})
title('Delamination duration from transcription onset')
ylabel('Duration (mins)')
ylim([0 50])

saveas(fig111,'NBDelaminatonDurationBoxplotsAllJustBay.pdf')

[t,p] = ttest2(DelamDurationBayControl,DelamDurationBay)

mean(DelamDurationBayControl)
std(DelamDurationBayControl)
mean(DelamDurationBay)
std(DelamDurationBay)

%% Transcription Duration just Bay
xpos = [1 1.5];
Y = [TranscriptionDurationBayControl TranscriptionDurationBay];

X = [xpos(1)*repmat(1,size(TranscriptionDurationBayControl,2),1)'...
    xpos(2)*repmat(1,size(TranscriptionDurationBay,2),1)'];


fig121 = figure(121)
boxchart(X(:),Y(:),'BoxWidth',0.4, 'BoxEdgeColor','none','BoxMedianLineColor','k')
hold on
scatter(ones(size(TranscriptionDurationBayControl)).*(xpos(1)+(rand(size(TranscriptionDurationBayControl))-0.5)/10),TranscriptionDurationBayControl,15,'k','filled')
scatter(ones(size(TranscriptionDurationBay)).*(xpos(2)+(rand(size(TranscriptionDurationBay))-0.5)/10),TranscriptionDurationBay,15,'k','filled')

xticks([xpos(1) xpos(2)])
xticklabels({'BayControl','Bay'})
title('Transcription Onset (from first MS2 spot appearing to last spot disappearing')
ylabel('Duration (mins)')
ylim([0 75])
saveas(fig121,'TranscriptionDurationBoxplotsAllJustBay.pdf')
%ttests
[t,p] = ttest2(TranscriptionDurationBayControl,TranscriptionDurationBay)

mean(TranscriptionDurationBayControl)
std(TranscriptionDurationBayControl)
mean(TranscriptionDurationBay)
std(TranscriptionDurationBay)

% % % %% Correlate NB Delamination duration with Transcription duration just Bay
% % % f131 = figure(131)
% % % hold on
% % % scatter(DelamDurationBayControl,TranscriptionDurationBayControl,'y','filled')
% % % scatter(DelamDurationBay,TranscriptionDurationBay,'m','filled')
% % % 
% % % xlabel('Delamination Duration (mins)')
% % % ylabel('Transcription Duration (mins)')
% % % legend
% % % axis tight
% % % ylim([0 50])
% % % xlim([0 30])
% % % 
% % % DelamDurationsCombined = [DelamDurationBayControl DelamDurationBay];
% % % TranscriptionDurationsCombined = [TranscriptionDurationBayControl TranscriptionDurationBay];
% % % saveas(f131,'CorrelationsAllJustBay.pdf')
% % % 
% % % [R,P] = corrcoef(DelamDurationsCombined,TranscriptionDurationsCombined)

%%
f131 = figure(131);
hold on

% Scatter plots
scatter(DelamDurationBayControl,TranscriptionDurationBayControl,'y','filled')
scatter(DelamDurationBay,TranscriptionDurationBay,'m','filled')

xlabel('Delamination Duration (mins)')
ylabel('Transcription Duration (mins)')
axis tight


% Combine data
DelamDurationsCombined = [DelamDurationBayControl DelamDurationBay];
TranscriptionDurationsCombined = [TranscriptionDurationBayControl TranscriptionDurationBay];

% --- Correlation ---
[R,P] = corrcoef(DelamDurationsCombined,TranscriptionDurationsCombined);
r = R(1,2);
p = P(1,2);

% --- Fit regression line ---
coeffs = polyfit(DelamDurationsCombined,TranscriptionDurationsCombined,1);
xFit = linspace(15,45,100);
yFit = polyval(coeffs,xFit);

% Plot regression line
plot(xFit,yFit,'k-','LineWidth',2)

% Add correlation text to plot
text(15,65,sprintf('r = %.2f\np = %.3f',r,p),...
    'FontSize',12,'BackgroundColor','w')

legend({'Control','Bay','Linear fit'},'Location','best')
ylim([5 75])
xlim([10 50])
saveas(f131,'CorrelationsAllJustBay.pdf')