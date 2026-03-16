% Example script that extracts apical areas, contact lengths and
% ms2 intensities for each cluster

close all
clear all

%% Set Parameters and Select Window Folder
cluster_name = 'Window1Cluster4';
WindowFolder = uigetdir; 
cd(WindowFolder)
res = 3.3453;
Interval =  24.36223; % seconds
mkdir([cluster_name,'/Absolute'])
mkdir([cluster_name,'/ROIs'])

%% load cell table
TACell = readmatrix('TA cell area/cell data.csv');

%% load bond table
TABond = readmatrix('TA cell area/bond data.csv');

%% load cluster spot table
[fn pt] = uigetfile('*.csv');
SpotTracks = readmatrix(fn);

%% PICK WINDOW AND CLUSTER

                 %%%%%%%%%%% WINDOW 1 %%%%%%%%%%%

%%  Cluster 1
% NB = 10782642;
% SC = [1723024 15321000 6893646 14501438 11560212 11895081 1763659];
% SC_Bond = [11401580 657217 9164000 3721360 10765573 14267099 2192091];
% SCPlusOne = [];
% SC_COMBINED = [SC SCPlusOne];
% ClusterFrames = 12:68;
% SignalStart = [5 7 14];
% ScSpotTrackIds = [0 2 4 18 NaN 3 25];

%%  Cluster 2
% NB = 11560212;
% SC = [15321000 10782642 1071479 14501438 16412547];
% SC_Bond = [4957480 10765573 11068878 3783697 8683198];
% SCPlusOne = [];
% SC_COMBINED = [SC SCPlusOne];
% ClusterFrames = 14:73;
% % SignalStart = [5 7 14];
% % ScSpotTrackIds = [Cl1_2 NaN 2 Cl1_18 NaN]; %RENAMED CL1_2 TO 18
% ScSpotTrackIds = [3 NaN 2 18 NaN];


%%  Cluster 3
% NB = 14535140;
% SC = [13565656 11832885 2359369 8667400 13483239];
% SC_Bond = [10347807 3590701 4492328 12194670 11924701];
% SCPlusOne = [];
% SC_COMBINED = [SC SCPlusOne];
% ClusterFrames = 13:65;
% % SignalStart = [5 7 14];
% ScSpotTrackIds = [0 8 NaN NaN NaN];

% 
%%  Cluster 4
NB = 16412547;
% SC = [12029041 14741985 (15-18) + 2536253 (19-56) + 13536444 (58-75) 11560212 6695658];
SC = [12029041 123456789 11560212 6695658];
% SC_Bond = [12693403 6027875 + 2338030 + 1100751 8683198 16604277];
SC_Bond = [12693403 123456789 8683198 16604277];
SCPlusOne = [];
SC_COMBINED = [SC SCPlusOne];
ClusterFrames = 16:76;
% SignalStart = [5 7 14];
ScSpotTrackIds = [2 4 NaN NaN];

if cluster_name == 'Window1Cluster4'

    idChange1 = find(TACell(:,71) == 14741985);
    framesChange1 = TACell(idChange1,2)+1;

    for j = 16:19
        idFrameChange = find(framesChange1==j);
        TACell(idChange1(idFrameChange),71) = 123456789;
    end

    idChange2 = find(TACell(:,71) == 2536253);
    framesChange2 = TACell(idChange2,2)+1;

    for j = 20:57
        idFrameChange = find(framesChange2==j);
        TACell(idChange2(idFrameChange),71) = 123456789;
    end

    idChange3 = find(TACell(:,71) == 13536444);
    framesChange3 = TACell(idChange3,2)+1;

    for j = 59:76
        idFrameChange = find(framesChange3==j);
        TACell(idChange3(idFrameChange),71) = 123456789;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    idChange1Bond = find(TABond(:,24) == 6027875);
    framesChange1Bond = TABond(idChange1Bond,2)+1;

    for j = 16:19
        idFrameChange = find(framesChange1Bond==j);
        TABond(idChange1Bond(idFrameChange),24) = 123456789;
    end

    idChange2Bond = find(TABond(:,24) == 2338030);
    framesChange2Bond = TABond(idChange2Bond,2)+1;

    for j = 20:57
        idFrameChange = find(framesChange2Bond==j);
        TABond(idChange2Bond(idFrameChange),24) = 123456789;
    end

    idChange3Bond = find(TABond(:,24) == 1100751);
    framesChange3 = TABond(idChange3Bond,2)+1;

    for j = 59:76
        idFrameChange = find(framesChange3==j);
        TABond(idChange3Bond(idFrameChange),24) = 123456789;
    end

end


%% Get data from cell tracking

idNB = find(TACell(:,71) == NB);
localidNB = TACell(idNB,3);
framesNB = TACell(idNB,2)+1;
areaNB = TACell(idNB,4);
nbdur = framesNB(end) - framesNB(1) + 1;

for i = 1:length(SC_COMBINED)
    idSC{i} = find(TACell(:,71) == SC_COMBINED(i));
    localidSC{i} = TACell(idSC{i},3);
    framesSC{i} = TACell(idSC{i},2)+1;
    FirstFramesSC(i) = framesSC{i}(1);
    LastFramesSC(i) = framesSC{i}(end);
    areaSC{i} = TACell(idSC{i},4);
end

FramesToAnalyse = ClusterFrames;
MasterMatLocalIds(:,1) = FramesToAnalyse; %Frame no
MasterMatArea(:,1) = FramesToAnalyse; %Frame no

for i = 1:length(framesNB)
    idFrame = find(framesNB(i)==FramesToAnalyse);
    MasterMatLocalIds(idFrame,2) = localidNB(i);
    MasterMatArea(idFrame,2) = areaNB(i);
end

for j = 1:length(SC_COMBINED)
    for i = 1:length(framesSC{j})
        idFrame = find(framesSC{j}(i)==FramesToAnalyse);
        MasterMatLocalIds(idFrame,2+j) = localidSC{j}(i);
        MasterMatArea(idFrame,2+j) = areaSC{j}(i);
    end
end

MasterMatLocalIds(find(MasterMatLocalIds==0))=NaN;
MasterMatArea(find(MasterMatArea==0))=NaN;

%% load bond table

for i = 1:length(SC_Bond)
    idSCBond{i} = find(TABond(:,24) == SC_Bond(i));
    localidSCBond{i} = TABond(idSCBond{i},3);
    framesSCBond{i} = TABond(idSCBond{i},2)+1;
    FirstFramesSCBond(i) = framesSCBond{i}(1);
    LastFramesSCBond(i) = framesSCBond{i}(end);
    lengthSCBond{i} = TABond(idSCBond{i},5);
end

MasterMatLocalIdsBond(:,1) = FramesToAnalyse; %Frame no
MasterMatLengthBond(:,1) = FramesToAnalyse; %Frame no


for j = 1:length(SC_Bond)
    for i = 1:length(framesSCBond{j})
        idFrame = find(framesSCBond{j}(i)==FramesToAnalyse);
        MasterMatLocalIdsBond(idFrame,1+j) = localidSCBond{j}(i);
        MasterMatLengthBond(idFrame,1+j) = lengthSCBond{j}(i);
    end
end

MasterMatLocalIdsBond(find(MasterMatLocalIdsBond==0))=NaN;
MasterMatLengthBond(find(MasterMatLengthBond==0))=NaN;


% % % % % %% Cell Size Line plot
% % % % % 
turnMin = Interval/60;
% % % % % 
LineStyles = {'-k','--b',':m','-.y','-c','--g','--k',':b','-.m','-y','--c',':g'};
% % % % % 

%% Cell Size Heatmap
f3 = figure(3);
sizeCdata = length(FramesToAnalyse);
% cdata = NaN(length(SC)+1,sizeCdata)

cdata = MasterMatArea(:,2:end)';

for i = 1:sizeCdata
    xvalues{i} = num2str(i * turnMin);
end

yvalues{1} = 'NB';
for i = 1:length(SC)
    yvalues{i+1} = ['Cell',num2str(i)];
end
for i = 1:length(SCPlusOne)
    yvalues{length(SC)+1+i} = ['SurrCell',num2str(i)]
end

h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none');
h.GridVisible = 'off';
h.Title = ['CellAreaMovie2',num2str(cluster_name)];
h.XLabel = 'Frame';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
clim([min(cdata(:)) max(cdata(:))])

%% Absolute Bond Length Heatmap

f4 = figure(4);

cdataBond = MasterMatLengthBond(:,2:end)';

for i = 1:length(SC_Bond)
    yvaluesBond{i} = ['BondCell',num2str(i)];
end

h = heatmap(xvalues,yvaluesBond,cdataBond,'CellLabelColor','none');
h.GridVisible = 'off';
h.Title = ['BondLengthMovie2',num2str(cluster_name)];
h.XLabel = 'Frame';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
clim([min(cdataBond(:)) max(cdataBond(:))])


%% Get data from Spot tracking

for i = 1:length(ScSpotTrackIds)
    if ~isnan(ScSpotTrackIds(i))
    idOfInt{i} = find((SpotTracks(:,3)==ScSpotTrackIds(i)));
    [B,K] = sort(SpotTracks(idOfInt{i},9));
    idOfInt{i} = idOfInt{i}(K);
    MS2Frames{i} = SpotTracks(idOfInt{i},9);
    MS2MeanInt{i} = SpotTracks(idOfInt{i},13);

    % keep only for frames we have all cells
    idFrameSTART(i) = min(find(MS2Frames{i} >= FramesToAnalyse(1)));
    idFrameEND(i) = max(find(MS2Frames{i} < FramesToAnalyse(end)));
    
    MS2Frames{i} = MS2Frames{i}(idFrameSTART(i):idFrameEND(i));
    MS2MeanInt{i} = MS2MeanInt{i}(idFrameSTART(i):idFrameEND(i));
    end
end

IntensitiesPerFrameMatrix(:,1) = FramesToAnalyse;
for i = 1:length(ScSpotTrackIds)
    if ~isnan(ScSpotTrackIds(i))
        RelevantFrames = intersect(FramesToAnalyse,MS2Frames{i});
        for j = 1:length(RelevantFrames)
            IntensitiesPerFrameMatrix(find(FramesToAnalyse==RelevantFrames(j)),i+1) = ...
                MS2MeanInt{i}(find(MS2Frames{i}==RelevantFrames(j)));
        end
    else IntensitiesPerFrameMatrix(:,i+1) = NaN;
    end
end


%% Spot Tracking Line plot

f5 = figure(5);
for j = 1:length(SC)
    plot([1:length(FramesToAnalyse)] * turnMin,...
        IntensitiesPerFrameMatrix(:,j+1),LineStyles{j},...
        'LineWidth',2); hold on
end
hold on
xlabel('Time (min)')
ylabel('Mean Intensity')

if length(SC) == 6
   legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5','Cell 6')
elseif length(SC) == 7
   yvalues = {'NB','SC1','SC2','SC3','SC4','SC5','SC6','SC7'};
elseif length(SC) == 5
   legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5')
elseif length(SC) == 4
    legend('Cell 1','Cell 2','Cell 3','Cell 4')
elseif length(SC) == 3
    legend('Cell 1','Cell 2','Cell 3')
end

%% Spot Tracking Heatmap

clear cdata xvalues yvalues
f6 = figure(6);
sizeCdata = length(FramesToAnalyse);
% cdata = NaN(length(SC)+1,sizeCdata)

cdata(1,1:length(FramesToAnalyse)) = NaN;
for j = 1:length(SC)
    cdata(j+1,1:length(FramesToAnalyse)) = IntensitiesPerFrameMatrix(:,j+1);
end

for i = 1:sizeCdata
    xvalues{i} = num2str(i* turnMin);
end

if length(SC) == 6
   yvalues = {'NB','SC1','SC2','SC3','SC4','SC5','SC6'};
elseif length(SC) == 7
   yvalues = {'NB','SC1','SC2','SC3','SC4','SC5','SC6','SC7'};
elseif length(SC) == 5
   yvalues = {'NB','SC1','SC2','SC3','SC4','SC5'};
elseif length(SC) == 4
    yvalues = {'NB','SC1','SC2','SC3','SC4'};
elseif length(SC) == 3
    yvalues = {'NB','SC1','SC2','SC3'};
end


h = heatmap(xvalues,yvalues,cdata,'CellLabelColor','none');
h.GridVisible = 'off';
h.Title = ['INTENSITIESMovie2',num2str(cluster_name)];
h.XLabel = 'Frame';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
clim([min(cdata(:)) max(cdata(:))])


save(['MasterMatAreaABSOLUTE_',cluster_name,'.mat'],'MasterMatArea')
save(['MasterMatLengthBondABSOLUTE_',cluster_name,'.mat'],'MasterMatLengthBond')
save(['IntensitiesPerFrameMatrix_',cluster_name,'.mat'],'IntensitiesPerFrameMatrix')
