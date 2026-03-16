% Script that carries out cross-correlation analyses shown in Figure 5: 
% NB apical area vs NC transcription levels and
% NC apical area vs NC transcription levels


close all
clear all

%% load isotropic data
load("combined_data_for_correlationsNEW.mat")
% for the last 18 columns dt = 30s, for all the rest dt = 25s

% select which population to analyse: 1 for [1:end-18], 2 for
% [end-18+1:end]

% population = 1;
maxLagMinutes = 23;
analysis_number = 3;
% 1 is NB area vs SC Intensity
% 2 is SC Contact Length vs SC Intensity
% 3 is SC Area vs SC Intensity
filter_data = 1;
chopMS2 = 0;

%% 

if analysis_number == 1

    Variable1 = CombinedNBArea./CombinedAreaNormalisation;
    Variable2 = CombinedIntensities - CombinedBackgrounds;
    AnalysisName = 'NB area vs SC Intensity'

elseif analysis_number == 2
    
    % Variable1 = CombinedContactLength./CombinedContactNormalisation;
    Variable2 = CombinedIntensities - CombinedBackgrounds;
    AnalysisName = 'SC Contact Length vs SC Intensity'

    Variable1 = (CombinedContactLength./CombinedContactNormalisation)./...
        (CombinedSCArea./CombinedAreaNormalisation);

else
    Variable1 = CombinedSCArea./CombinedAreaNormalisation;
    Variable2 = CombinedIntensities - CombinedBackgrounds;
    AnalysisName = 'SC Area vs SC Intensity'

end

Variable2(find(Variable2<0)) = 0;
dt = 30/60; % minutes

% Parameters
Variable1 = movmedian(Variable1, 2, 1);
Variable2 = movmedian(Variable2, 2, 1);

%replace 0s in instensities with NaNs to avoid 'fake' correlation values
Variable2(find(Variable2==0)) = NaN;

maxLag = round(maxLagMinutes / dt); % lag in number of timepoints
lags = (-maxLag:maxLag) * dt;
numLags = length(lags);

% Ensure MamNorm and MS2Norm exist and have the same size
numTracks = size(Variable1, 2);

% Preallocate cross-correlation storage
xcorrs = nan(numLags, 1, numTracks); % for compatibility with SEMplotwithNaN
% Compute cross-correlation for each track
for i = 1:numTracks
    c = crosscorr_nan_safe_prachi(Variable1(:, i), Variable2(:, i), maxLag);
    xcorrs(:, 1, i) = c;
end

%% Plotting correlations
fig = figure;

subplot(2,1,1);
hold on;

% Compute mean and SEM
meanCorr = mean(xcorrs(:,1,:), 3, 'omitnan');
semCorr = std(xcorrs(:,1,:), 0, 3, 'omitnan') ./ sqrt(sum(~isnan(xcorrs(:,1,:)), 3));
% Make sure mean and SEM are row vectors
meanCorr = meanCorr(:)';
semCorr = semCorr(:)';
lagsRow = lags(:)';
plotcolor = [193/255 37/255 101/255];
% Plot SEM as shaded area
fill([lagsRow fliplr(lagsRow)], ...
     [meanCorr + semCorr, fliplr(meanCorr - semCorr)], ...
     plotcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Plot mean line
plot(lags, meanCorr, 'Color', plotcolor, 'LineWidth', 2);
% Final touches
yline(0, '--k');
xlabel('lag (min)');
ylabel('correlation');
title('Individual Cross-Correlations + Mean ± SEM');
ylim([-0.5 1]);
box on;
colorbar
clim([-1 1])
xlim([min(lags), max(lags)]);

xlim([-5,5]);
xticks([-5:5])
 
% xlim([-10, 10]);

% xlim([-10, 10]);
ylim([-0.5,0.5])

subplot(2,1,2);
data = squeeze(xcorrs(:,1,:))';  % Transpose so tracks are rows

h = imagesc(lags, 1:numTracks, data); 

cmap = hot;
cmap = [cmap; 0 0 0];        % add white as an extra color
colormap(cmap);

% Color axis stays for [-1,1] range
caxis([-1 1]);

% Mark NaNs so they use the last colormap row
nanMask = isnan(data);
data(nanMask) = max(data(:)) + 1;  % assign NaNs outside the data range
set(h, 'CData', data);

colorbar;
clim([-1 1])

xlabel('lag (min)');
ylabel('Track #');
title('Cross-Correlation Heatmap');

xlim([-5,5]);
xticks([-5:5])
 
% xlim([-10, 10]);

% xlim([-10, 10]);


saveas(fig,['COMBINED_CorrelationHeatmap',AnalysisName,'_','_time_lag',num2str(maxLagMinutes), '.pdf'])
set(gca,'visible','off')
saveas(fig,['COMBINED_CorrelationHeatmap',AnalysisName,'_','_time_lag',num2str(maxLagMinutes), '.png'])
