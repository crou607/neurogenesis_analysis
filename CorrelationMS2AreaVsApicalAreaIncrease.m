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
% Variable2(find(Variable2==0)) = NaN;

Variable2(find(isnan(Variable2))) = 0;

%% calculate area under ms2 curve
Variable2 = cumtrapz(Variable2,1);


%%
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


% figure
for i = 1 : size(Variable1,2)

    CorrVar1 = Variable1(:,i);
    CorrVar2 = Variable2(:,i);

    idNan1 = isnan(CorrVar1);
    CorrVar1(idNan1) = [];
    CorrVar2(idNan1) = [];
    idNan2 = isnan(CorrVar2);
    CorrVar1(idNan2) = [];
    CorrVar2(idNan2) = [];
    
    if ~isempty(CorrVar1) && ~isempty(CorrVar2)
        corrVar1Var2(i) = corr(CorrVar1,CorrVar2);
        % scatter(CorrVar1,CorrVar2,'k','filled')
        % hold on
        corrVar1Var2Grad(i) = corr(diff(CorrVar1),CorrVar2(1:end-1));
    else
        corrVar1Var2(i) = NaN;
    end

    clear CorrVar1 CorrVar2 idNan1 idNan2

end

figure
%boxchart(corrVar1Var2)
violinplot(corrVar1Var2)
hold on
scatter(ones(size(corrVar1Var2)).*(1+(rand(size(corrVar1Var2))-0.5)/2),corrVar1Var2,'k','filled')
nanmean(corrVar1Var2)
ylim([-1.2 1.2])
std(corrVar1Var2)
mean(corrVar1Var2)

figure
[f,xf] = kde(corrVar1Var2,Bandwidth=0.01);
hold on
scatter(ones(size(corrVar1Var2)).*(1+(rand(size(corrVar1Var2))-0.5)/5),corrVar1Var2,'k','filled')
nanmean(corrVar1Var2)
ylim([-2 2])

violinplot(EvaluationPoints=xf,DensityValues=f)