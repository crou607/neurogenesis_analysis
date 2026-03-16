% Example script that extracts pre-transcription apical areas for NB vs 
% non-NB analysis shown in Figure 3

% Inputs are MasterMatAreaABSOLUTE, MasterMatLengthBondABSOLUTE and 
% IntensitiesPerFrameMatrix .mats obtained from GetSizeAndBondTables.m

close all
clear all

movie_number = '3';
window_number = '1';

if str2double(movie_number) == 9 && str2double(window_number) == 1

    all_clusters_in_window = [1 2 3];
    Resolution = 4.3957;

elseif str2double(movie_number) == 9 && str2double(window_number) == 2
    
    all_clusters_in_window = [1 2 3 4 5 6];
    Resolution = 4.3957;

elseif str2double(movie_number) == 1 && str2double(window_number) == 1

    all_clusters_in_window = [1 2 3 4];
    Resolution = 3.3453;

elseif str2double(movie_number) == 8 && str2double(window_number) == 1

    all_clusters_in_window = [1 2 3 4 5 6];
    Resolution = 3.1649;   

elseif str2double(movie_number) == 11 && str2double(window_number) == 1

    all_clusters_in_window = [1 2];
    Resolution = 3.3407;   

elseif str2double(movie_number) == 3 && str2double(window_number) == 1

    all_clusters_in_window = [1 2 3 4 5 6 7 8 9];
    Resolution = 3.3407;   
end

j1 = 1;
j2 = 1;

for w = all_clusters_in_window
    load(['MasterMatAreaABSOLUTE_Window',window_number,'Cluster',num2str(w),'.mat']);
    load(['MasterMatLengthBondABSOLUTE_Window',window_number,'Cluster',num2str(w),'.mat']);
    load(['IntensitiesPerFrameMatrix_Window',window_number,'Cluster',num2str(w),'.mat']);
   
    IntensitiesPerFrameMatrix(find(IntensitiesPerFrameMatrix==0))=NaN;

    SignallingFrame = NaN(size(IntensitiesPerFrameMatrix,2)-1,1);

    for i = 1 : size(IntensitiesPerFrameMatrix,2)-1
        if ~isempty(min(find(~isnan(IntensitiesPerFrameMatrix(:,i+1)))))
            SignallingFrame(i) = min(find(~isnan(IntensitiesPerFrameMatrix(:,i+1))));
        end
    end
   
    FirstSignallingFrame = min(SignallingFrame);

     %% find SCs area & contact length on the signalling frame %% redo removing the cells already signalling
     
     
     NBPreSignalArea(j1) = nanmean(MasterMatArea(1:FirstSignallingFrame,2));
     j1 = j1 + 1;

     for i = 1 : size(MasterMatArea,2)-2        
         OtherCellsPreSignalArea(j2) = nanmean(MasterMatArea(1:FirstSignallingFrame,i+2));
         j2 = j2 + 1;
     end

end

NBPreSignalArea = NBPreSignalArea./Resolution^2;
OtherCellsPreSignalArea = OtherCellsPreSignalArea./Resolution^2;


save(['NBPreSignalArea',movie_number,'_Window',window_number,'.mat'],'NBPreSignalArea')
save(['OtherCellsPreSignalArea',movie_number,'_Window',window_number,'.mat'],'OtherCellsPreSignalArea')