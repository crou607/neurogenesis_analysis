%Example script that extracts areas for fold change analysis shown in
%Figure 5

% Inputs are MasterMatAreaABSOLUTE, MasterMatLengthBondABSOLUTE and 
% IntensitiesPerFrameMatrix .mats obtained from GetSizeAndBondTables.m

close all
clear all

movie_number = '3';
window_number = '1';
window_size = 10;

if str2double(movie_number) == 9 && str2double(window_number) == 1

    isotropic_no_clusters_in_window = [1 2 3];
    ClusterStartFrame = [1 1 14];
    Resolution = 4.3957;

elseif str2double(movie_number) == 9 && str2double(window_number) == 2
    
    isotropic_no_clusters_in_window = [1 2 3 4 5 6];
    ClusterStartFrame = [10 10 20 20 20 41];
    Resolution = 4.3957;

elseif str2double(movie_number) == 1 && str2double(window_number) == 1

    isotropic_no_clusters_in_window = [1 2 3 4];
    ClusterStartFrame = [12 14 13 16];
    Resolution = 3.3453;

elseif str2double(movie_number) == 8 && str2double(window_number) == 1

    isotropic_no_clusters_in_window = [1 2 3 4 5 6];
    ClusterStartFrame = [25 23 16 21 4 20];
    Resolution = 3.1649;   

elseif str2double(movie_number) == 11 && str2double(window_number) == 1

    isotropic_no_clusters_in_window = [1 2];
    ClusterStartFrame = [44 24];
    Resolution = 3.3407;   


elseif str2double(movie_number) == 3 && str2double(window_number) == 1

isotropic_no_clusters_in_window = [1 2 3 4 5 6 7 8 9];
ClusterStartFrame = [15 20 24 20 22 21 11 26 9];
Resolution = 3.3407;   

end

j1 = 1;
j2 = 1;
j3 = 1;

j4 = 1;
j5 = 1;
j6 = 1;

for w = isotropic_no_clusters_in_window
    load(['MasterMatAreaABSOLUTE_Window',window_number,'Cluster',num2str(w),'.mat']);
    load(['MasterMatLengthBondABSOLUTE_Window',window_number,'Cluster',num2str(w),'.mat']);
    load(['IntensitiesPerFrameMatrix_Window',window_number,'Cluster',num2str(w),'.mat']);
   

    IntensitiesPerFrameMatrix(find(IntensitiesPerFrameMatrix==0)) = NaN;

    SignallingFrame = NaN(size(IntensitiesPerFrameMatrix,2)-1,1);

    for i = 1 : size(IntensitiesPerFrameMatrix,2)-1
        if ~isempty(min(find(~isnan(IntensitiesPerFrameMatrix(:,i+1)))))
            SignallingFrame(i) = min(find(~isnan(IntensitiesPerFrameMatrix(:,i+1))));
        end
    end

    SignallingCells = find(~isnan(SignallingFrame));
    NonSignallingCells = find(isnan(SignallingFrame));
    FirstSignallingFrame = min(SignallingFrame);


   for i = 1 : size(SignallingCells,1)        
       frame_to_analyse = SignallingFrame(SignallingCells(i));
        
       % if frame_to_analyse > 1
          AreaOnSignalingFrameWindow1SIGNALLING(j1) = nanmean(MasterMatArea(max(frame_to_analyse-4,1):frame_to_analyse,SignallingCells(i)+2));
          ContactOnSignalingFrameWindow1SIGNALLING(j1) = nanmean(MasterMatLengthBond(max(frame_to_analyse-4,1):frame_to_analyse,SignallingCells(i)+1));
          j1 = j1 + 1;
          
           if isempty(find(isnan(IntensitiesPerFrameMatrix(frame_to_analyse+1:frame_to_analyse+window_size,SignallingCells(i)+1))))
                AreaOnSignalingFrameWindow2SIGNALLING(j2) = nanmean(MasterMatArea(frame_to_analyse+1:frame_to_analyse+window_size,SignallingCells(i)+2));
                ContactOnSignalingFrameWindow2SIGNALLING(j2) = nanmean(MasterMatLengthBond(frame_to_analyse+1:frame_to_analyse+window_size,SignallingCells(i)+1));
                j2 = j2 + 1;
                
                if isempty(find(isnan(IntensitiesPerFrameMatrix(frame_to_analyse+window_size+1:min(size(IntensitiesPerFrameMatrix,1),frame_to_analyse+2*window_size+1),SignallingCells(i)+1))))
                   AreaOnSignalingFrameWindow3SIGNALLING(j3) = nanmean(MasterMatArea(frame_to_analyse+window_size+1:min(size(IntensitiesPerFrameMatrix,1),frame_to_analyse+2*window_size+1),SignallingCells(i)+2));
                   ContactOnSignalingFrameWindow3SIGNALLING(j3) = nanmean(MasterMatLengthBond(frame_to_analyse+window_size+1:min(size(IntensitiesPerFrameMatrix,1),frame_to_analyse+2*window_size+1),SignallingCells(i)+1));
                   j3 = j3 + 1;
                
                else
                    AreaOnSignalingFrameWindow3SIGNALLING(j3) = nanmean(MasterMatArea(frame_to_analyse+window_size+1:max(find(~isnan(IntensitiesPerFrameMatrix(:,SignallingCells(i)+1)))),SignallingCells(i)+2));
                    ContactOnSignalingFrameWindow3SIGNALLING(j3) = nanmean(MasterMatLengthBond(frame_to_analyse+window_size+1:max(find(~isnan(IntensitiesPerFrameMatrix(:,SignallingCells(i)+1)))),SignallingCells(i)+1));
                    j3 = j3 +1;
                end
          
           else 
                AreaOnSignalingFrameWindow2SIGNALLING(j2) = nanmean(MasterMatArea(frame_to_analyse+1:max(find(~isnan(IntensitiesPerFrameMatrix(:,SignallingCells(i)+1)))),SignallingCells(i)+2));
                ContactOnSignalingFrameWindow2SIGNALLING(j2) = nanmean(MasterMatLengthBond(frame_to_analyse+1:max(find(~isnan(IntensitiesPerFrameMatrix(:,SignallingCells(i)+1)))),SignallingCells(i)+1));                
                j2 = j2 + 1;
    
                AreaOnSignalingFrameWindow3SIGNALLING(j3) = NaN;
                ContactOnSignalingFrameWindow3SIGNALLING(j3) = NaN;                
                j3 = j3 + 1;
           end
    
       % end    

   end

   for i = 1 : size(NonSignallingCells,1)        
        
          AreaOnSignalingFrameWindow1NONSIGNALLING(j4) = nanmean(MasterMatArea(max(FirstSignallingFrame-4,1):FirstSignallingFrame,NonSignallingCells(i)+2));
          ContactOnSignalingFrameWindow1NONSIGNALLING(j4) = nanmean(MasterMatLengthBond(max(FirstSignallingFrame-4,1):FirstSignallingFrame,NonSignallingCells(i)+1));         
          j4 = j4 + 1;
          
          AreaOnSignalingFrameWindow2NONSIGNALLING(j5) = nanmean(MasterMatArea(FirstSignallingFrame+1:FirstSignallingFrame+window_size,NonSignallingCells(i)+2));
          ContactOnSignalingFrameWindow2NONSIGNALLING(j5) = nanmean(MasterMatLengthBond(FirstSignallingFrame+1:FirstSignallingFrame+window_size,NonSignallingCells(i)+1));         
          j5 = j5 + 1;

          AreaOnSignalingFrameWindow3NONSIGNALLING(j6) = nanmean(MasterMatArea(FirstSignallingFrame+window_size+1:min(size(IntensitiesPerFrameMatrix,1),FirstSignallingFrame+2*window_size+1),NonSignallingCells(i)+2));
          ContactOnSignalingFrameWindow3NONSIGNALLING(j6) = nanmean(MasterMatLengthBond(FirstSignallingFrame+window_size+1:min(size(IntensitiesPerFrameMatrix,1),FirstSignallingFrame+2*window_size+1),NonSignallingCells(i)+1));
          j6 = j6 + 1;

   end
    
end

AreasInWindowsSignalling = [AreaOnSignalingFrameWindow1SIGNALLING' AreaOnSignalingFrameWindow2SIGNALLING' AreaOnSignalingFrameWindow3SIGNALLING'];
AreasInWindowsNonSignalling = [AreaOnSignalingFrameWindow1NONSIGNALLING' AreaOnSignalingFrameWindow2NONSIGNALLING' AreaOnSignalingFrameWindow3NONSIGNALLING'];

AreasInWindowsSignalling = AreasInWindowsSignalling/Resolution^2;
AreasInWindowsNonSignalling = AreasInWindowsNonSignalling/Resolution^2;

ContactsInWindowsSignalling = [ContactOnSignalingFrameWindow1SIGNALLING' ContactOnSignalingFrameWindow2SIGNALLING' ContactOnSignalingFrameWindow3SIGNALLING'];
ContactsInWindowsNonSignalling = [ContactOnSignalingFrameWindow1NONSIGNALLING' ContactOnSignalingFrameWindow2NONSIGNALLING' ContactOnSignalingFrameWindow3NONSIGNALLING'];

ContactsInWindowsSignalling = ContactsInWindowsSignalling/Resolution;
ContactsInWindowsNonSignalling = ContactsInWindowsNonSignalling/Resolution;



save(['AreasInWindowsSignallingNanMean',movie_number,'_Window',window_number,'.mat'],'AreasInWindowsSignalling')
save(['AreasInWindowsNonSignallingNanMean',movie_number,'_Window',window_number,'.mat'],'AreasInWindowsNonSignalling')


save(['ContactsInWindowsSignallingNanMean',movie_number,'_Window',window_number,'.mat'],'ContactsInWindowsSignalling')
save(['ContactsInWindowsNonSignallingNanMean',movie_number,'_Window',window_number,'.mat'],'ContactsInWindowsNonSignalling')


