% tCFS_prepE2_data


%this script concats group effects for b-CFS and re-CFS, and exports in
%JASP friendly format.


% reorient load.

job=[];
job.concatGFX= 1;
% job.plotPFX=1;
% job.plotGFX=0;
%% set up some directories.
%PC
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir = 'C:\Users\mdav0285\Documents\GitHub\tCFS-ver1';

%macbook:
% homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1/';
cd(homedir)
%
cd(homedir)
%
cd(['DATA - anon' filesep 'Tracking_Exp2']);% filesep 'MaskTracking_v1'])
datadir = pwd;
%% datafols = {'DownThresholds', 'UpThresholds', 'Tracking'};

if job.concatGFX== 1
    %%
    GFX_table = table();
    GFX_histograms=[];
        binsAt= [0:.25:10]; % bin locations for histograms across participants.

    cd(datadir)
 
allfols = dir;
%%
 
        pfols= striphiddenFiles(allfols);
%%
    for ippant = 1:length(pfols)

    cd(datadir);
%     cd(pfols(ippant).name);
    
%     lfile = dir([pwd filesep pfols(ippant).name '*.mat']);
    load(pfols(ippant).name);

    % now we have the relevant data (see Data_notes for details). 

    %key variables are 
    % contrFlipPts = contrast at responses (10 x 20, ims x measures) for
    % discrete. (8 x 16 for continuous).
    
    % allProfileCnts = frame at response  
   
    %contrProfile = profile of contrast change within a trial.

    %convert contrFlipPts to dB
    FaceIms= 1:2;    
    ObjIms= 3:4;
    GratIms = 5:6;
    PhaseIms= 7:8;
    PolarIms = 9:10;

    contrFlipPts = 20*log10(contrFlipPts);
    [nIms,nResps] = size(contrFlipPts);



    % we can also compute te difference since last flip:
    contrFlip_diffs= zeros(size(contrFlipPts));
    contrFrame_diffs= zeros(size(allProfileCnts));
    for itrial= 1:size(contrFlipPts,1)
        for iblock = 1:size(contrFlipPts,3)

           
                tmpT = [0,squeeze(contrFlipPts(itrial,:,iblock))];
                contrFlip_diffs(itrial,:,iblock)= abs(diff(tmpT)); % have to unwrap like this to perform diff within trial.
                
                tmpT = [0,squeeze(allProfileCnts(itrial,:,iblock))];
                contrFrame_diffs(itrial,:,iblock)= diff(tmpT);
        end
    end

%category rows.
Facerows = find(ismember(imOrder, FaceIms));
Objrows= find(ismember(imOrder, ObjIms));

Gratrows = find(ismember(imOrder, GratIms));
Phaserows= find(ismember(imOrder, PhaseIms));

Polarrows = find(ismember(imOrder, PolarIms));


TargVisible= 2:2:size(contrFlipPts,2); % bCFS.
TargInvisible= 1:2:size(contrFlipPts,2); % reCFS.

% code now updated to store both types of analysis (average of thresholds,
% average of change in dB since last flip).

useData = {contrFlipPts, contrFlip_diffs, contrFrame_diffs};
% include a suffix to distinguish data types
sfx= {'', '_diff', '_rts'};
for idata= 1:3

    cfsData = useData{idata};

Contr_perIm_Vis= squeeze(mean(cfsData(:,TargVisible),2));
Contr_perIm_Invis= squeeze(mean(cfsData(:,TargInvisible),2));
%>> store:
GFX_table.(['UpThresholds_tCFS_all'  sfx{idata}])(ippant) = mean(Contr_perIm_Vis);
GFX_table.(['DownThresholds_tCFS_all'  sfx{idata}])(ippant) = mean(Contr_perIm_Invis);

%now store per type:
% we want the intersection of each type.
allImrows = {Facerows, Objrows, Gratrows, Phaserows, Polarrows};
%save as fields:
imFields={'Face', 'Obj', 'Grat', 'Phase', 'Polar'};

for iIm= 1:5
%Upthresholds (bCFS)
GFX_table.(['UpThresholds_tCFS_' imFields{iIm}  sfx{idata}])(ippant) = mean(Contr_perIm_Vis(allImrows{iIm}));
%Downthresholds (reCFS)
GFX_table.(['DownThresholds_tCFS_' imFields{iIm}  sfx{idata}])(ippant) = mean(Contr_perIm_Invis(allImrows{iIm}));

%Difference 
GFX_table.(['DiffThresholds_tCFS_' imFields{iIm}  sfx{idata}])(ippant) = mean(Contr_perIm_Invis(allImrows{iIm})) -mean(Contr_perIm_Vis(allImrows{iIm})) ;

%Contrast mean:
GFX_table.(['meanContrast_tCFS_' imFields{iIm}  sfx{idata}])(ippant) =   (mean(Contr_perIm_Invis(allImrows{iIm}))  + mean(Contr_perIm_Vis(allImrows{iIm})) )/2 ;

end




% we need to calculate the duration relative to previous.
perceptDurs = zeros(size(allProfileCnts));
for itrial = 1:size(allProfileCnts,1);
   

        trialPs = diff([1, squeeze(allProfileCnts(itrial,:))])./60;
        perceptDurs(itrial,:,iblock)= trialPs;
  
end % itrial


% summary stats:
pCounts= histcounts(perceptDurs(:), binsAt);
GFX_histograms.(['tCFS_all'])(ippant,:) = pCounts;
GFX_histograms.(['tCFS_all_mean'])(ippant) = mean(perceptDurs(:));
GFX_histograms.(['tCFS_all_median'])(ippant) = median(perceptDurs(:));


 % sort by type:
 bCFS_durs = perceptDurs(:,TargInvisible,:);
  reCFS_durs = perceptDurs(:,TargVisible,:);% note the index is switched, since the duration corresponds to previous percept!
    
   GFX_histograms.(['reCFS_tCFS_all'])(ippant,:) = histcounts(reCFS_durs(:), binsAt);
      GFX_histograms.(['reCFS_tCFS_all_mean'])(ippant) = mean(reCFS_durs(:));
            GFX_histograms.(['reCFS_tCFS_all_median'])(ippant) = median(reCFS_durs(:));


   GFX_histograms.(['bCFS_tCFS_all'])(ippant,:) = histcounts(bCFS_durs(:), binsAt);
  GFX_histograms.(['bCFS_tCFS_all_mean'])(ippant) = mean(reCFS_durs(:));
    GFX_histograms.(['bCFS_tCFS_all_median'])(ippant) = median(reCFS_durs(:));

      
end % idata


end % per participant.

cd(datadir)

save('GFX_table', 'GFX_table');
save('GFX_histograms', 'GFX_histograms');
%%
writetable(GFX_table, 'GFX_table.csv')
end % concat GFX job.

%% 