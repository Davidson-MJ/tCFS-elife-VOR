% tCFS_prepE3_data


%this script concats group effects for b-CFS and re-CFS, and exports in
%JASP friendly format.


% reorient load.

job=[];
job.concatGFX= 1;
% job.plotPFX=1;
% job.plotGFX=0;

% adding data for average contrProfile per speed (review figure)
%% set up some directories.
%PC
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir = 'C:\Users\mdav0285\Documents\GitHub\tCFS-ver1';

%macbook:
% homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1/';
cd(homedir)
%
cd(['DATA - anon' filesep 'Tracking_Exp3']);% filesep 'MaskTracking_v1'])
datadir = pwd;
%% datafols = {'DownThresholds', 'UpThresholds', 'Tracking'};

if job.concatGFX== 1
    %%
    GFX_table = table();
    GFX_histograms=[];
    GFX_eCDF=[];
    GFX_ramps_x_time=[];
        binsAt= [0:.25:10]; % bin locations for histograms across participants.

    cd(datadir)
 
allfols = dir;
%%

        pfols= striphiddenFiles(allfols);
%%
    for ippant = 1:length(pfols)

%     cd(datadir);
%     cd(pfols(ippant).name);
    
%     lfile = dir([pwd filesep pfols(ippant).name '*.mat']);
%     load(lfile.name);
load(pfols(ippant).name)

    % now we have the relevant data (see Data_notes for details). 

    %key variables are 
    % contrFlipPts = contrast at responses (12 x 20, ims x measures) for
    
    % allProfileCnts = frame at response  
 
    %contrProfile = profile of contrast change within a trial.

    % randDecrStepsIms (2x 12). Top row is which step rate (1-3). bottom
    % row is which image.

    %convert contrFlipPts to dB
    FaceIms= 1;    
    ObjIms= 2;
    GratIms = 3;
    PhaseIms= 4;

    contrFlipPts = 20*log10(contrFlipPts);
    [nIms,nResps] = size(contrFlipPts);

% note that we also have the contrast profile, which can be displayed in dB
[contrProfile_stacked, contrProfile_stacked_diff]= deal(nan(nIms, 5000)); % long length to catch longer trials.
% we can also compute te difference since last flip:
contrFlip_diffs= zeros(size(contrFlipPts));
for itrial= 1:size(contrFlipPts,1)

    tmpT = [0,squeeze(contrFlipPts(itrial,:))];
    contrFlip_diffs(itrial,:)= abs(diff(tmpT)); % have to unwrap like this to perform diff within trial.
% tmpTrack = contrProfile{itrial};
%     contrProfile_stacked(itrial,1:length(tmpTrack)) = 20*log10(tmpTrack); % dB over time. 

%     contrProfile_stacked_diff(itrial,1:length(tmpTrack)) = diff([0, 20*log10(tmpTrack)]); % dB over time. 

end

%category rows.
rocOrder = randDecrStepsIms(1,:);
imOrder = randDecrStepsIms(2,:);

Facerows = find(ismember(imOrder, FaceIms));
Objrows= find(ismember(imOrder, ObjIms));

Gratrows = find(ismember(imOrder, GratIms));
Phaserows= find(ismember(imOrder, PhaseIms));

slowRows = find(ismember(rocOrder, 1));
normRows = find(ismember(rocOrder, 2));
fastRows = find(ismember(rocOrder, 3));



% we want the intersection of each type.
allImrows = {Facerows, Objrows, Gratrows, Phaserows};
allSpeedrows= {slowRows, normRows, fastRows};
%save as fields:
imFields={'Face', 'Obj', 'Grat', 'Phase'};

speedFields={'Slow', 'Normal', 'Fast'};


TargVisible= 2:2:size(contrFlipPts,2); % bCFS.
TargInvisible= 1:2:size(contrFlipPts,2); % reCFS.

% store large conditions:

useData= {contrFlipPts, contrFlip_diffs};
  sfx= {'', '_diff'};
for idata=1:2

    cfsdata= useData{idata};

Contr_perIm_Vis= squeeze(mean(cfsdata(:,TargVisible),2));
Contr_perIm_Invis= squeeze(mean(cfsdata(:,TargInvisible),2));
%>> store:

%bCFS, reCFS
GFX_table.(['UpThresholds_tCFS_all' sfx{idata} ])(ippant) = mean(Contr_perIm_Vis);
GFX_table.(['DownThresholds_tCFS_all' sfx{idata}])(ippant) = mean(Contr_perIm_Invis);

%speed types:
GFX_table.(['UpThresholds_tCFS_Slow' sfx{idata}])(ippant) = mean(Contr_perIm_Vis(slowRows));
GFX_table.(['UpThresholds_tCFS_Normal' sfx{idata}])(ippant) = mean(Contr_perIm_Vis(normRows));
GFX_table.(['UpThresholds_tCFS_Fast' sfx{idata}])(ippant) = mean(Contr_perIm_Vis(fastRows));

GFX_table.(['DownThresholds_tCFS_Slow' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(slowRows));
GFX_table.(['DownThresholds_tCFS_Normal' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(normRows));
GFX_table.(['DownThresholds_tCFS_Fast' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(fastRows));

GFX_table.(['DiffThresholds_tCFS_Slow' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(slowRows)) - mean(Contr_perIm_Vis(slowRows));
GFX_table.(['DiffThresholds_tCFS_Normal' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(normRows))- mean(Contr_perIm_Vis(normRows));
GFX_table.(['DiffThresholds_tCFS_Fast' sfx{idata}])(ippant) = mean(Contr_perIm_Invis(fastRows))-mean(Contr_perIm_Vis(fastRows));





% also include diff by flip number
for iflip = 1:size(cfsdata,2)
    GFX_table.(['Slow_dB_flip_' num2str(iflip) sfx{idata}])(ippant) = mean(cfsdata(slowRows,iflip),1);
        GFX_table.(['Medium_dB_flip_' num2str(iflip) sfx{idata}])(ippant) = mean(cfsdata(normRows,iflip),1);
        GFX_table.(['Fast_dB_flip_' num2str(iflip) sfx{idata}])(ippant) = mean(cfsdata(fastRows,iflip),1);

end


% now their intersection:
for imtype= 1:4
    for ispd= 1:3

        subrows = intersect(allImrows{imtype}, allSpeedrows{ispd});

GFX_table.(['UpThresholds_tCFS_' imFields{imtype} '_' speedFields{ispd} sfx{idata}])(ippant) = mean(Contr_perIm_Vis(subrows));
GFX_table.(['DownThresholds_tCFS_' imFields{imtype} '_' speedFields{ispd} sfx{idata}])(ippant) = mean(Contr_perIm_Invis(subrows));
GFX_table.(['DiffThresholds_tCFS_' imFields{imtype} '_' speedFields{ispd} sfx{idata}])(ippant) = mean(Contr_perIm_Invis(subrows)) - mean(Contr_perIm_Vis(subrows));

    end % speeds
end % imtypes


% also include the av per trialID



% we need to calculate the duration relative to previous.
perceptDurs = zeros(size(allProfileCnts));
for itrial = 1:size(allProfileCnts,1);
  

        trialPs = diff([1, squeeze(allProfileCnts(itrial,:))])./60;
        perceptDurs(itrial,:)= trialPs;
    
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
   GFX_histograms.(['bCFS_tCFS_all'])(ippant,:) = histcounts(bCFS_durs(:), binsAt);
    

   GFX_histograms.(['reCFS_tCFS_all_mean'])(ippant,:) = mean(reCFS_durs(:));
   GFX_histograms.(['bCFS_tCFS_all_mean'])(ippant,:) = mean(bCFS_durs(:));


   GFX_histograms.(['reCFS_tCFS_all_median'])(ippant,:) = median(reCFS_durs(:));
   GFX_histograms.(['bCFS_tCFS_all_median'])(ippant,:) = median(bCFS_durs(:));
  
ecdfbins = [0:.01:60];
   for ispd=1:3
    
%        pdurs= perceptDurs(allSpeedrows{ispd},:);
       GFX_histograms.(['tCFS_' speedFields{ispd}])(ippant,:) = histcounts(perceptDurs(allSpeedrows{ispd},:), binsAt);
       
       tmp= perceptDurs(allSpeedrows{ispd},:);
       GFX_histograms.(['tCFS_' speedFields{ispd} '_mean'])(ippant,:) = mean(tmp(:));
       GFX_histograms.(['tCFS_' speedFields{ispd} '_median'])(ippant,:) = median(tmp(:));

        % for reviewers comments, also calculate participant level ecdf:
        %round first to make binning easier.
        tmp = round(tmp(:),1);
        
        % manually create ecdf at each ecdfbin location.
        ecdfcount=[];
        for ibin = 1:length(ecdfbins)
            
            ecdfcount(ibin) = length(find(tmp==ecdfbins(ibin)));
        end
        Y = cumsum(ecdfcount);
        Y = Y./ max(Y); % ecdf.
        
%         [f,x]= ecdf(tmp(:)); % bugs
        GFX_eCDF.(['tCFS_' speedFields{ispd} '_ecdf_x'])(ippant,:) =ecdfbins;
        GFX_eCDF.(['tCFS_' speedFields{ispd} '_ecdf_y'])(ippant,:) =Y;
   end
   
end % idata

%% quick store:
% > store the average ramps over time 
for ispd=1:3
userows= allSpeedrows{ispd};

% need to unpack from cell array.
tmp_ts= nan(length(userows), 60*60); % first 30 seconds should be enough?
for irow =1:length(userows);
    tmpD= contrProfile{userows(irow)};
    if length(tmpD) <= length(tmp_ts)
    tmp_ts(irow,1:length(tmpD)) = 20*log10(tmpD);
    else
    tmp_ts(irow,:) = 20*log10(tmpD(1:length(tmp_ts)));
    end    
end
GFX_ramps_x_time.(['contrRamps_x_time_' speedFields{ispd}])(ippant,:) = nanmean(tmp_ts,1);

end


end % per participant.



%%

cd(datadir)

save('GFX_table', 'GFX_table');
save('GFX_histograms', 'GFX_histograms','GFX_eCDF');
%%
writetable(GFX_table, 'GFX_table.csv')

HistTable=table();
HistTable.SlowDurs_mean = GFX_histograms.tCFS_Slow_mean;

HistTable.NormalDurs_mean = GFX_histograms.tCFS_Normal_mean;

HistTable.FastDurs_mean = GFX_histograms.tCFS_Fast_mean;


HistTable.SlowDurs_median = GFX_histograms.tCFS_Slow_median;

HistTable.NormalDurs_median = GFX_histograms.tCFS_Normal_median;

HistTable.FastDurs_median = GFX_histograms.tCFS_Fast_median;
writetable(HistTable, 'HistTable.csv')


end % concat GFX job.




%% 