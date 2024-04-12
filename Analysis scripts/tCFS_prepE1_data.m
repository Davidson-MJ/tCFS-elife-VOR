% tCFS_prepE1_data


%this script concats group effects for b-CFS and re-CFS, and exports in
%JASP friendly format.


% reorient load.

job=[];
job.concatGFX= 1;
% job.plotPFX=1;
% job.plotGFX=0;
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir = 'C:\Users\mdav0285\Documents\GitHub\tCFS-ver1';
% homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1';

cd(homedir)
% cd(['DATA' filesep 'E1 mat files']);% 

cd(['DATA - anon' filesep 'Tracking_Exp1']);% 
datadir = pwd;
% datafols = {'DownThresholds', 'UpThresholds', 'Tracking'};
datafols = {'reCFS_discrete', 'bCFS_discrete', 'tCFS_continuous'};
%%
if job.concatGFX== 1
    %%
    GFX_table=table();
    GFX_histograms=[];
    binsAt= [0:.25:10]; % bin locations for histograms across participants.
    for id=1:3

        cd(datadir)
        cd(datafols{id});
        allfols = dir;

        %%
        pfols= striphiddenFiles(allfols);
        %%
        for ippant = 1:length(pfols)

%             cd(datadir);
%             cd(datafols{id});
%             cd(pfols(ippant).name);

%             lfile = dir([pwd filesep pfols(ippant).name '*.mat']);
%             load(lfile.name);

load(pfols(ippant).name)
            % now we have the relevant data (see Data_notes for details).

            %key variables are
            % contrFlipPts = contrast at responses (8 x 8, ims x measures) for
            % discrete. (8 x 16 for continuous).

            % allProfileCnts = frame at response

            %contrProfile = profile of contrast change within a trial.

            %convert contrFlipPts to dB
            FaceIms= 1:4;
            ObjIms= 5:8;

            contrFlipPts = 20*log10(contrFlipPts);
            [nIms,nResps] = size(contrFlipPts);

            %compare all contrFlipPts to mask contrast (0.1). (dB= -20)?


            % we can also compute te difference since last flip:
            contrFlip_diffs= zeros(size(contrFlipPts));
            for itrial= 1:size(contrFlipPts,1)
                for iblock = 1:size(contrFlipPts,3)

                    if id==1 || id==2 % show change since starting contrast.
                    tmpT = squeeze(contrFlipPts(itrial,:,iblock));
                    subD= repmat(20*log10(contrProfile{itrial, iblock}(1)), [1, length(tmpT)]);
                    contrFlip_diffs(itrial,:,iblock) = abs(tmpT-subD);
                    else

                    tmpT = [0,squeeze(contrFlipPts(itrial,:,iblock))];
                    contrFlip_diffs(itrial,:,iblock)= abs(diff(tmpT)); % have to unwrap like this to perform diff within trial.
                    end
                end
            end
            


            %category rows.
            Facerows = find(ismember(imOrder, FaceIms));
            Objrows= find(ismember(imOrder, ObjIms));


            % code now updated to store both types of analysis (average of thresholds,
            % average of change in dB since last flip).

            useData = {contrFlipPts, contrFlip_diffs};
            % include a suffix to distinguish data types
            sfx= {'', '_diff'};
            for idata= 1:2

                cfsData = useData{idata};
             

                % note that contrasts start at different directions.
                if id==1 || id==2 % downthresholds. (re-CFS), or Upthresholds (b-CFS)

                    % note that all Flip points are same direction.
                    Contr_perIm= squeeze(nanmean(cfsData,2));
                    Contr_Face = Contr_perIm(Facerows);
                    Contr_Objs = Contr_perIm(Objrows);

                    %>> store: (DOwn/Up)
                    GFX_table.([datafols{id} '_all' sfx{idata}])(ippant) = mean(Contr_perIm);
                    GFX_table.([datafols{id} '_Face' sfx{idata}])(ippant) = mean(Contr_Face);
                    GFX_table.([datafols{id} '_Obj' sfx{idata}])(ippant) = mean(Contr_Objs);

                    if id==2
                        %% include difference when data available.
                        GFX_table.(['Diff_all' sfx{idata}])(ippant) =  GFX_table.([datafols{2} '_all' sfx{idata}])(ippant) - GFX_table.([datafols{1} '_all' sfx{idata}])(ippant) ;
                        GFX_table.(['Diff_Face' sfx{idata}])(ippant) =  GFX_table.([datafols{2} '_Face' sfx{idata}])(ippant) - GFX_table.([datafols{1} '_Face' sfx{idata}])(ippant) ;
                        GFX_table.(['Diff_Obj' sfx{idata}])(ippant) =  GFX_table.([datafols{2} '_Obj'  sfx{idata}])(ippant) - GFX_table.([datafols{1} '_Obj' sfx{idata}])(ippant) ;
                    end

                elseif id==3 % tracking
                    TargVisible= 2:2:size(cfsData,2); % bCFS.
                    TargInvisible= 1:2:size(cfsData,2); % reCFS.

                    Contr_perIm_Vis= squeeze(nanmean(cfsData(:,TargVisible),2));
                    Contr_perIm_Invis= squeeze(nanmean(cfsData(:,TargInvisible),2));
                    %>> store:


                    %Upthresholds (bCFS)
                    GFX_table.(['UpThresholds_tCFS_all' sfx{idata}])(ippant) = mean(Contr_perIm_Vis);
                    GFX_table.(['UpThresholds_tCFS_Face' sfx{idata}])(ippant) = mean(Contr_perIm_Vis(Facerows));
                    GFX_table.(['UpThresholds_tCFS_Obj' sfx{idata}])(ippant) = mean(Contr_perIm_Vis(Objrows));


                    %Downthresholds (reCFS)
                    GFX_table.(['DownThresholds_tCFS_all' sfx{idata} ])(ippant) = mean(Contr_perIm_Invis);
                    GFX_table.(['DownThresholds_tCFS_Face' sfx{idata} ])(ippant) = mean(Contr_perIm_Invis(Facerows));
                    GFX_table.(['DownThresholds_tCFS_Obj' sfx{idata} ])(ippant) = mean(Contr_perIm_Invis(Objrows));


                    %adding diff suppression depth (for Bayesian stats).
                    GFX_table.(['DiffThresholds_tCFS_all' sfx{idata} ])(ippant) = mean(Contr_perIm_Vis)-  mean(Contr_perIm_Invis);
                    GFX_table.(['DiffThresholds_tCFS_Face' sfx{idata} ])(ippant) = mean(Contr_perIm_Vis(Facerows))- mean(Contr_perIm_Invis(Facerows));
                    GFX_table.(['DiffThresholds_tCFS_Obj' sfx{idata} ])(ippant) = mean(Contr_perIm_Vis(Objrows))-mean(Contr_perIm_Invis(Objrows));



                    % also include diff by flip number
                    for iflip = 1:size(cfsData,2);
                    GFX_table.(['dB_flip_' num2str(iflip) sfx{idata}])(ippant) = mean(cfsData(:,iflip),1);
                    end
                end

                %%% also store histograms.


                if id<3

                    pCounts= histcounts(allProfileCnts(:)./60, binsAt);
                    GFX_histograms.([datafols{id} '_all'])(ippant,:) = pCounts;
                    GFX_histograms.([datafols{id} '_all_mean'])(ippant) = mean(allProfileCnts(:)./60);

                    GFX_histograms.([datafols{id} '_all_median'])(ippant) = median(allProfileCnts(:)./60);
                elseif id==3

                    % we need to calculate the duration relative to previous.
                    perceptDurs = zeros(size(allProfileCnts));
                    for itrial = 1:size(allProfileCnts,1);
                        for iblock = 1:size(allProfileCnts,3);

                            trialPs = diff([1, squeeze(allProfileCnts(itrial,:,iblock))])./60;
                            perceptDurs(itrial,:,iblock)= trialPs;
                        end
                    end % itrial

                    % sort by type:
                    bCFS_durs = perceptDurs(:,TargInvisible,:);
                    reCFS_durs = perceptDurs(:,TargVisible,:);% note the index is switched, since the duration corresponds to previous percept!

                    GFX_histograms.(['reCFS_tCFS_all'])(ippant,:) = histcounts(reCFS_durs(:), binsAt);
                    GFX_histograms.(['reCFS_tCFS_all_mean'])(ippant) = mean(reCFS_durs(:));
                    GFX_histograms.(['reCFS_tCFS_all_median'])(ippant) = median(reCFS_durs(:));

                    GFX_histograms.(['bCFS_tCFS_all'])(ippant,:) = histcounts(bCFS_durs(:), binsAt);
                    GFX_histograms.(['bCFS_tCFS_all_mean'])(ippant,:) =  mean(bCFS_durs(:));
                    GFX_histograms.(['bCFS_tCFS_all_median'])(ippant,:) =  median(bCFS_durs(:));

                    % and overall:
                    GFX_histograms.(['tCFS_all'])(ippant,:) = histcounts(perceptDurs(:), binsAt);
                    GFX_histograms.(['tCFS_all_mean'])(ippant,:) = mean(perceptDurs(:));
                    GFX_histograms.(['tCFS_all_median'])(ippant,:) = median(perceptDurs(:));


                end %id
            end % idata
        end % ippant
    end % per id
    cd(datadir)

    save('GFX_table', 'GFX_table');
    save('GFX_histograms', 'GFX_histograms');
    %%
    writetable(GFX_table, 'GFX_table.csv')
end % concat GFX job.

%%