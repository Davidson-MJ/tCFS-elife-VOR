% plot histograms

job.plot_GFXhistograms=1; % MS ver

job.plot_GFXhistograms_cdf=1; % review ver
%% set up some directories.
%PC
homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
%macbook:
% homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1/';
cd(homedir)
folsare = {'E1 mat files', 'Tracking_Exp2', 'Tracking_Exp3'};

%including option to plot cdf instead:

clf;
%%


binsAt= [0:.25:10]; % bin locations for histograms across participants.

cmap= cbrewer('seq', 'RdPu', 100);

speedCols= {cmap(33,:), cmap(66,:), cmap(99,:)};

if job.plot_GFXhistograms==1
%
for iExp =1:3
    cd(homedir);
    cd(['DATA' filesep folsare{iExp}]);% filesep 'MaskTracking_v1'])

    load('GFX_histograms.mat');

%     subplot(2,2,iExp); hold on

    subplot(1,2,1); hold on
    cla

xlabel('Percept duration (sec)')
ylabel('Frequency')
xlim([0 10]);
ylim([0 17]); hold on;
set(gca,'fontsize', 12, 'Xtick', [0:2:max(binsAt)])

title(['Experiment ' num2str(iExp)])
colsAre = {'b', 'r', 'k'};

if iExp==1

    fldsAre= {'UpThresholds_all', 'DownThresholds_all', 'bCFS_tCFS_all', 'reCFS_tCFS_all', 'tCFS_all'};
    colsAre= {'b', 'r', 'b', 'r', 'k'};
    lstls={':', ':', '-', '-', '-'};
else
    fldsAre= {'reCFS_tCFS_all', 'bCFS_tCFS_all', 'tCFS_all'};
 lstls={'-', '-', '-'};

end
lgh=[]
for ifield =1 :length(fldsAre);
    tmpData = GFX_histograms.(fldsAre{ifield});
    mData = squeeze(mean(tmpData,1));
    errData = CousineauSEM(tmpData);

    sh=shadedErrorBar(binsAt(2:end), mData, errData, {'color', colsAre{ifield}}, 1)
    sh.mainLine.LineWidth= 2;
sh.mainLine.LineStyle= lstls{ifield};
    lgh(ifield) = sh.mainLine;

     % add mean plot:
    mAt = mean(GFX_histograms.([fldsAre{ifield} '_median']));
%     xAt =dsearchn(binsAt', mAt);
%     plot([mAt, mAt], [0 mData(xAt-1)], 'color', colsAre{ifield}, 'LineWidth', 1);
%     plot(mAt, 0.2, 'v', 'color', colsAre{ifield}, 'MarkerFaceColor',colsAre{ifield} )

end
        if iExp==1
legend([lgh(1), lgh(2), lgh(3), lgh(4), lgh(5)],{'bCFS (discrete)', 'reCFS (discrete)', 'bCFS (tCFS)', 'reCFS(tCFS)', 'all tCFS'});
        else

legend([lgh(1), lgh(2), lgh(3)],{'bCFS (tCFS)', 'reCFS (tCFS)', 'all tCFS'})
        end


    if iExp==3; % add final plot by speed        
%        subplot(2,2,4);
subplot(1,2,1); cla
xlabel('Percept duration (sec)')
ylabel('Frequency')
xlim([0 10]);
ylim([0 17]); hold on;
set(gca,'fontsize', 12)
title('Experiment 3 by RCC')
fldsAre = {'tCFS_Slow', 'tCFS_Normal', 'tCFS_Fast'};
colsAre ={cmap(33,:); cmap(66,:); cmap(99,:)};
lgh=[]
condMeans=[];
for ifield =1 :3
    tmpData = GFX_histograms.(fldsAre{ifield});
    mData = squeeze(mean(tmpData,1));
    errData = CousineauSEM(tmpData);

    sh=shadedErrorBar(binsAt(2:end), mData, errData, {'color', colsAre{ifield}}, 1)
    sh.mainLine.LineWidth= 2;
    lgh(ifield) = sh.mainLine;

    % add mean plot:
    mAt = mean(GFX_histograms.([fldsAre{ifield} '_median']));
    xAt =dsearchn(binsAt', mAt);
    plot([mAt, mAt], [0 mData(xAt-1)], 'color', colsAre{ifield}, 'LineWidth', 2);
    plot(mAt, 0.2, 'v', 'color', colsAre{ifield}, 'MarkerFaceColor',colsAre{ifield} )
condMeans(ifield) = mAt;
end
axis tight
legend([lgh(1), lgh(2), lgh(3)],{'slow', 'medium', 'fast'}, 'autoupdate', 'off');
    end
        %% export stats:
        set(gcf, 'color', 'w')
% xlabel('Duration (sec)')
end
end % job
%%

if job.plot_GFXhistograms_cdf==1
    %%
%     clf
for iExp =3
    cd(homedir);
    cd(['DATA' filesep folsare{iExp}]);% filesep 'MaskTracking_v1'])

    load('GFX_histograms.mat', 'GFX_eCDF');

% set(gca,'fontsize', 20)
subplot(122); cla
title(['Experiment ' num2str(iExp)])
colsAre = {'b', 'r', 'k'};
set(gca,'fontsize',12)
    
xlabel('Percept duration (sec)')
ylabel('F(x)')

title({['Cumulative density function by rate of contrast change']})
fldsAre = {'tCFS_Slow', 'tCFS_Normal', 'tCFS_Fast'};
colsAre ={cmap(33,:); cmap(66,:); cmap(99,:)};
lgh=[]
%
for ifield =1:3
    % overlay each ppants data:
    hold on
    
    xvec = GFX_eCDF.([fldsAre{ifield} '_ecdf_x'])(1,:);        
    tmpD= GFX_eCDF.([fldsAre{ifield} '_ecdf_y']);
      
    stE = CousineauSEM(tmpD);
    meanCDF= squeeze(mean(tmpD,1));
    
   sh= shadedErrorBar(xvec, meanCDF, stE, {'linew', 3,'color', [colsAre{ifield}]},1);
    lgh(ifield)= sh.mainLine;
    xlim([0 10])
    ylim([0 1])

    % add mean:
      % add mean plot:
    mAt = condMeans(ifield); % calculated above
    xAt =dsearchn(xvec', mAt);
    plot([mAt, mAt], [0 meanCDF(xAt-1)], 'color', colsAre{ifield}, 'LineWidth', 2);
    plot(mAt, 0, 'v', 'color', colsAre{ifield}, 'MarkerFaceColor',colsAre{ifield} )

end
% axis tight
legend([lgh(1), lgh(2), lgh(3)],{'slow', 'medium', 'fast'}, 'autoupdate', 'off', 'location', 'SouthEast', 'fontsize',20);
    
        %% export stats:
        set(gcf, 'color', 'w')
% xlabel('Duration (sec)')
end
end
    