% tCFS_E2_figure3

% job=[];
% job.concatGFX= 1;
% job.plotPFX=0;
% job.plotGFX=0;


% plotType==1; % as per other MS figure 
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir = 'C:\Users\mdav0285\Documents\GitHub\tCFS-ver1';

% homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1';
cd(homedir)
cd(['DATA' filesep 'Tracking_Exp2']);% filesep 'MaskTracking_v1'])
datadir = pwd;


load('GFX_table.mat');
%%
use_dBorRT=1;
%%
% compare each type (of 5).
if use_dBorRT==1
barD = [ GFX_table.DownThresholds_tCFS_Polar, GFX_table.UpThresholds_tCFS_Polar,...
    GFX_table.DownThresholds_tCFS_Face, GFX_table.UpThresholds_tCFS_Face,...
    GFX_table.DownThresholds_tCFS_Obj, GFX_table.UpThresholds_tCFS_Obj,...
    GFX_table.DownThresholds_tCFS_Grat, GFX_table.UpThresholds_tCFS_Grat,...
        GFX_table.DownThresholds_tCFS_Phase, GFX_table.UpThresholds_tCFS_Phase];
else
    barD = [ GFX_table.DownThresholds_tCFS_Polar_rts, GFX_table.UpThresholds_tCFS_Polar_rts,...
    GFX_table.DownThresholds_tCFS_Face_rts, GFX_table.UpThresholds_tCFS_Face_rts,...
    GFX_table.DownThresholds_tCFS_Obj_rts, GFX_table.UpThresholds_tCFS_Obj_rts,...
    GFX_table.DownThresholds_tCFS_Grat_rts, GFX_table.UpThresholds_tCFS_Grat_rts,...
        GFX_table.DownThresholds_tCFS_Phase_rts, GFX_table.UpThresholds_tCFS_Phase_rts];
    
    barD= barD./60
end
%% plot pretty MS figure:

% first scatter the individual data points.

clf;
xM= [1 1 2 2 3 3 4 4 5 5];
xM= xM +[ -.2 ,.2 ,-.2, .2 ,-.2,.2,-.2,.2,-.2,.2];

% xM = [.8, 1.2, 1.8, 2.2]; % space bad boys
nppants= size(barD,1);
mB = mean(barD,1);
% cla
subplot(2,2,[1,2]);
%reCFS?
reCFSpos = 1:2:size(barD,2);
rD=plot(xM(1,[reCFSpos]),mB(1,[reCFSpos]), 'r:d', 'MarkerSize', 10,'LineWidth',2);
hold on;
%
bCFSpos = 2:2:size(barD,2);

bD=plot(xM(1,[bCFSpos]),mB(1,[bCFSpos]), 'b:d', 'MarkerSize', 10, 'LineWidth',2);
hold on;
% xlim([.5 2.5])
shg
%
% add errorbars
% cla
stE= CousineauSEM(barD);
errorbar(xM, mean(barD,1), stE, 'k','LineStyle','none', 'LineWidth',3)
if use_dBorRT==1
ylim([-34 0])
else
    axis tight
end
shg
%
jit_width=.03;
offset= .1;
offset= offset.* [1, -1, 1, -1,1, -1,1, -1,1, -1];
scCols={'r', 'b', 'r','b','r','b','r','b','r','b'};
Xcollect=[];
%%
for iX=1:10

    xAt= repmat(xM(iX), [1,nppants]);
    %random jitter per datapoint.
    
    jit= (-.5 + rand(1, nppants)) * jit_width;
    
    %scatter:
%   sc=  scatter(xAt+jit+offset(iX), barD(:,iX))
% sc.MarkerFaceColor= scCols{iX};
% sc.MarkerEdgeColor= 'k';

Xcollect(:,iX) = xAt+jit+offset(iX);
end
%%
% now connect the dots:
for idot= 1:nppants

plot([Xcollect(idot,1), Xcollect(idot,2)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


plot([Xcollect(idot,5), Xcollect(idot,6)],  [barD(idot,5), barD(idot,6)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
plot([Xcollect(idot,7), Xcollect(idot,8)],  [barD(idot,7), barD(idot,8)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


plot([Xcollect(idot,9), Xcollect(idot,10)],  [barD(idot,9), barD(idot,10)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
% plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


end
%%

set(gca,'Xtick', [1:5], 'XTickLabel', {'Polar','Faces', 'Objects', 'Gratings', 'Scrambled'});
legend([bD,rD], {'bCFS', 'reCFS'},'Location', 'South', 'Orientation', 'Horizontal');
set(gca,'fontsize',15);
% axis tight
xlim([.5 5.5])

if use_dBorRT==1
ylim([-40 0])
ylabel('Target contrast (dB)');
else
    axis tight
    ylabel('mean duration (sec)');
end

shg
%%

% display depth:
subplot(2,2,3:4);

Face_diff = barD(:,1) - barD(:,2);
Obj_diff = barD(:,3 ) - barD(:,4);
Grat_diff= barD(:,5 ) - barD(:,6);
Phase_diff= barD(:,7 ) - barD(:,8);
Polar_diff= barD(:,9 ) - barD(:,10);

barDiff= [Face_diff, Obj_diff, Grat_diff, Phase_diff, Polar_diff];

barM= abs(mean(barDiff,1))
barErr = CousineauSEM(barDiff);

%add scatter (first)
jit_width=.06
for iX=1:5

    xAt= repmat(iX, [1,nppants]);
    %random jitter per datapoint.
    
    jit= (-.5 + rand(1, nppants)) * jit_width;
    
    %scatter:
  sc=  scatter(xAt+jit, abs(barDiff(:,iX))); 
  hold on
sc.MarkerFaceColor= [.8 .8 .8];
sc.MarkerEdgeColor= 'w';

% Xcollect(:,iX) = xAt+jit+offset(iX);
end
hold on;
bh=bar(1:5, abs(barM)); 
bh.FaceColor= [.8 .8 .8];
bh.FaceAlpha= .4;

hold on;
errorbar(1:5, barM, barErr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth',2)
%


set(gca,'Xtick', [1:5], 'XTickLabel', {, 'Polar','Faces', 'Objects', 'Gratings', 'Scrambled'});
xlim([.5 5.5]);

if use_dBorRT==1
ylim([0 40])
ylabel('Target contrast (dB)');
else
    ylim([-1 1])
    ylabel('mean duration (sec)');
end

shg

set(gca,'fontsize',15);
set(gcf, 'color', 'w', 'Units', 'normalized','Position', [.1 .1 .4 .7])
%%
hold on;
plot([1 5], [30 30], 'k-', 'linew',1);
text(3, 32, '\itns', 'fontsize',12, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
shg

% 
% 
% 
% disp(['Discrete difference: M=' sprintf('%.2f', mean(Discrete_diff)) ',SD=' sprintf('%.2f', std(Discrete_diff))])
% disp(['Continuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])
% [h,p,CI,stat] = ttest(Discrete_diff, Continuous_diff);
% %%
% disp(['tContinuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])


