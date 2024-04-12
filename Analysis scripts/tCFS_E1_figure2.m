% tCFS_E1_figure2

% job=[];
% job.concatGFX= 1;
% job.plotPFX=0;
% job.plotGFX=0;
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir='/Users/matthewdavidson/Documents/GitHub/tCFS-ver1/'
cd(homedir)
cd(['DATA' filesep 'E1 mat files']);% filesep 'MaskTracking_v1'])
datadir = pwd;


load('GFX_table.mat');


barD = [GFX_table.DownThresholds_all, GFX_table.UpThresholds_all,...
    GFX_table.DownThresholds_tCFS_all, GFX_table.UpThresholds_tCFS_all];

%% plot pretty MS figure:

% first scatter the individual data points.

% clf;
xM = [.8, 1.2, 1.8, 2.2];
nppants= size(barD,1);
mB = mean(barD,1);
% cla
subplot(3,2,3:6);
%reCFS?
rD=plot(xM(1,[1,3]),mB(1,[1,3]), 'rd', 'MarkerSize', 10,'LineWidth',2);
hold on;
bD=plot(xM(1,[2,4]),mB(1,[2,4]), 'bd', 'MarkerSize', 10, 'LineWidth',2);
hold on;
xlim([.5 2.5])
shg
%
% add errorbars
% cla
stE= CousineauSEM(barD);
errorbar(xM, mean(barD,1), stE, 'k','LineStyle','none', 'LineWidth',3)
ylim([-34 0])
shg
%
jit_width=.03;
offset= .1;
offset= offset.* [1, -1, 1, -1];
scCols={'r', 'b', 'r','b'}
Xcollect=[];
for iX=1:4

    xAt= repmat(xM(iX), [1,nppants]);
    %random jitter per datapoint.
    
    jit= (-.5 + rand(1, nppants)) * jit_width;
    
    %scatter:
  sc=  scatter(xAt+jit+offset(iX), barD(:,iX))
sc.MarkerFaceColor= scCols{iX};
sc.MarkerEdgeColor= 'k';

Xcollect(:,iX) = xAt+jit+offset(iX);
end

% now connect the dots:
for idot= 1:nppants

plot([Xcollect(idot,1), Xcollect(idot,2)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


end
ylabel('Target contrast (dB)');
set(gca,'Xtick', [1,2], 'XTickLabel', {'discrete', 'tCFS'});
legend([bD,rD], {'bCFS', 'reCFS'},'Location', 'South', 'Orientation', 'horizontal', 'autoupdate', 'off');
set(gca,'fontsize',18);
% % axis tight

xlim([.5 2.5])
ylim([-40 0]);
shg

%% display depth:
Discrete_diff = barD(:,1) - barD(:,2);
Continuous_diff = barD(:,3 ) - barD(:,4);

disp(['Discrete difference: M=' sprintf('%.2f', mean(Discrete_diff)) ',SD=' sprintf('%.2f', std(Discrete_diff))])
disp(['Continuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])
[h,p,CI,stat] = ttest(Discrete_diff, Continuous_diff);
%%
disp(['tContinuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])


