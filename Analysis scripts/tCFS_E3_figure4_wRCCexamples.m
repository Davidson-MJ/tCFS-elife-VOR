% tCFS_E3_figure4_v2


% this version condenses the plots to include example trials:
% job=[];
% job.concatGFX= 1;
% job.plotPFX=0;
% job.plotGFX=0;


% plotType==1; % as per other MS figure 
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1';

cd(homedir)
cd(['DATA' filesep 'Tracking_Exp3']);% filesep 'MaskTracking_v1'])
datadir = pwd;

fntsize= 14;


load('GFX_table.mat');
%
% compare each type (of 5).
barD = [GFX_table.DownThresholds_tCFS_Slow, GFX_table.UpThresholds_tCFS_Slow,...
    GFX_table.DownThresholds_tCFS_Normal, GFX_table.UpThresholds_tCFS_Normal,...
    GFX_table.DownThresholds_tCFS_Fast, GFX_table.UpThresholds_tCFS_Fast];

% plot pretty MS figure:
% first plot example slow, normal, and fast trials.

pp= 'AP';
cd(pp);
lfile = dir([pwd filesep pp '*']);
load(lfile.name, 'contrProfile',...containst the raw data:
    'allProfileCnts','randDecrStepsIms'); %contains condition info.



%category rows.
rocOrder = randDecrStepsIms(1,:);

slowRows = find(ismember(rocOrder, 1));
normRows = find(ismember(rocOrder, 2));
fastRows = find(ismember(rocOrder, 3));

% plot the second example of each type:
%%
clf
cb=cbrewer('seq', 'Blues', 10);
blueCols = cb(7:10,:);
usetrials= [slowRows(1), normRows(1), fastRows(1)];
linetps= {':', '-.', '-'};
lgs= {'Slow RCC', 'Medium RCC', 'Fast RCC'};
for itrial= 1:3; 
subplot(6,1,itrial+3);
hold on;
dbData= 20*log10(contrProfile{usetrials(itrial)});
ps=plot([1:length(dbData)]./60, dbData, '.', 'color', blueCols(itrial,:), 'linew',2);
ps=plot([1:length(dbData)]./60, dbData, 'linestyle','-', 'color', blueCols(itrial,:), 'linew',2);
    xlim([1  120]);
    ylim([-30 0]);
    if  itrial<3
        set(gca,'XTickLabel', [])
    end
    
    if itrial==2
        ylabel('Target contrast (dB)');
    end

set(gca,'FontSize',fntsize)
legend(ps,lgs{itrial}, 'Location', 'northeast', 'autoupdate', 'off');

% add flip points?
flps= allProfileCnts(usetrials(itrial),:);
bcfs= 2:2:length(flps);
recfs= 1:2:length(flps);


for ifl= 1:length(bcfs);

    xAt = allProfileCnts(usetrials(itrial), bcfs(ifl));
plot(xAt./60,dbData(xAt), 'bo', 'LineWidth',1, 'MarkerFaceColor', 'b');

    xAt = allProfileCnts(usetrials(itrial), recfs(ifl));
    plot(xAt./60,dbData(xAt), 'ro' ,'MarkerFaceColor', 'r');

end
end
%%
% ylabel('Target contrast (dB)');
xlabel('Time (sec)');
set(gca,'FontSize',fntsize);
%
    % Next first scatter the individual data points.


xM= [1 1 2 2 3 3 ];
xM= xM +[ -.2 ,.2 ,-.2, .2 ,-.2,.2];

% xM = [.8, 1.2, 1.8, 2.2]; % space bad boys
nppants= size(barD,1);
mB = mean(barD,1);
% cla
subplot(2,2,1)%[1,2]);
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
ylim([-34 0])
shg
%
jit_width=.03;
offset= .1;
offset= offset.* [1, -1, 1, -1,1, -1];
scCols={'r', 'b', 'r','b','r','b','r','b','r','b'};
Xcollect=[];
for iX=1:6

    xAt= repmat(xM(iX), [1,nppants]);
    %random jitter per datapoint.
    
    jit= (-.5 + rand(1, nppants)) * jit_width;
    
    %scatter:
%   sc=  scatter(xAt+jit+offset(iX), barD(:,iX))
% sc.MarkerFaceColor= scCols{iX};
% sc.MarkerEdgeColor= 'k';

Xcollect(:,iX) = xAt+jit+offset(iX);
end

% now connect the dots:
for idot= 1:nppants

plot([Xcollect(idot,1), Xcollect(idot,2)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


plot([Xcollect(idot,5), Xcollect(idot,6)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
% plot([Xcollect(idot,7), Xcollect(idot,8)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


% plot([Xcollect(idot,9), Xcollect(idot,10)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,.4]);
% plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, .4]);


end
ylabel('Target contrast (dB)');
set(gca,'Xtick', [1:3], 'XTickLabel', {'Slow', 'Medium', 'Fast'});
xlabel('Rate of contrast change');
legend([bD,rD], {'bCFS', 'reCFS'},'Location', 'South', 'Orientation','horizontal');
set(gca,'fontsize',fntsize);
% axis tight
xlim([.5 3.5])
ylim([-40 0])
shg


% display depth:
subplot(2,2,2);

Slow_diff = barD(:,1) - barD(:,2);
Normal_diff = barD(:,3 ) - barD(:,4);
Fast_diff= barD(:,5 ) - barD(:,6);
% Phase_diff= barD(:,7 ) - barD(:,8);
% Polar_diff= barD(:,9 ) - barD(:,10);

barDiff= [Slow_diff, Normal_diff, Fast_diff];%, Phase_diff, Polar_diff];

barM= abs(mean(barDiff,1));
barErr = CousineauSEM(barDiff);

%add scatter (first)
jit_width=.06
for iX=1:3

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
bh=bar(1:3, abs(barM)); 
bh.FaceColor= [.8 .8 .8];
bh.FaceAlpha= .4;

hold on;
errorbar(1:3, barM, barErr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth',2)
%

set(gca,'Xtick', [1:3], 'XTickLabel', {'Slow ', 'Medium', 'Fast'});
xlim([.5 3.5]);
ylim([0 35]);
ylabel('Suppression depth (dB)');
set(gca,'fontsize',fntsize);
% xlabel('Rate of contrast change')
%
% legend(bh,'1','Location', 'NorthEastOutside', 'visible', 'off')
set(gcf, 'color', 'w', 'Units', 'normalized','Position', [.1 .1 .5 .8])
shg
xlabel('Rate of contrast change');

hold on;
plot([1 3], [30 30], 'k-', 'linew',1);
text(2, 32, '\itp \rm< .001', 'fontsize',12, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
shg

% 
% 
% 
% disp(['Discrete difference: M=' sprintf('%.2f', mean(Discrete_diff)) ',SD=' sprintf('%.2f', std(Discrete_diff))])
% disp(['Continuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])
% [h,p,CI,stat] = ttest(Discrete_diff, Continuous_diff);
% %%
% disp(['tContinuous difference: M=' sprintf('%.2f', mean(Continuous_diff)) ',SD=' sprintf('%.2f', std(Continuous_diff))])


