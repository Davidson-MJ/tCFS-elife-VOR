% tCFS fig 1

% mock data:
ContrRamp = 0:.01:1;
timeRamp = 0:1:10;
FData = 35;
nFData = 55;

FData = 65;
nFData = 85;

barData = [ContrRamp(1,[FData, nFData])];
% 
% subplot(2,4,1);
% figure(1); clf;
% set(gcf,'color', 'w', 'units', 'normalized','Position', [.2 .2 .8 .4]);
% % subplot(1,4,1)
% b = bar3(barData, .55);
% 
% % ylim([.2 1.2])
cmap =cbrewer('seq', 'Greys', 100);
% colormap(cmap);
% b.FaceColor = "interp";
% 
% for a = 1:length(b)
%     b(a).CData = b(a).ZData;
% end
% caxis([0 1])
% view(-90,0)
% grid off
% 
% zlim([0 1]); %ylim actually.
% zlabel('contrast');
% set(gca,'fontsize', 12, 'YtickLabel', []);

%%
figure(2); clf; 
subplot(1,4,1);

b = bar(barData, .55);
hold on; errorbar(1:2, barData, [.02 .03], 'k','LineStyle','none', 'LineWidth',1)
b.FaceColor= 'k';
b.BarWidth= .4;
b.FaceAlpha= .2;
b.EdgeColor= 'w';
ylim([0 1]);
% ylabel('RT (or Contrast)');
set(gca,'ytick',[0.1,.15 .85,.9], 'YTickLabel', [], 'ytick', [])
set(gcf,'color', 'w', 'units', 'normalized','Position', [.1 .1 .9 .4])
set(gca,'xticklabels', {'Face', 'non-Face'});
% ylabel('reaction times', 'VerticalAlignment','top');
% set(gca, 'YTick', [.2 .9]);
xlim([.5 2.5])
box off
colormap(cmap)
c= colorbar;
set(c, 'Location', 'Westoutside', 'ytick',[]);
% c.Ticks = [.1, .9];
% c.TickLabels = {'(Invisible)', '(Visible)'};

% set(c, 'YTickLabelRotation', 0)

ylabel(c,'contrast ramp', 'fontsize', 15 , 'VerticalAlignment', 'bottom');
set(gca,'fontsize', 15)

%
hold on
orngCol= [1, .67, .25];
plot([1, 1], [.1,barData(1)/2], ':','linew', 2,'color', orngCol)
plot([1, 1], [barData(1)/2+.1 barData(1)-.05], ':','linew', 2,'color', orngCol)

plot([2, 2], [.1,barData(2)/2], ':','linew', 2,'color', orngCol)
plot([2, 2], [barData(2)/2+.1 barData(2)-.05], ':','linew', 2,'color', orngCol)

text(1, barData(1)/2+.01, '\delta f', 'fontsize', 15, 'FontWeight', 'normal',...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom');
text(2, barData(2)/2+.01, '\delta nf', 'fontsize', 15, 'FontWeight', 'normal',...
    'HorizontalAlignment', 'center', 'VerticalAlignment','bottom');




text(1, barData(1)+.05, 'Visible', 'Color', 'b', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold')
text(2, barData(2)+.05, 'Visible', 'Color', 'b', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold')

text(1, 0, 'Invisible', 'Color', 'r', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold','BackgroundColor', [1 1 1 .5], 'EdgeColor', 'w')

tt=text(2, 0, 'Invisible', 'Color', 'r', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold','BackgroundColor', [1 1 1 .5], 'EdgeColor', 'w')

% now give hypothetical examples.
oppD= [-.8 -.8;... % equal disappearance points (A< B)
    -.9 -.45;... % 
    -.7 -.5];   % shifted (A= B)

% no diff
%
%
cmapD= flipud(cmap);
% colormap(cmapD)
% clf
for isub= 2:4
subplot(1,4,isub);
cla
b = bar([.8, 1.8],barData, .55);
b.BarWidth= .4;
b.EdgeColor='w';
% hold on; errorbar([.8 1.8], barData, [.03 .05], 'k','LineStyle','none', 'LineWidth',1)
box off

b.FaceAlpha= .05; 
b.FaceColor= 'b';
ylim([0 1]);
set(gca,'ytick', [])
set(gca, 'xtick', [1,2],'xticklabels', {'Face', 'non-Face'});

% c= colorbar;
% set(c, 'Location', 'Westoutside', 'ytick',[]);
xlim([.5 2.5])
box off
% ylabel(c,'contrast ramp', 'fontsize', 15 );

% %add text.
% text(.8, barData(1), 'Vis.', 'Color', 'b', 'HorizontalAlignment', 'center',...
%     'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','normal');
% 
% text(1.75, barData(2), 'Vis.', 'Color', 'b', 'HorizontalAlignment', 'center',...
%     'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','normal')
hold on;
% add supp depth:

xb= [.8, 1.8];
stngs = {'\delta f', '\delta nf'};
for ib= 1:2
strt = 1+oppD(isub-1,ib);
stp= barData(ib);
d1=diff([strt,stp]);

plot([xb(ib), xb(ib)], [strt,(strt+d1/2 - .05)], ':','linew', 2,'color', orngCol)
plot([xb(ib), xb(ib)], [(strt+d1/2 +.05), stp], ':','linew', 2,'color', orngCol)

% add connection"

plot([xb(ib)-.1, xb(ib)+.1], [stp,stp], '-','linew', 2,'color', 'b')
plot([xb(ib)-.1, xb(ib)+.1], [strt,strt], '-','linew', 2,'color', 'r')

text(xb(ib),strt+d1/2-.03, stngs{ib}, 'Color', 'k', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold')

end

% plot([1.75, 1.75], [1+oppD(isub-1,2),barData(2)], ':','linew', 2,'color', orngCol)

%%
% plot([.55, .95], [1+oppD(isub-1,1),1+oppD(isub-1,1)], '-','linew', 1,'color', 'r')

% plot([1.55, 1.95], [1+oppD(isub-1,2),1+oppD(isub-1,2)], '-','linew', 1,'color', 'r')



hold on;
yyaxis right;
bI= bar([1.2, 2.2],oppD(isub-1,:), .55); ylim([-1 0]);
bI.BarWidth= .4;
bI.FaceAlpha= .2;
bI.EdgeColor='w';
hold on; errorbar([1.2,2.2], oppD(isub-1,:), [.01 .02], 'k','LineStyle','none', 'LineWidth',1)
box off
bI.FaceColor= 'r';
set(gca,'YTick', [],'YColor', 'k');
box off
xlim([.5 2.5])

% 
text(1.2, -.05, 'reCFS', 'Rotation', -90, 'color', 'r', 'fontweight', 'bold',...
    'fontsize', 15)
% 
text(2.2, -.05, 'reCFS', 'Rotation', -90, 'color', 'r', 'fontweight', 'bold',...
    'fontsize', 15)


% add text:
text(1.25, oppD(isub-1,1)-.1, 'Invis.', 'Color', 'r', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold');

text(2.25, oppD(isub-1,2)-.1, 'Invis.', 'Color', 'r', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom', 'FontSize', 15, 'FontWeight','bold')

set(gca,'fontsize', 15)

end
%%


