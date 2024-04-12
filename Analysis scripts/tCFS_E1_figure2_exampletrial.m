% this script plots example trials for Jacobs E1.
clear all;
close all;


%% which jobs to complete?
job=[];
job.concatGFX= 1;
job.plotPFX=1;
job.plotGFX=1;
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir='/Users/matthewdavidson/Documents/GitHub/tCFS-ver1/'
cd(homedir)
cd(['DATA' filesep 'E1 mat files' filesep 'Tracking'])
datadir = pwd;
%

ippant = 'AK';

cd(ippant);

lfile = dir([pwd filesep '*.mat']);
load(lfile.name, 'contrProfile', 'allProfileCnts');

%%
% clf;
% for itrial= 1:8
%     subplot(2,4,itrial)
% pData = contrProfile{1,itrial};
% tax = [1:length(pData)]./60;
% plot(tax, pData);
% end
clf;
usetrial = 6; % good one for AK
contrData = contrProfile{1,usetrial};
tax = [1:length(contrData)]./60;

hold on;
turnsAt = allProfileCnts(usetrial,:);
bCFS = 2:2:length(turnsAt);
reCFS = 1:2:length(turnsAt);

bCFS_turns = turnsAt(bCFS);
reCFS_turns = turnsAt(reCFS);

plotData = {contrData , 20*log10(contrData)};
ysare = {'Target contrast', 'Target contrast (dB)'};
for id=1:2
subplot(3,2,id);

pData= plotData{id};
plot(tax, pData, 'k', 'LineWidth', 1); hold on;
ylabel(ysare{id});
xlabel('Time (sec)');

plot(tax(reCFS_turns), pData(reCFS_turns), 'rx', 'LineWidth',1)

rl=plot(tax(reCFS_turns), pData(reCFS_turns), 'ro','LineWidth',1);

bl=plot(tax(bCFS_turns), pData(bCFS_turns), 'bo','LineWidth',1);

if id==1
legend([bl,rl], {'b-CFS', 're-CFS'}, 'fontsize', 12)
end
shg;
set(gcf, 'color', 'w','units','normalized','position', [.1 .1 .5 .85]);
set(gca,'fontsize',18 )

end
%%

shg
% tax = 1:length(contrPr)