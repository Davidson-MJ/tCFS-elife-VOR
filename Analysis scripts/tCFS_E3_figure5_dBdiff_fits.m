% tCFS_E3_figure4_dBdiff

% job=[];
% job.concatGFX= 1;
% job.plotPFX=0;
job.plotGFX=1;
job.compareFits=1;
job.plotGFX_oscillator=1;

% plotType==1; % as per other MS figure 
%% set up some directories.
% homedir = 'C:\Users\mdav0285\OneDrive - The University of Sydney (Staff)\Desktop\tracking-CFS';
homedir= '/Users/matthewdavidson/Documents/GitHub/tCFS-ver1';

cd(homedir)
cd(['DATA' filesep 'Tracking_Exp3']);% filesep 'MaskTracking_v1'])
datadir = pwd;
fntsize=15;

load('GFX_table.mat');

spdCols= {'Slow_', 'Medium_', 'Fast_'};
 sfx= {'', '_diff'};

%%

if job.plotGFX==1;
% compare each type (of 5).
barD = [GFX_table.DownThresholds_tCFS_Slow_diff, GFX_table.UpThresholds_tCFS_Slow_diff,...
    GFX_table.DownThresholds_tCFS_Normal_diff, GFX_table.UpThresholds_tCFS_Normal_diff,...
    GFX_table.DownThresholds_tCFS_Fast_diff, GFX_table.UpThresholds_tCFS_Fast_diff];

%% plot pretty MS figure:
nSubs= size(barD,1);
% first scatter the individual data points.

clf;
xM= [1 1 2 2 3 3 ];

% xM= xM +[ -.2 ,.2 ,-.2, .2 ,-.2,.2,-.2,.2,-.2,.2];

% xM = [.8, 1.2, 1.8, 2.2]; % space bad boys
nppants= size(barD,1);

mB = mean(barD,1);
mBar= [mB(1,[1,2]);...
   mB(1,[3,4]);
    mB(1,[5,6])];

bh=bar(1:length(mBar), mBar);
bh(1).FaceColor= 'r';
bh(1).FaceAlpha= .4;
bh(2).FaceColor= 'b';
bh(2).FaceAlpha= .4;
% xlim([.5 2.5])
shg

stE= CousineauSEM(barD);
stEbar= [stE(1,[1,2]);...
    stE(1,[3,4]);...
    stE(1,[5,6])]; %
errorbar_groupedfit(mBar, stEbar);

legend(bh, {'reCFS', 'bCFS'}, 'AutoUpdate','off', 'location', 'North')
%%
%
% clf
jit_width=.03;
offset= .1;
offset= offset.* [1, -1, 1, -1,1, -1,1, -1,1, -1];
scCols={'r', 'b', 'r','b','r','b','r','b','r','b'}
Xcollect=[];
for iX=1:6

    xAt= repmat(xM(iX), [1,nppants]);
    %random jitter per datapoint.
    
    jit= (-.5 + rand(1, nppants)) * jit_width;
    
    %scatter:
%   sc=  scatter(xAt+jit+offset(iX), barD(:,iX))
% sc.MarkerFaceColor= scCols{iX};
% sc.MarkerEdgeColor= 'k';

Xcollect(:,iX) = xAt+jit-offset(iX);
end

% now connect the dots:
dalpha= .9;
for idot= 1:nppants

plot([Xcollect(idot,1), Xcollect(idot,2)],  [barD(idot,1), barD(idot,2)], 'linew', 2 ,'color', [.8 .8 .8,dalpha]);
plot([Xcollect(idot,3), Xcollect(idot,4)],  [barD(idot,3), barD(idot,4)], 'linew', 2 ,'color', [.8 .8 .8, dalpha]);

plot([Xcollect(idot,5), Xcollect(idot,6)],  [barD(idot,5), barD(idot,6)], 'linew', 2 ,'color', [.8 .8 .8,dalpha]);


end
%%
ylabel('Change to target contrast (|dB|)');
set(gca,'Xtick', [1:3], 'XTickLabel', {'Slow', 'Medium', 'Fast'});
% legend([bD,rD], {'bCFS', 'reCFS'},'Location', 'South', 'Orientation', 'Horizontal');
set(gca,'fontsize',15);
% axis tight
xlim([.5 3.5])
% ylim([-40 0])
shg


for id=1:(size(barD,2)/2)

    useCols= [1,2] + 2*(id-1);
    [h,p,ci, stat]= ttest(barD(:, useCols(1)), barD(:, useCols(2)));
if p<.05
    txtmsg1 = ['\itt = ' sprintf('%.2f', stat.tstat) ];
    txtmsg2 = ['\itp \rm=' sprintf('%.2f', p)];
else
        txtmsg1 = 'ns';
        txtmsg2 = '';
end
text(id, 28, txtmsg1, 'fontsize',11,'HorizontalAlignment','center')
text(id, 27, txtmsg2, 'fontsize',11, 'HorizontalAlignment','center')
end


end


if job.compareFits==1

    % here we will compare the DHO, simple Oscillator, and Cubic
    % polynomial. 
    % trying to justify selection of the DHO.

    diffinBIC=[];
    diffinAIC=[];
for itype=2%1:2
for ispd=1:3

    figure(2);
subplot(4,3,ispd + 3*(itype-1)); % show flips by trial pos.
cla
barD = [GFX_table.([spdCols{ispd} 'dB_flip_1' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_2' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_3' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_4' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_5' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_6' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_7' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_8' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_9' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_10' sfx{itype}]),...
       GFX_table.([spdCols{ispd} 'dB_flip_11' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_12' sfx{itype}]),...
   GFX_table.([spdCols{ispd} 'dB_flip_13' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_14' sfx{itype}]),...
      GFX_table.([spdCols{ispd} 'dB_flip_15' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_16' sfx{itype}])];

   

%% set up input:
y= mean(barD,1); % data to fit
interpAt = linspace(1, length(y), 1000);
yInterp= interp1(1:length(y), y, interpAt);

yInterp= detrend(yInterp);
t= interpAt;

%% fit various models:
% first model: (simple hO).
 % fit hO
%     function yhat = harmonic_oscillator(beta, t)
%     omega = beta(1);
%     amplitude = beta(2);
%     phase = beta(3);
%     yhat = amplitude * sin(omega * t + phase);
%     end

beta0 = [2*pi, 6, 0];  % initial guess for parameters

[betahat_hO, resnorm_hO] = lsqcurvefit(@harmonic_oscillator, beta0, t, yInterp);
% predicted values at time t:
yhat_hO = harmonic_oscillator(betahat_hO, t); %(needed to compute residuals)


% next is Cubic:
% function yhat = Cubic_poly(beta, x)
%     a = beta(1);
%     b = beta(2);
%     c = beta(3);
%     yhat = a*x.^2 + b*x + c;
% end

beta0 = [1, 1, 1];  % initial guess for parameters
[betahat_qP, resnorm_qP] = lsqcurvefit(@cubic_poly, beta0, t, yInterp);
yhat_qP = cubic_poly(betahat_qP, t);

% next the DHO:
%  amp, damp coeffic, freq, phase shift 
x0 = [6, .01, 4,1];
[betahat_dhO, resnorm_dhO] = lsqcurvefit(@dampedHarmonic_oscillator, x0, t, yInterp);
yhat_dhO = dampedHarmonic_oscillator(betahat_dhO, t);

% compute the R2 of each fit (note this isnt ideal for non-linear funcs,
% but good enough for now).

r2_poly = 1 - sum((yInterp - yhat_qP).^2) / sum((yInterp - mean(yInterp)).^2);  % coefficient of determination for qp 
r2_ho = 1 - sum((yInterp - yhat_hO).^2) / sum((yInterp - mean(yInterp)).^2);  % coefficient of determination for harmonic oscillator
r2_dhO = 1 - sum((yInterp - yhat_dhO).^2) / sum((yInterp - mean(yInterp)).^2);  % coefficient of determination for qp 


fprintf('R^2 for Cubic Polynomial Model: %.3f\n', r2_poly)
fprintf('R^2 for Harmonic Oscillator Model: %.3f\n', r2_ho)
fprintf('R^2 for damped Harmonic oscillatorModel: %.3f\n', r2_dhO);

% lastly, compute AIC and BIC, to account for model complexity and GoF.
n = length(t);  % number of observations
k_Cubic = 3;  % number of parameters for Cubic polynomial
k_hO = 3;  % number of parameters for oscillator
k_dhO = 4; % params for dHO


resid_Cubic = yInterp- yhat_qP;  % residuals for Cubic polynomial
resid_hO = yInterp - yhat_hO;  % residuals for oscillator
resid_dhO = yInterp - yhat_dhO;  % residuals for dh oscillator

SSE_Cubic = sum(resid_Cubic.^2);  % sum of squared errors for Cubic polynomial
SSE_hO= sum(resid_hO.^2);  % 
SSE_dhO = sum(resid_dhO.^2); 

AIC_Cubic = n*log(SSE_Cubic/n) + 2*k_Cubic;  % AIC for Cubic polynomial
AIC_hO= n*log(SSE_hO/n) + 2*k_hO;  % AIC for hO
AIC_dhO= n*log(SSE_dhO/n) + 2*k_dhO;  % AIC for dhO

BIC_Cubic = n*log(SSE_Cubic/n) + k_Cubic*log(n);  % BIC for Cubic polynomial
BIC_hO= n*log(SSE_hO/n) + k_hO*log(n);  % BIC for h O
BIC_dhO= n*log(SSE_dhO/n) + k_dhO*log(n);  % BIC for d hO

% fprintf('AIC for Cubic Polynomial Model: %.3f\n', AIC_Cubic)
% fprintf('AIC for harmonic oscillator : %.3f\n', AIC_hO)
% fprintf('AIC for damped harmonic oscillator : %.3f\n', AIC_dhO);

fprintf('BIC for Cubic Polynomial Model: %.3f\n', BIC_Cubic)
fprintf('BIC for harmonic oscillator : %.3f\n', BIC_hO)
fprintf('BIC for damped harmonic oscillator : %.3f\n', BIC_dhO);


fprintf('BIC difference for (3-1) dHO vs Cubic: %.3f\n', ((BIC_dhO- BIC_Cubic)));
fprintf('BIC difference for (3-2)dHO vs oscillator: %.3f\n', ((BIC_dhO- BIC_hO)));
fprintf('BIC difference for (2-1) oscillator vs Cubic: %.3f\n', ((BIC_hO- BIC_Cubic)));


% fprintf('AIC difference for dHO vs Cubic: %.3f\n', ((AIC_dhO- AIC_Cubic)));
% fprintf('AIC difference for dHO vs oscillator: %.3f\n', ((AIC_dhO- AIC_hO)));



% store for later plots

% diffinAIC(ispd,1)= (AIC_dhO- AIC_hO);
% diffinAIC(ispd,2)= (AIC_dhO- AIC_Cubic);


diffinBIC(ispd,1)= (BIC_dhO- BIC_hO);
diffinBIC(ispd,2)= (BIC_dhO- BIC_Cubic);
diffinBIC(ispd,3)= (BIC_hO- BIC_Cubic);
% % display each:
% 
% clf
% plot(t, yInterp, 'k', 'linew',2);
% hold on;
% %
% plot(t, yhat_hO, 'r:', 'linew',2)
% 
% plot(t, yhat_qP, 'b:', 'linew',2)
% 
% plot(t, yhat_dhO, 'k:', 'linew',2)

end % ispeed.
end % itype


end






if job.plotGFX_oscillator==1
%% % resp # and harmonic oscillator fits.
figure(2); clf
set(gcf,'color', 'w', 'Units', 'normalized', 'position', [.1 .1 .7 .8])
figure(3); clf
set(gcf,'color', 'w', 'Units', 'normalized', 'position', [.1 .1 .7 .8])
 sfx= {'', '_diff'};
ylabsare= {{['Target'];['contrast (dB)']},{['Change to Target'];['contrast (|dB|)']}};
ylimsare=[-35 0;...
    0 25];

spdCols= {'Slow_', 'Medium_', 'Fast_'};
tsare = {'Slow RCC', 'Medium RCC', 'Fast RCC'};
for itype=1:2 % raw then diff.
for ispd=1:3

    figure(2);
barD = [GFX_table.([spdCols{ispd} 'dB_flip_1' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_2' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_3' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_4' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_5' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_6' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_7' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_8' sfx{itype}]),...
    GFX_table.([spdCols{ispd} 'dB_flip_9' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_10' sfx{itype}]),...
       GFX_table.([spdCols{ispd} 'dB_flip_11' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_12' sfx{itype}]),...
   GFX_table.([spdCols{ispd} 'dB_flip_13' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_14' sfx{itype}]),...
      GFX_table.([spdCols{ispd} 'dB_flip_15' sfx{itype}]), GFX_table.([spdCols{ispd} 'dB_flip_16' sfx{itype}])];

nflips= size(barD,2);
TargVisible= 2:2:nflips; % bCFS.
TargInvisible= 1:2:nflips; % reCFS.
tmpB = nan(1,nflips);
tmpB(TargVisible)= mean(barD(:, TargVisible),1);
tmpR = nan(1,nflips);
tmpR(TargInvisible)= mean(barD(:, TargInvisible),1);


if itype==1
    
subplot(4,3,ispd + 3*(itype-1)); % show flips by trial pos.
cla

br= plot(1:nflips, tmpR, 'rd','MarkerSize',10, 'LineWidth',2); hold on
bb= plot(1:nflips, tmpB, 'bd','MarkerSize',10, 'LineWidth',2);
hold on;
pc= plot(0:nflips, [0,mean(barD,1)], 'LineWidth',2);

else
    
    % 
    
subplot(4,3,ispd + 3*(1-1)); % show flips by trial pos.
% add subplot immediately beneath.
prevpos= get(gca,'position');
% width and height of prev sub:
wdth = prevpos(3);
hght= prevpos(4);
newpos = [prevpos(1), prevpos(2)-hght-.04, wdth, hght+.04];
subplot('position', newpos)
    
br=bar(1:nflips, tmpR,'FaceColor', 'r', 'FaceAlpha',.4); hold on;
bb=bar(1:nflips,tmpB,'FaceColor', 'b', 'FaceAlpha', .4); hold on;
hold on;
p=plot(1:nflips,  mean(barD,1), 'color', pc.Color, 'LineWidth',2);
xlabel('Response #')

end


set(gca,'fontsize', fntsize)
if itype==2 && ispd ==1
    legend([br, bb, p], {'reCFS', 'bCFS', 'Trend'}, 'autoupdate', 'off');
end
    
if itype==1
    title(tsare{ispd}, 'fontsize',15);
    set(gca,'xtick', [])
else
    set(gca,'xtick', 1:nflips, 'fontsize', fntsize-2)
end


stE= CousineauSEM(barD);
errorbar(1:nflips, mean(barD,1), stE, 'LineStyle', 'none', 'color', 'k', 'LineWidth',1)
xlim([-1  nflips+1])

ylim(ylimsare(itype,:))
if ispd==1
ylabel(ylabsare{itype})
end

% 
% text(-.75, 4, 'Starting', 'FontWeight','bold')
% text(-.75, 2, 'contrast:', 'FontWeight','bold')
% text(1, 2, '0 dB', 'HorizontalAlignment','center', 'FontWeight','bold');
% types= {'reCFS', 'bCFS'};
% for isp= 2:nflips
% 
% idx= mod(isp,2) +1;
% 
% text(isp, 4, 'prev', 'HorizontalAlignment','center', 'FontWeight','bold');
% text(isp, 2, types{idx}, 'HorizontalAlignment','center', 'FontWeight','bold');
% 
% end


% now do the fits:

% fir per ppant:
subFits = zeros(nSubs, 1000);
subRes= zeros(nSubs,1);
for isub= 1:size(barD,1);

    % t= 1:nflips;
    if itype==1
        y= [0,barD(isub,:)];
    else
        y=barD(isub,:);
    end
    [fitted,resnorm,t, yInterp] = myDHOfit(y);

    %evaluate fit:
  subFits(isub,:)=  yInterp;
  subRes(isub)= resnorm;
end








if itype==1
y= [0,mean(barD,1)];
else
y=mean(barD,1);
end

[fitted,resnorm,t, yInterp] = myDHOfit(y);


% figure(3)
% subplot(4,3,ispd+ 3*(itype-1))

% Plot the data and the fitted curve
% plot(t, propInterp, 'o');
% hold on;
% plot(t, oscillator(fitted, t), 'r-', 'linew',2);
% legend('Data', 'Fitted DHO');
% title({['eqn: x(1)*exp(-x(2)*t).*sin(x(3)*t + x(4))'];['coeffs: ' num2str(fitted)]})
% 
% xlabel('Response #');
% ylabel(ylabsare{itype})
% set(gca,'fontsize',fntsize);
% % also add the DHO to figure 2 for MS:

if itype==2
    figure(2);
loc= ispd+ 6*(itype-1);
    subplot(4,3,[loc, loc+3])

%plot the m and stE across subjs:
pM = mean(subFits,1);
stE= CousineauSEM(subFits);
sh=shadedErrorBar(t, pM, stE, {'color', pc.Color, 'linew',2},1)

% Plot the data and the fitted curve
% plot(t, yInterp, 'o', 'Color', pc.Color);
hold on;

% if ispd<3
    linets= {'-', ':'};
% else
%     linets= {':', '-'};
% end
    
pD=plot(t, dampedHarmonic_oscillator(fitted, t),'color', 'k', 'linew',2, 'linestyle', linets{1});

% add the others:
%sHO:
% beta0 = [2*pi, 6, 0,0];  % initial guess for parameters
% [betahat_hO, resnorm_hO] = lsqcurvefit(@harmonic_oscillator, beta0, t, yInterp);
% % predicted values at time t:
% yhat_hO = harmonic_oscillator(betahat_hO, t); %(needed to compute residuals)
% pO=plot(t, yhat_hO, 'b-', 'linew',2);

%cubic:

beta0 = [1, 1, 1];  % initial guess for parameters
[betahat_qP, resnorm_qP] = lsqcurvefit(@cubic_poly, beta0, t, yInterp);
yhat_qP = cubic_poly(betahat_qP, t);
pC=plot(t, yhat_qP ,'color','r', 'linew',2, 'linestyle', linets{2});



legend([sh.mainLine,pD, pC], {'Data', 'Fitted DHO', 'Fitted Cubic'});
% title({['eqn: x(1)*exp(-x(2)*t).*sin(x(3)*t + x(4))'];['coeffs: ' num2str(fitted)]})

xlabel('Response #');
if ispd==1
ylabel({['Detrended Change to'];['Target contrast (|dB|)']});
end
ylim([-4 4]);
% show delta AIC.
% dAIC = min(diffinAIC(ispd,:));
% dBIC = min(diffinBIC(ispd,:));
dBIC = diffinBIC(ispd,2);
%%
% text(8, -3.5, ['\DeltaAIC = ' sprintf('%.1f', dAIC)], 'EdgeColor', 'k');
text(4, -3.5, ['\DeltaBIC (DHO - Cubic) = ' sprintf('%.1f', dBIC)], 'EdgeColor', 'k');
set(gca,'fontsize',fntsize);
% xlim([0 10])
end
end % ispd
end % itype
%%
end





%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>
%%%%% FUNCTIONS CALLED
function [fitted,resnorm,t,yInterp] = myDHOfit(y)

% interpd ver:
interpAt = linspace(1, length(y), 1000);
yInterp= interp1(1:length(y), y, interpAt);

yInterp= detrend(yInterp);

% plot(interpAt, propInterp, 'o-', 'Color', pc.Color);

t = interpAt;
% y = exp(-0.5*t).*sin(2*pi*5*t);
oscillator=@(x,t)x(1)*exp(-x(2)*t).*sin(x(3)*t +x(4));

% Define initial parameter values
% amp, damp coeffic, freq, phase shift 
% x0 = [1, .1, 10, 90];

x0 = [6, .01, 4,1];

% Fit the oscillator equation to the data
[fitted,resnorm] = lsqcurvefit(oscillator, x0, interpAt, yInterp);

% return fitted, resnorm, t;

end
%%
function yhat= dampedHarmonic_oscillator(x0,t)
%  amp, damp coeffic, freq, phase shift 
% Define oscillator equation
yhat= x0(1)*exp(-x0(2)*t).*sin(x0(3)*t +x0(4));

end
 function yhat = harmonic_oscillator(beta0, t)
    omega = beta0(1);
    amplitude = beta0(2);
    phase = beta0(3);
%     d= beta0(4);
    yhat = amplitude * sin(omega * t + phase);% + d;
 end
%%
  function yhat = cubic_poly(beta0, x)
    a = beta0(1);
    b = beta0(2);
    c = beta0(3);
%     d = beta0(4);
    yhat = a*x.^3 + b*x.^2 + c*x;% + d;
end