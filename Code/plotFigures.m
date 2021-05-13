%% plot figures tutorial paper

% clear all
% close all
% clc

%% load results
plotR(1) = load('../Results/Results_SEIR_MPC_SINDy_C0.mat','Results');
plotR(2) = load('../Results/Results_SEIR_MPC_SINDy_C1.mat','Results');
plotR(3) = load('../Results/Results_SEIR_MPC_DMD_C0.mat','Results');
plotR(4) = load('../Results/Results_SEIR_MPC_DMD_C1.mat','Results');

% define plot properties
lw = 3; % linewidth 1
lwV = 2; % linewidth 2
fontsizeLabel = 16;
fontsizeLegend = 12;
fontsizeTicks = 12;
fontsizeTitel = 16;
sizeX = 1800; % size figure
sizeY = 900; % size figure

% color propreties
C1 = [0 119 187]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255;
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;

beta0 = 0.5; % SEIR model beta parameter 
trE = plotR(1).Results.tTrain(end); % end of training time
teE = plotR(1).Results.tTest(end); % end of testing time
tME = plotR(1).Results.t(end); % end of MPC time
constraint = 0.05; % MPC constraint


%% plot system ID: comparison true vs SINDy vs DMD

% plot number of infectious cases
yLimM = 0.3;
figure('Position', [10 10 sizeX sizeY])
subplot(2,2,1)
box on
plot([plotR(3).Results.tTrain, trE+plotR(3).Results.tTest],[plotR(3).Results.xTrueTrain(:,3);plotR(3).Results.xTrueTest(:,3)],'-','Color',C5,'LineWidth',lw*1.5); hold on
plot([plotR(3).Results.tTrain, trE+plotR(3).Results.tTest],[plotR(3).Results.xModelTrain(:,3);plotR(3).Results.xModelTest(:,3)],'-','Color',C9,'LineWidth',lw);
plot([plotR(1).Results.tTrain, trE+plotR(1).Results.tTest],[plotR(1).Results.xModelTrain(:,3);plotR(1).Results.xModelTest(:,3)],'--','Color',C2,'LineWidth',lw);
plot([trE trE],[0 yLimM],'k--','Linewidth',lwV);
area([0 trE],[yLimM yLimM],'Facecolor','k','Facealpha',0.05);

xlim([0 trE+teE])
ylim([0 yLimM])
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('time, days','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$\#$ of infectious cases','interpreter','latex','fontsize',fontsizeLabel)
title({'a) SEIR system identification: DMD vs. SINDy';''},'fontsize',fontsizeTitel,'interpreter','latex')
% set(gca,'Xticklabel',[])
legend('true dynamics','DMD','SINDy','interpreter','latex','fontsize',fontsizeLegend,'Location','west')


% plot PRBS forcing
yLimM = beta0*0.8;
subplot(2,2,3)
box on
plot([plotR(1).Results.tTrain, trE+plotR(1).Results.tTest],beta0-[plotR(1).Results.uTrain; plotR(1).Results.uTest],'-','Color',C9,'LineWidth',lw); hold on
plot([trE trE],[0 yLimM],'k--','Linewidth',lwV);
area([0 trE],[yLimM yLimM],'Facecolor','k','Facealpha',0.05)

xlim([0 trE+teE])
ylim([0 yLimM])
yticks([0 0.5*yLimM yLimM])
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('time, days','interpreter','latex','fontsize',fontsizeLabel)
ylabel('control input','interpreter','latex','fontsize',fontsizeLabel)
legend('$\beta_0$ - u(t)','interpreter','latex','fontsize',fontsizeLegend,'Location','Southeast')
title({'';''},'fontsize',fontsizeTitel,'interpreter','latex')


%% plot comparison MPC SINDy vs DMD (with and without constraint)

% plot number of infectious cases
subplot(2,2,2)
box on
plot(plotR(1).Results.tUnforced,plotR(1).Results.xUnforced(:,3),'-','Color',C5,'LineWidth',lw); hold on
plot(plotR(3).Results.t,plotR(3).Results.x(3,:),'-','Color',C9,'LineWidth',lw);
plot(plotR(4).Results.t,plotR(4).Results.x(3,:),'-','Color',C7,'LineWidth',lw);
plot(plotR(1).Results.t,plotR(1).Results.x(3,:),'-','Color',C2,'LineWidth',lw);
plot(plotR(2).Results.t,plotR(2).Results.x(3,:),'-','Color',C10,'LineWidth',lw);
plot([0 tME],[constraint constraint],'--','Color',C8,'LineWidth',lw);
xlim([0 tME])
ylim([0 0.15])

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('time, days','interpreter','latex','fontsize',fontsizeLabel)
ylabel('$\#$ of infectious cases','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend('no control','DMD-MPC (no constraint)','DMD-MPC$_c$ (constraint)','SINDy-MPC (no constraint)','SINDy-MPC$_c$ (constraint)','constraint: I $\leq$ 0.05','interpreter','latex','fontsize',fontsizeLegend)
title({'b) MPC: DMD vs. SINDy';''},'fontsize',fontsizeTitel,'interpreter','latex')


% plot MPC control SINDy vs DMD (control u)
subplot(2,2,4)
box on
plot(plotR(3).Results.t,beta0-plotR(3).Results.u,'-','Color',C9,'LineWidth',lw); hold on
plot(plotR(4).Results.t,beta0-plotR(4).Results.u,'-','Color',C7,'LineWidth',lw);
plot(plotR(1).Results.t,beta0-plotR(1).Results.u,'-','Color',C2,'LineWidth',lw);
plot(plotR(2).Results.t,beta0-plotR(2).Results.u,'-','Color',C10,'LineWidth',lw);
xlim([0 tME])
ylim([0 yLimM])
yticks([0 0.5*yLimM yLimM])

set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
xlabel('time, days','interpreter','latex','fontsize',fontsizeLabel)
ylabel('control input','interpreter','latex','fontsize',fontsizeLabel)
title({'';''},'fontsize',fontsizeTitel,'interpreter','latex')
