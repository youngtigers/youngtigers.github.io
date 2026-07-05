

%% Load saved files
load Params
load store_Kstar_L1

p = Params.p;
tt = 1:1:80;
    if Kstar_L1 == 1
        load(['hat_alphag1_Application.mat']);
        load(['hat_beta1g1_Application.mat']);
        load(['hat_beta2g1_Application.mat']);
        load(['hat_beta3g1_Application.mat']);
        load(['hat_beta4g1_Application.mat']);
        load(['hat_beta5g1_Application.mat']);

    elseif Kstar_L1 == 2
        load(['hat_alphag1_Application.mat']);
        load(['hat_alphag2_Application.mat']);

        load(['hat_beta1g1_Application.mat']);
        load(['hat_beta1g2_Application.mat']);

        load(['hat_beta2g1_Application.mat']);
        load(['hat_beta2g2_Application.mat']);

        load(['hat_beta3g1_Application.mat']);
        load(['hat_beta3g2_Application.mat']);
        
        load(['hat_beta4g1_Application.mat']);
        load(['hat_beta4g2_Application.mat']);

        load(['hat_beta5g1_Application.mat']);
        load(['hat_beta5g2_Application.mat']);

    elseif Kstar_L1 == 3
        load(['hat_alphag1_Application.mat']);
        load(['hat_alphag2_Application.mat']);
        load(['hat_alphag3_Application.mat']);

        load(['hat_beta1g1_Application.mat']);
        load(['hat_beta1g2_Application.mat']);
        load(['hat_beta1g3_Application.mat']);

        load(['hat_beta2g1_Application.mat']);
        load(['hat_beta2g2_Application.mat']);
        load(['hat_beta2g3_Application.mat']);

        load(['hat_beta3g1_Application.mat']);
        load(['hat_beta3g2_Application.mat']);
        load(['hat_beta3g3_Application.mat']);

        load(['hat_beta4g1_Application.mat']);
        load(['hat_beta4g2_Application.mat']);
        load(['hat_beta4g3_Application.mat']);

        load(['hat_beta5g1_Application.mat']);
        load(['hat_beta5g2_Application.mat']);
        load(['hat_beta5g3_Application.mat']);
    elseif Kstar_L1 == 4
        load(['hat_alphag1_Application.mat']);
        load(['hat_alphag2_Application.mat']);
        load(['hat_alphag3_Application.mat']);
        load(['hat_alphag4_Application.mat']);

        load(['hat_beta1g1_Application.mat']);
        load(['hat_beta1g2_Application.mat']);
        load(['hat_beta1g3_Application.mat']);
        load(['hat_beta1g4_Application.mat']);

        load(['hat_beta2g1_Application.mat']);
        load(['hat_beta2g2_Application.mat']);
        load(['hat_beta2g3_Application.mat']);
        load(['hat_beta2g4_Application.mat']);

        load(['hat_beta3g1_Application.mat']);
        load(['hat_beta3g2_Application.mat']);
        load(['hat_beta3g3_Application.mat']);
        load(['hat_beta3g4_Application.mat']);

        load(['hat_beta4g1_Application.mat']);
        load(['hat_beta4g2_Application.mat']);
        load(['hat_beta4g3_Application.mat']);
        load(['hat_beta4g4_Application.mat']);

        load(['hat_beta5g1_Application.mat']);
        load(['hat_beta5g2_Application.mat']);
        load(['hat_beta5g3_Application.mat']);
        load(['hat_beta5g4_Application.mat']);
    end

%% Frontiers with standard errors 
dpi = 300;             % Resolution
ssz = [0 0 2400 1500]; % Image size in pixels
figure(...
  'PaperUnits','Centimeters', ...
  'PaperPositionMode','manual', ...
  'PaperPosition', [1.97 3.30 24.00 15.00], ...
  'PaperOrientation','landscape', ...
  'Visible', 'on')


% alphag1
subplot(2,6,1)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_alphag1_low' fliplr(hat_alphag1_high') ];
g1=plot(tt, hat_alphag1','b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
hline = refline([0 0]);
hline.Color = 'k';
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\alpha}_{(1)}(\tau_{t})$','interpreter','latex')

% alphag2
subplot(2,6,7)
qy=[ hat_alphag2_low' fliplr(hat_alphag2_high') ];
g2=plot(tt, hat_alphag2','r','LineWidth',1.2);
hold on
patch(qx, qy, 'r','LineStyle', 'none', 'EdgeColor','none')
hline = refline([0 0]);
hline.Color = 'k';
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\alpha}_{(2)}(\tau_{t})$','interpreter','latex')

% beta1g1
subplot(2,6,2)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta1g1_low' fliplr(hat_beta1g1_high') ];
g1=plot(tt, hat_beta1g1,'b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(1)1}(\tau_{t})$','interpreter','latex')

% beta1g2
subplot(2,6,8)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta1g2_low' fliplr(hat_beta1g2_high') ];
g2=plot(tt, hat_beta1g2,'r','LineWidth',1.2);
hold on
patch(qx, qy, 'r',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(2)1}(\tau_{t})$','interpreter','latex')

% beta2g1
subplot(2,6,3)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta2g1_low' fliplr(hat_beta2g1_high') ];
g1=plot(tt, hat_beta2g1,'b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(1)2}(\tau_{t})$','interpreter','latex')

% beta2g2
subplot(2,6,9)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta2g2_low' fliplr(hat_beta2g2_high') ];
g2=plot(tt, hat_beta2g2,'r','LineWidth',1.2);
hold on
patch(qx, qy, 'r',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(2)2}(\tau_{t})$','interpreter','latex')

% beta3g1
subplot(2,6,4)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta3g1_low' fliplr(hat_beta3g1_high') ];
g1=plot(tt, hat_beta3g1,'b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(1)3}(\tau_{t})$','interpreter','latex')

% beta3g2
subplot(2,6,10)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta3g2_low' fliplr(hat_beta3g2_high') ];
plot(tt, hat_beta3g2,'r','LineWidth',1.2);
hold on
patch(qx, qy, 'r',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(2)3}(\tau_{t})$','interpreter','latex')

% beta4g1
subplot(2,6,5)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta4g1_low' fliplr(hat_beta4g1_high') ];
g1=plot(tt, hat_beta4g1,'b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(1)4}(\tau_{t})$','interpreter','latex')

%beta4g2
subplot(2,6,11)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta4g2_low' fliplr(hat_beta4g2_high') ];
g2=plot(tt, hat_beta4g2,'r','LineWidth',1.2);
hold on
patch(qx, qy, 'r',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(2)4}(\tau_{t})$','interpreter','latex')

% beta5g1
subplot(2,6,6)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta5g1_low' fliplr(hat_beta5g1_high') ];
g1=plot(tt, hat_beta5g1,'b','LineWidth',1.2);
hold on
patch(qx, qy, 'b',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g1,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(1)5}(\tau_{t})$','interpreter','latex')

%beta5g2
subplot(2,6,12)
qx=[ tt fliplr(tt)]; % make closed patch
qy=[ hat_beta5g2_low' fliplr(hat_beta5g2_high') ];
g2=plot(tt, hat_beta5g2,'r','LineWidth',1.2);
hold on
patch(qx, qy, 'r',  'LineStyle', 'none', 'EdgeColor','none')
alpha(0.2)
uistack(g2,'top')
set(gca,'Box','on');
xlim([0 80])
xticks([0 40 80])
ytickformat('%.2f')
xlabel('$t$','interpreter','latex')
title('$\hat{\beta}_{(2)5}(\tau_{t})$','interpreter','latex')

print(gcf,'Grouped_Frontiers_Application_SE.pdf','-dpdf','-r300');

%%
vthetahat_d = vthetahat';
%% Grouping results based on pi hat
dpi = 300;             % Resolution
ssz = [0 0 2400 1500]; % Image size in pixels
figure(...
  'PaperUnits','Centimeters', ...
  'PaperPositionMode','manual', ...
  'PaperPosition', [1.97 3.30 24.00 15.00], ...
  'PaperOrientation','landscape', ...
  'Visible', 'on')

% pi1 vs sigv
subplot(2,3,1)
scatter(sqrt(vthetahat_d(Ksieves2==1,end)),vthetahat_d(Ksieves2==1,1),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(sqrt(vthetahat_d(Ksieves2==2,end)),vthetahat_d(Ksieves2==2,1),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
ylabel('$\hat{\pi}_{i1}$: coeff on $B_{1}(\tau_{t})$','interpreter','latex')
xlabel('$\hat{\sigma}_{vi}$','interpreter','latex')
legend('Group 1','Group 2','Location','southeast')
title('(a)','interpreter','latex')
ytickformat('%.1f')
xtickformat('%.2f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

% pi2 vs pi3
subplot(2,3,2)
scatter(vthetahat_d(Ksieves2==1,2),vthetahat_d(Ksieves2==1,3),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(vthetahat_d(Ksieves2==2,2),vthetahat_d(Ksieves2==2,3),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
xlabel('$\hat{\pi}_{i2}$: coeff on $w_{it1}$','interpreter','latex')
ylabel('$\hat{\pi}_{i3}$: coeff on $w_{it1}B_{1}(\tau_{t})$','interpreter','latex')
title('(b)','interpreter','latex')

ytickformat('%.2f')
xtickformat('%.1f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

% pi4 vs pi5
subplot(2,3,3)
scatter(vthetahat_d(Ksieves2==1,4),vthetahat_d(Ksieves2==1,5),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(vthetahat_d(Ksieves2==2,4),vthetahat_d(Ksieves2==2,5),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
xlabel('$\hat{\pi}_{i4}$: coeff on $w_{it2}$','interpreter','latex')
ylabel('$\hat{\pi}_{i5}$: coeff on $w_{it2}B_{1}(\tau_{t})$','interpreter','latex')
title('(c)','interpreter','latex')
ytickformat('%.2f')
xtickformat('%.1f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

% pi6 vs pi7
subplot(2,3,4)
scatter(vthetahat_d(Ksieves2==1,6),vthetahat_d(Ksieves2==1,7),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(vthetahat_d(Ksieves2==2,6),vthetahat_d(Ksieves2==2,7),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
xlabel('$\hat{\pi}_{i6}$: coeff on $y_{it1}$','interpreter','latex')
ylabel('$\hat{\pi}_{i7}$: coeff on $y_{it1}B_{1}(\tau_{t})$','interpreter','latex')
title('(d)','interpreter','latex')
ytickformat('%.1f')
xtickformat('%.1f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

% pi8 vs pi9
subplot(2,3,5)
scatter(vthetahat_d(Ksieves2==1,8),vthetahat_d(Ksieves2==1,9),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(vthetahat_d(Ksieves2==2,8),vthetahat_d(Ksieves2==2,9),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
xlabel('$\hat{\pi}_{i8}$: coeff on $y_{it2}$','interpreter','latex')
ylabel('$\hat{\pi}_{i9}$: coeff on $y_{it2}B_{1}(\tau_{t})$','interpreter','latex')
title('(e)','interpreter','latex')
ytickformat('%.1f')
xtickformat('%.1f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

% pi10 vs pi11
subplot(2,3,6)
scatter(vthetahat_d(Ksieves2==1,10),vthetahat_d(Ksieves2==1,11),12,'b','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
hold on
scatter(vthetahat_d(Ksieves2==2,10),vthetahat_d(Ksieves2==2,11),12,'r','filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4);
set(gca,'Box','on');
xlabel('$\hat{\pi}_{i10}$: coeff on $y_{it3}$','interpreter','latex')
ylabel('$\hat{\pi}_{i11}$: coeff on $y_{it3}B_{1}(\tau_{t})$','interpreter','latex')
title('(f)','interpreter','latex')
ytickformat('%.1f')
xtickformat('%.1f')
grid on 
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;

print(gcf,'Grouped_Frontiers_pihat.pdf','-dpdf','-r300');

%% Economies of scale
ttt = 1986:0.25:2005.75;
ES_g1 = 1./(hat_beta3g1 + hat_beta4g1 + hat_beta5g1) - 1;
ES_g2 = 1./(hat_beta3g2 + hat_beta4g2 + hat_beta5g2) - 1;

dpi = 300;             % Resolution
ssz = [0 0 2400 1500]; % Image size in pixels
figure(...
  'PaperUnits','Centimeters', ...
  'PaperPositionMode','manual', ...
  'PaperPosition', [1.97 3.30 16.00 10.00], ...
  'PaperOrientation','landscape', ...
  'Visible', 'on')

plot(ttt, ES_g1,'b','LineWidth',2);
hold on
plot(ttt, ES_g2,'r','LineWidth',2);
ytickformat('%.2f')
xlim([1986 2005.75])
xlabel('Year','interpreter','latex')
ylabel('Economies of scale','interpreter','latex')
legend('Group 1', 'Group 2','Location','northwest')
print(gcf,'ES_Application.pdf','-dpdf','-r300');
