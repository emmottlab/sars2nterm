% Fig 3. Replotting data as boxplot/scatter combinations

clear
clc

% Fig3CnormalisedDat % Lacking TMPRSS2
% Fig3DnormalisedDat % TMPRSS2

dat = struct();
path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Import data
dat.d3c = readtable([path,'Fig3CnormalisedDat.csv']);
dat.Md3c = table2array(dat.d3c);

dat.d3d = readtable([path,'Fig3DnormalisedDat.csv']);
dat.Md3d = table2array(dat.d3d);

dat.d3f = readtable([path,'3FnormalisedDensitometry.csv']);
dat.Md3f = table2array(dat.d3f([1:6],[2:7]))';
%% Stats
% Comparisons relative to WT spike

% HEK-Ace2 stats
% Two sample t-tests, without assuming equal variance (using Satterthwaite?s
% approximation for the effective degrees of freedom)
[~,p(1)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,2),'Vartype','unequal');
[~,p(2)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,3),'Vartype','unequal');
[~,p(3)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,4),'Vartype','unequal');
[~,p(4)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,5),'Vartype','unequal');
[~,p(5)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,6),'Vartype','unequal');
[~,p(6)] = ttest2(dat.Md3c(:,1),dat.Md3c(:,7),'Vartype','unequal');

pGrp = p<= 0.05;
pGrp = pGrp + 2;
pGrp = [1, pGrp];

% HEK-Ace2-TMPRSS2
% Two samplet-tests, without assuming equal variance (using Satterthwaite?s
% approximation for the effective degrees of freedom) 
[~,p2(1)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,2),'Vartype','unequal');
[~,p2(2)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,3),'Vartype','unequal');
[~,p2(3)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,4),'Vartype','unequal');
[~,p2(4)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,5),'Vartype','unequal');
[~,p2(5)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,6),'Vartype','unequal');
[~,p2(6)] = ttest2(dat.Md3d(:,1),dat.Md3d(:,7),'Vartype','unequal');

p2Grp = p2<= 0.05;
p2Grp = p2Grp + 2;
p2Grp = [1, p2Grp];

% WB densitometry
% Two samplet-tests, without assuming equal variance (using Satterthwaite?s
% approximation for the effective degrees of freedom) 
[~,p3(1)] = ttest2(dat.Md3f(:,1),dat.Md3f(:,2),'Vartype','unequal');
[~,p3(2)] = ttest2(dat.Md3f(:,1),dat.Md3f(:,3),'Vartype','unequal');
[~,p3(3)] = ttest2(dat.Md3f(:,1),dat.Md3f(:,4),'Vartype','unequal');
[~,p3(4)] = ttest2(dat.Md3f(:,1),dat.Md3f(:,5),'Vartype','unequal');
[~,p3(5)] = ttest2(dat.Md3f(:,1),dat.Md3f(:,6),'Vartype','unequal');


p3Grp = p3<= 0.05;
p3Grp = p3Grp + 2;
p3Grp = [1, p3Grp];
%% Plotting
% Note the split axis function (breakyaxis) replicates the sample names over the plot.
% These are readily removed in inkscape or equivalent.

labs1 = dat.d3c.Properties.VariableNames;
labs1{1,6} = 'V635G/C671G';
labs1{1,7} = 'VSV-G';

figure
subplot(2,2,1)
boxplot(dat.Md3c,'PlotStyle','compact','Labels',labs1,'Colors','krb','ColorGroup',pGrp,'Symbol','');
hold on
line([0,8],[1,1],'Color','k','LineStyle','--')
x = repmat(1:7,6,1);
scatter(x(:),dat.Md3c(:),'filled','k','MarkerFaceAlpha',0.25','jitter','on','jitterAmount',0.25);
ylim([0,35])
text(2,33,'p < 0.05 vs. WT','Color','b','FontSize',14);
text(2,32,'p = n.s. vs. WT','Color','r','FontSize',14);
hold off
ylabel('Normalised Infectivity')
title('Infectivity in HEK-Ace2 cells')
xtickangle(45)
breakyaxis([11,30],0.005)

subplot(2,2,2)
boxplot(dat.Md3d,'PlotStyle','compact','Labels',labs1,'Colors','krb','ColorGroup',p2Grp,'Symbol','');
hold on
line([0,8],[1,1],'Color','k','LineStyle','--')
x = repmat(1:7,6,1);
scatter(x(:),dat.Md3d(:),'filled','k','MarkerFaceAlpha',0.25','jitter','on','jitterAmount',0.25);
ylim([0,35])
text(2,33,'p < 0.05 vs. WT','Color','b','FontSize',14);
text(2,32,'p = n.s. vs. WT','Color','r','FontSize',14);
hold off
ylabel('Normalised Infectivity')
title('Infectivity in HEK-Ace2-TMPRSS2 cells')
xtickangle(45)
breakyaxis([11,30],0.005)

labs2 = dat.d3f.Var1(1:6);
labs2{6,1} = 'V635G/C671G';

subplot(2,2,4)
boxplot(dat.Md3f,'PlotStyle','compact','Labels',labs2,'Colors','krb','ColorGroup',p3Grp,'Symbol','');
hold on
line([0,7],[1,1],'Color','k','LineStyle','--')
x = repmat(1:6,6,1);
scatter(x(:),dat.Md3f(:),'filled','k','MarkerFaceAlpha',0.25,'jitter','on','jitterAmount',0.25);
ylim([0,4.5])
text(2,4.3,'p < 0.05 vs. WT','Color','b','FontSize',14);
text(2,4,'p = n.s. vs. WT','Color','r','FontSize',14);
hold off
ylabel('S1/S0 ratio')
title({'Ratio of cleaved (S1) to ';'uncleaved (S0) SARS-CoV-2 spike'})
xtickangle(45)

