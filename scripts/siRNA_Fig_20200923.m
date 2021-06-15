% siRNA data analysis

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
rnaDat = readtable([path , 'BjornData/siRNA_RNApfuEquiv.csv']); % RNA data (5uM)
pfuDat = readtable([path , 'BjornData/siRNA_pfu.csv']); % Pfu data (5um)

% Convert RNA copies/pfu portion of table to matrix. Transpose so columns
% are groups
rnaMat = table2array(rnaDat(:,2:7))'; 
pfuMat = table2array(pfuDat(:,2:7))'; 

grpNames = table2cell(  rnaDat(:,1)); % Identical for both

% One-way ANOVA on RNA data
[p , tbl , stats] = anova1(log10(rnaMat),grpNames,'off');
rnaMC = multcompare(stats,'Display','off');
clear p tbl stats

% One-way ANOVA on Pfu data
[p , tbl , stats] = anova1(log10(pfuMat),grpNames,'off');
pfuMC = multcompare(stats,'Display','off');
clear p tbl stats

% Generate a grouping variable based on if meets p-value threshold
% ... for RNA
rnaSig = rnaMC([1:14],6) <= 0.01; % Which have p-values <= 0.01
rnaSig = rnaSig + 2;
rnaSig = [1; rnaSig];

% ... for Pfu
pfuSig = pfuMC([1:14],6) <= 0.01; % Which have p-values <= 0.01
pfuSig = pfuSig + 2;
pfuSig = [1; pfuSig];

%% Plot results.
% Color bars by ANOVA p-value. blue = sig, red = NS.

figure
% Plot RNA data
subplot(2,1,1) % tiledlayout appears to not cooperate with boxplot
boxplot(log10(rnaMat),'PlotStyle','compact','Colors','krb', 'ColorGroup',rnaSig);
set(gca,'XTickLabel',[]);
hold on
x = repmat(1:15,6,1);
scatter(x(:),log10(rnaMat(:)),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);


text(11,1.75,'p \leq 0.01 vs. Ctrl','Color','b','FontSize',14);
text(11,1.25,'p > 0.01 vs. Ctrl','Color','r','FontSize',14);
ylabel({'Log_1_0 RNA copies'; 'PFU equivalent/mL'})
ylim([0,7])
hold off
set(gca,'FontSize',14);

% Plot Pfu Data
subplot(2,1,2)
boxplot(log10(pfuMat),'PlotStyle','compact','Colors','kb', 'ColorGroup',pfuSig,'Labels',grpNames);
% Note: as no non-significant titre drops, had to adjust colors to remove
% the red non-significant category from as boxplot otherwise treats the 3
% as 2 and plots with red rather than blue.
hold on
x2 = repmat(1:15,6,1);
scatter(x2(:),log10(pfuMat(:)),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);

text(11,5.5,'p \leq 0.01 vs. Ctrl','Color','b','FontSize',14);
text(11,5,'p > 0.01 vs. Ctrl','Color','r','FontSize',14);
text(11,4.5,'...Limit of detection','Color','k','FontSize',14)
ylabel({'Log_1_0 PFU/mL'})
xlabel('siRNA')
set(gca,'FontSize',14);
line([0 15],[log10(40),log10(40)],'Color','k','LineStyle',':')
ylim([1,6])
% Save figure
print([path , '/Figures/Fig_siRNA.pdf'],'-dpdf');




