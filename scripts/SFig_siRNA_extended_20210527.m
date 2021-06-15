% siRNA cell viability and Antiviral protease substrate SFig

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';


% Import the normalised (to scrambled control) cell viability data
uDat = readtable([path , 'BjornData/cellCounts_Uninfected.csv']);
iDat = readtable([path , 'BjornData/cellCounts_Infected.csv']);

uMat = readmatrix([path , 'BjornData/cellCounts_Uninfected.csv']);
iMat = readmatrix([path , 'BjornData/cellCounts_Infected.csv']);

% Transpose
uMat = uMat';
iMat = iMat';

uMat(1,:) = [];
iMat(1,:) = [];

% Names
grpNames = table2cell(uDat(:,1));

%% Import the RNA/Pfu/KO data

%
RNAdat = readtable([path , 'BjornData/siRNA_ext_RNA.csv']);
PFUdat = readtable([path , 'BjornData/siRNA_ext_PFU.csv']);

RNAmat = readmatrix([path , 'BjornData/siRNA_ext_RNA.csv']);
PFUmat = readmatrix([path , 'BjornData/siRNA_ext_PFU.csv']);

RNAmat = RNAmat';
PFUmat = PFUmat';

RNAmat(1,:) = [];
PFUmat(1,:) = [];

grpNames2 = table2cell(RNAdat(:,1));



%% Stats
% One-way ANOVA on RNA
[p , tbl , stats] = anova1(log10(RNAmat),grpNames2,'off');
RNAMC = multcompare(stats,'Display','off');
clear p tbl stats

% One-way ANOVA on PFU
[p , tbl , stats] = anova1(log10(PFUmat),grpNames2,'off');
PFUMC = multcompare(stats,'Display','off');
clear p tbl stats

% One-way ANOVA on Uninfected viability data
[p , tbl , stats] = anova1(uMat,grpNames,'off');
uMC = multcompare(stats,'Display','off');
clear p tbl stats

% One-way ANOVA on Infected viability data
[p , tbl , stats] = anova1(iMat,grpNames,'off');
iMC = multcompare(stats,'Display','off');
clear p tbl stats
%%
% Generate a grouping variable based on if meets p-value threshold (0.01)
% ... for RNA
rSig = RNAMC([1:4],6) <= 0.01; % Which have p-values <= 0.01
rSig = rSig + 2;
rSig = [1; rSig];

pSig = PFUMC([1:4],6) <= 0.01; % Which have p-values <= 0.01
pSig = pSig + 2;
pSig = [1; pSig];


uSig = uMC([1:16],6) <= 0.05; % Which have p-values <= 0.01
uSig = uSig + 2;
uSig = [1; uSig];

% ... for infected
iSig = iMC([1:16],6) <= 0.05; % Which have p-values <= 0.01
iSig = iSig + 2;
iSig = [1; iSig];

%% siRNA KO efficiency calculations
dat.treated = readtable([path , 'BjornData/siRNA_ext_treated.csv']);
% Data structure is:
% col 1: siRNA name
% col 2-7: siRNA-treated, gene of interest level
% col 8:13: siRNA-treated, actin level (col 6 matches col 2, col 7/3 etc);

% Note individual dilute samples used for control cells (so repeated actin
% measurements
dat.control = readtable([path , 'BjornData/siRNA_ext_untreated.csv']);
% Data structure repeated as above, but in duplicate

% convert table to array
dat.mat.treated = table2array(dat.treated(:,2:13));
dat.mat.control = table2array(dat.control(:,2:5));

% Calculate from Ct values

% Work out per sample dCT (sample minus actin)
dat.calcs.dCTCtrl  = dat.mat.control(:,1:2) - dat.mat.control(:,3:4);
dat.calcs.dCTsiRNA = dat.mat.treated(:,1:6) - dat.mat.treated(:,7:12);

% Wnt siRNA omitted(11,3) as clear outlier
%dat.calcs.dCTsiRNA(11,3) = NaN;

% Work out mean
dat.calcs.dCTCtrlMean  = mean(dat.calcs.dCTCtrl , 2 , 'omitnan');
dat.calcs.dCTsiRNAMean = mean(dat.calcs.dCTsiRNA, 2 , 'omitnan');

% Work out SDs - each is paired so no SD
dat.calcs.dCTCtrlSD  = std(dat.calcs.dCTCtrl , 0 , 2 , 'omitnan');
dat.calcs.dCTsiRNASD = std(dat.calcs.dCTsiRNA, 0 , 2 , 'omitnan');

% Work out ddCT
dat.calcs.ddCTmean = dat.calcs.dCTsiRNAMean - dat.calcs.dCTCtrlMean;
dat.calcs.ddCTSD   = sqrt(dat.calcs.dCTCtrlSD.^2 + dat.calcs.dCTsiRNASD.^2);

% Calculate fold-change
% = 2^ddCt
dat.calcs.ddCT_FC = 2.^dat.calcs.ddCTmean;

% SD = 2^SD
dat.calcs.ddCT_FCmSD = 2.^(dat.calcs.ddCTmean - dat.calcs.ddCTSD);
dat.calcs.ddCT_FCpSD = 2.^(dat.calcs.ddCTmean + dat.calcs.ddCTSD);

% Convert FC to KO %

dat.KO.mean = (1 - (1./dat.calcs.ddCT_FC))    .* 100;
dat.KO.mSD  = (1 - (1./dat.calcs.ddCT_FCmSD)) .* 100;
dat.KO.pSD  = (1 - (1./dat.calcs.ddCT_FCpSD)) .* 100;

dat.KO.mSD = dat.KO.mean - dat.KO.mSD;
dat.KO.pSD = dat.KO.pSD - dat.KO.mean;

% Now plot using bar and errorbar functions
siNames = dat.treated.Var1

% Want to reorder genes for consistency with other figures:

reorder = [1,2,4,3];

 for ii = 1:numel(siNames)
 siNames2{ii} = siNames{reorder(ii)};
 end
% 
 dat.KO.mean = dat.KO.mean(reorder);
 dat.KO.mSD = dat.KO.mSD(reorder);
 dat.KO.pSD = dat.KO.pSD(reorder);

%% Plotting

figure

%%%%% Todo: add KO, and RNA, Pfu
% RNA
subplot(3,3,1)
boxplot(log10(RNAmat),'PlotStyle','compact','Colors','krb','ColorGroup',rSig,'Labels',grpNames2,'Symbol','')
hold on
x = repmat(1:5,6,1);
scatter(x(:),log10(RNAmat(:)),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);
ylabel({'Log_1_0 RNA copies';'PFU equivalent/mL'})

text(4,7,'p \leq 0.01 vs. Ctrl','Color','b','FontSize',14);
text(4,6,'p > 0.01 vs. Ctrl','Color','r','FontSize',14);

hold off
set(gca,'FontSize',14);

% PFU
subplot(3,3,2)
boxplot(log10(PFUmat),'PlotStyle','compact','Colors','krb','ColorGroup',pSig,'Labels',grpNames2,'Symbol','')
hold on
x = repmat(1:5,6,1);
scatter(x(:),log10(PFUmat(:)),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);
ylabel({'Log_1_0 PFU/mL'})
% Add LOD
line([0 6.5],[log10(40),log10(40)],'Color','k','LineStyle',':')
text(4,7,'p \leq 0.01 vs. Ctrl','Color','b','FontSize',14);
text(4,6,'p > 0.01 vs. Ctrl','Color','r','FontSize',14);
hold off
set(gca,'FontSize',14);


% KO efficiency
subplot(3,3,3)
x = 1:numel(siNames);
bar(x , dat.KO.mean,'k','FaceAlpha',0.2)
hold on

errorbar(x , dat.KO.mean , dat.KO.mSD , dat.KO.pSD ,'Color','k', 'LineStyle' , 'none','LineWidth',1);
xticklabels(siNames2)
xtickangle(90)
ylabel({'mRNA knockdown efficiency','Relative to untreated'})
ylim([0 120])
%xlabel('siRNA')
set(gca,'FontSize',14);

% Now viability
subplot(3,3,4:6)
boxplot(uMat,'PlotStyle','compact','Colors','krb', 'ColorGroup',uSig,'Symbol','')
set(gca,'XTickLabel',[]);
hold on
x = repmat(1:17,6,1);
scatter(x(:),uMat(:),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);
ylim([0,250])

text(16,225,'p \leq 0.05 vs. Ctrl','Color','b','FontSize',14);
text(16,200,'p > 0.05 vs. Ctrl','Color','r','FontSize',14);
ylabel({'Relative cell count vs. Ctrl';'Uninfected'})
hold off
set(gca,'FontSize',14);

subplot(3,3,7:9)
boxplot(iMat,'PlotStyle','compact','Colors','krb', 'ColorGroup',iSig,'Labels',grpNames,'Symbol','')
hold on
x = repmat(1:17,5,1);
scatter(x(:),iMat(:),'filled','k','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.15);
ylim([0, 400])
xlabel('siRNA')

text(16,370,'p \leq 0.05 vs. Ctrl','Color','b','FontSize',14);
text(16,330,'p > 0.05 vs. Ctrl','Color','r','FontSize',14);
ylabel({'Relative cell count vs. Ctrl';'Infected'})
hold off
set(gca,'FontSize',14);

%% Supplemental analysis - IRF3
ctrl = PFUmat(:,1) ./ mean(PFUmat(:,1),'omitnan')

irf3 = PFUmat(:,2) ./ mean(PFUmat(:,1),'omitnan')

test = [ctrl,irf3] .* 100; % As % of Ctrl

[h,p,ci,stats] = ttest2(test(:,1),test(:,2));

%[p , tbl , stats] = anova1(test);
%tStats = multcompare(stats,'Display','off');
%clear p tbl stats


figure
boxplot(test,'PlotStyle','compact','Colors','kr','Labels',{'Ctrl','IRF3'});
hold on
x = repmat(1:2,6,1);
scatter(x(:),test(:),'filled','k','MarkerFaceAlpha',0.3,'jitter','on','jitterAmount',0.15);
ylabel('PFU as % of Ctrl')
text(1,300,['p = ',num2str(p)],'FontSize',14);
hold off

title('PFU in IRF3 siRNA-treated cells as % of Control')
xlabel('siRNA')

set(gca,'FontSize',14);