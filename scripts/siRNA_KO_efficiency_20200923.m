% siRNA KO efficiency analysis

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
dat = struct();
dat.treated = readtable([path , 'BjornData/siRNA-treatedData.csv']);
% Data structure is:
% col 1: siRNA name
% col 2-5: siRNA-treated, gene of interest level
% col 6:9: siRNA-treated, actin level (col 6 matches col 2, col 7/3 etc);

% Note individual dilute samples used for control cells (so repeated actin
% measurements
dat.control = readtable([path , 'BjornData/siRNA-controlData.csv']);
% Data structure repeated as above % triplicate not quadruplicate
% Note columns 6-9 contain identical rows for the actin values to keep
% comparable structure.

% convert table to array
dat.mat.treated = table2array(dat.treated(:,2:9));
dat.mat.control = table2array(dat.control(:,2:7));

%% Calculate from Ct values

% Work out per sample dCT (sample minus actin)
dat.calcs.dCTCtrl  = dat.mat.control(:,1:3) - dat.mat.control(:,4:6);
dat.calcs.dCTsiRNA = dat.mat.treated(:,1:4) - dat.mat.treated(:,5:8);

% Wnt siRNA omitted(11,3) as clear outlier
dat.calcs.dCTsiRNA(11,3) = NaN;

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

siNames = dat.treated.Var1;
x = 1:numel(siNames);

% Want to reorder genes for consistency with other figures:

reorder = [1,4,5,3,2,7,6,8,14,12,9,10,13,11];

for ii = 1:numel(siNames)
siNames2{ii} = siNames{reorder(ii)};
end

dat.KO.mean = dat.KO.mean(reorder);
dat.KO.mSD = dat.KO.mSD(reorder);
dat.KO.pSD = dat.KO.pSD(reorder);
%%
figure
scatter(1:1:14 , dat.KO.mean,'filled', 'k')
% Altered from bar to scatter formatting for compliance with journal
% policy.
%bar(x , dat.KO.mean,'k','FaceAlpha',0.2)
hold on

errorbar(x , dat.KO.mean , dat.KO.mSD , dat.KO.pSD ,'Color','k', 'LineStyle' , 'none','LineWidth',1);
xticklabels(siNames2)
xtickangle(45)
ylabel({'mRNA knockdown efficiency','Relative to untreated'})
ylim([0 120])
xlim([0,15])
xlabel('siRNA')
xticks([1:1:14])
set(gca,'FontSize',14);

print([path , '/Figures/siRNA_KO_efficiency.pdf'],'-dpdf');

writematrix([dat.KO.mean';dat.KO.mSD';dat.KO.pSD'],'FigS17.csv')