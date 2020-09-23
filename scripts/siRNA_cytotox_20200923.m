% siRNA cytotoxicity data analysis
 % For S. Figure

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
% Note - data analysed at Institut Pasteur - this script plots the
% previously analysed data
toxDat = readtable([path , 'BjornData/siRNA_cytotox.csv']); % RNA data (5uM)

% Convert table to matrix. Transpose so columns are groups
toxMat = table2array(toxDat(:,2:4))'; 


grpNames = table2cell(  toxDat(:,1)); % Identical for both

%% Plot results.
% Color bars by ANOVA p-value. blue = sig, red = NS.

% RNA data
figure
bar( mean(toxMat,'omitnan'),'k','FaceAlpha',0.2);

hold on
% Add error bar indicating the standard deviation
errorbar([1:14], mean(toxMat,'omitnan'), std(toxMat,'omitnan'),'k','LineStyle','none')
set(gca,'XTickLabel',[]);

% Overlay jittered scatter plot to display individual datapoints
x = repmat(1:14,3,1);
scatter(x(:),toxMat(:),'filled','r','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0.15);

ylabel({'Relative Cell viability'; '(% of untreated)'})
ylim([0 130])
xlabel('siRNA')
xticklabels(grpNames);
xtickangle(45);
hold off
set(gca,'FontSize',14);

% Save plot
print([path , '/Figures/siRNA_cytotox.pdf'],'-dpdf');






