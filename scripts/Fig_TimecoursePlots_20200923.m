% Figure 1 data: script plots titre, qRT-PCR and protein-level data as well
% as Volcano N-terminomics plots

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

%% Data Import and Processing for N-terminomics/Volcano Plots
dat = struct();

% Load the peptide files
dat.pep.A549 = readtable([path , 'tA549_Enrich/peptides.txt']);
dat.pep.Vero = readtable([path , 'tVero_Enrich/peptides.txt']);

samples = {'A549','Vero'};

% Load in the TMTpro design matrices containing the randomised layouts
% A549-ACE2
dat.tmt.A549 = readmatrix([path , 'SARS2a549tmtlabelling_20200507.csv']);
% VeroE6
dat.tmt.Vero = readmatrix([path , 'Verotmtlabelling_20200511.csv']);

% Rearrange the Reporter channels in dat.evi to obtain the correct order of
% Channels:

% Reorder the corrected RI intensity channels
dat.pep.A549(:,44:59) = dat.pep.A549(: , dat.tmt.A549(2,:) + 43);
dat.pep.Vero(:,44:59) = dat.pep.Vero(: , dat.tmt.Vero(2,:) + 43);

% Remove Reverse hits, Contaminants and low PEP identifications
for ii = 1:numel(samples)
    % Remove PEP <= 0.02
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).PEP <= 0.02 , :);
    % Remove Reverse hits
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(categorical(dat.pep.(samples{ii}).Reverse) ~= '+' , :);
    % Remove potential contaminants
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(categorical(dat.pep.(samples{ii}).PotentialContaminant) ~= '+' , :);
end

% Now Convert the TMT RI data to matrix form
for ii = 1:numel(samples)
   % Convert table columns to a matrix
    dat.mat.(samples{ii}) = table2array(dat.pep.(samples{ii})(: , [44:59]));
   
   % Convert 0 to NaN;
   dat.mat.(samples{ii})(dat.mat.(samples{ii}) == 0) = NaN; 
end

% Remove unquantified hits, column normalise by column median, knn impute and row normalise.
for ii = 1:numel(samples)
   % identify rows with >75% NaN
   allNaN = sum(isnan(dat.mat.(samples{ii})) , 2) > 13;
   dat.mat.(samples{ii}) = dat.mat.(samples{ii})(~allNaN , :);
   dat.pep.(samples{ii}) = dat.pep.(samples{ii})(~allNaN , :);
   
   clear allNaN
   
   % Median normalise to control for protein loading
    dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ nanmedian(dat.mat.(samples{ii}));

    % KNNimpute missing data
    dat.mat.(samples{ii}) = knnimpute(dat.mat.(samples{ii}));

    % Row normalise by mean (whole dataset) 
    dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ mean(dat.mat.(samples{ii}) , 2);

    % Generate matrix with 24h data only
    dat.h24.(samples{ii}) = dat.mat.(samples{ii})(:,[end-5: end]);
   
end

lia = contains(dat.pep.A549.Proteins , 'SARS2');

% Keep N-termini where the position in protein is greater than or equal 2.
 for ii = 1:numel(samples)
    dat.h24.(samples{ii}) = dat.h24.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
    dat.mat.(samples{ii}) = dat.mat.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
 end

% Perform t-test for significance
[~ , dat.p.A549 , ~ , ~] = ttest2(log2(dat.h24.A549(:,[1:3])) , log2(dat.h24.A549(:,[4:6])),'Dim',2);
[~ , dat.p.Vero , ~ , ~] = ttest2(log2(dat.h24.Vero(:,[1:3])) , log2(dat.h24.Vero(:,[4:6])),'Dim',2);

%Correct for multiple hypothesis testing
[dat.fdr.A549 , dat.q.A549] = mafdr(dat.p.A549);
[dat.fdr.Vero , dat.q.Vero] = mafdr(dat.p.Vero);

% Calculate mean fold-chage
dat.mean.A549 = mean(dat.h24.A549(:,[1:3]),2) ./ mean(dat.h24.A549(:,[4:6]) , 2);
dat.mean.Vero = mean(dat.h24.Vero(:,[1:3]),2) ./ mean(dat.h24.Vero(:,[4:6]) , 2);

% Plot settindgs
pmax = -5; % Y axis cutoffs: -Q values greater than log10(pmax) will show as pmax
qcut = 0.05; % Q value cutoff: minimum Q value for significance.
fccut = 0.5; % Fold-change cutoff: minimum fold change for hits
xcut = 6; % X axis cutoff: Fold-change of greater than xcut will show as xcut.

pcutoff = 10^pmax;

for ii = 1:numel(samples)
dat.q2.(samples{ii}) = dat.q.(samples{ii});
dat.q2.(samples{ii})(dat.q2.(samples{ii}) <= pcutoff) = pcutoff;
dat.mean.(samples{ii})(dat.mean.(samples{ii}) <= -2^xcut) = -2^xcut;
dat.mean.(samples{ii})(dat.mean.(samples{ii}) >= 2^xcut) = 2^xcut;
end

% Identify those hits which meet the fccut and qcut
for ii = 1:numel(samples)
   dat.qcut.(samples{ii}) = dat.q2.(samples{ii}) <= qcut;
   dat.xcutp.(samples{ii}) = log2(dat.mean.(samples{ii})) >=fccut;
   dat.xcutn.(samples{ii}) = log2(dat.mean.(samples{ii})) <= -fccut;
   dat.cut.(samples{ii}) = dat.xcutp.(samples{ii}) + dat.xcutn.(samples{ii}) >= 1;
   
   dat.hits.(samples{ii}) = dat.qcut.(samples{ii}) + dat.cut.(samples{ii}) == 2;
   dat.viral.(samples{ii}) = contains(dat.pep.(samples{ii}).Proteins , 'SARS2');
end

%% Data Import and processing for protein-level data.
datP = struct();

datP.prot.A549 = readtable([path , 'tA549_TotalFrac/proteinGroups.txt']);
datP.prot.Vero = readtable([path , 'tVero_TotalFrac/proteinGroups.txt']);

% RI in 19:34

% Rearrange the Reporter channels in datP.prot to obtain the correct order of
% Channels:

% Reorder the corrected RI intensity channels
datP.prot.A549(:,19:34) = datP.prot.A549(: , dat.tmt.A549(2,:) + 18);
datP.prot.Vero(:,19:34) = datP.prot.Vero(: , dat.tmt.Vero(2,:) + 18);

% Now Convert the TMT RI data to matrix form
for ii = 1:numel(samples)
   % Convert table columns to a matrix
    datP.mat.(samples{ii}) = table2array(datP.prot.(samples{ii})(: , [19:34]));
   
   % Convert 0 to NaN;
   datP.mat.(samples{ii})(datP.mat.(samples{ii}) == 0) = NaN; 
end

% Remove unquantified hits, column normalise, impute and row normalise.
for ii = 1:numel(samples)
   % identify rows with >75% NaN
   allNaN = sum(isnan(datP.mat.(samples{ii})) , 2) > 13;
   datP.mat.(samples{ii}) = datP.mat.(samples{ii})(~allNaN , :);
   datP.prot.(samples{ii}) = datP.prot.(samples{ii})(~allNaN , :);
   
   clear allNaN
   
   % Median normalise
    datP.mat.(samples{ii}) = datP.mat.(samples{ii}) ./ nanmedian(datP.mat.(samples{ii}));

    % KNNimpute missing data
    datP.mat.(samples{ii}) = knnimpute(datP.mat.(samples{ii}));

    % Row normalise by mean (whole dataset) 
    datP.mat.(samples{ii}) = datP.mat.(samples{ii}) ./ mean(datP.mat.(samples{ii}) , 2);
  
end

% Now want to search out the data for spike and nsp8 (High/Low expression)
datP.loc.AS = contains(datP.prot.A549.ProteinIDs,'SARS2|S|Spike');
datP.loc.VS = contains(datP.prot.Vero.ProteinIDs,'SARS2|S|Spike');
datP.loc.Ansp8 = contains(datP.prot.A549.ProteinIDs,'SARS2|nsp8|');
datP.loc.Vnsp8 = contains(datP.prot.Vero.ProteinIDs,'SARS2|nsp8|');

% Multiple the protein intensity by the normalised TMT channels
datP.s.A549 = datP.prot.A549.Intensity(datP.loc.AS) * datP.mat.A549(datP.loc.AS,:);
datP.s.Vero = datP.prot.Vero.Intensity(datP.loc.VS) * datP.mat.Vero(datP.loc.VS,:);

datP.nsp8.A549 = datP.prot.A549.Intensity(datP.loc.Ansp8) * datP.mat.A549(datP.loc.Ansp8,:);
datP.nsp8.Vero = datP.prot.Vero.Intensity(datP.loc.Vnsp8) * datP.mat.Vero(datP.loc.Vnsp8,:);

% Now reorder this data to make it easier for plotting
datP.Rs.A549 = log2(datP.s.A549([2,3,4;5,6,7;8,9,10;11,12,13]));
datP.Rs.Vero = log2(datP.s.Vero([2,3,4;5,6,7;8,9,10;11,12,13]));
datP.Rnsp8.A549 = log2(datP.nsp8.A549([2,3,4;5,6,7;8,9,10;11,12,13]));
datP.Rnsp8.Vero = log2(datP.nsp8.Vero([2,3,4;5,6,7;8,9,10;11,12,13]));


%% Load in timecourse data files (non-protein)

% Timecourse Data
timecourseTimes = [0,3,6,9,12,24];
timecourseTimesP = [0,6,12,24];
tcTicks = [0,3,6,9,12,15,18,21,24];
tcTicklabels = {'0','3','6','9','12','','18','','24'};

% Timecourse: qRT-PCR Data
tc_r_V = readmatrix([path , 'BjornData/TC_Vero_SupRNA.csv'])';
tc_r_A = readmatrix([path , 'BjornData/TC_A549_SupRNA.csv'])';

% Timecourse: Plaque Assay Data
tc_p_V = readmatrix([path , 'BjornData/TC_Vero_Plaque.csv'])';
tc_p_A = readmatrix([path , 'BjornData/TC_A549_Plaque.csv'])';

%% Generate text annotations for volcano plots (so only significant/viral hits)
dat.pep.A549.GN = dat.pep.A549.LeadingRazorProtein;
dat.pep.Vero.GN = dat.pep.Vero.LeadingRazorProtein;

% Identify SARS2 proteins
dat.S2.A549 = contains(dat.pep.A549.LeadingRazorProtein,'SARS2');
dat.S2.Vero = contains(dat.pep.Vero.LeadingRazorProtein,'SARS2');

% Shorten SARS2 proteins
dat.pep.A549.GN(dat.S2.A549) = extractAfter(dat.pep.A549.GN(dat.S2.A549),'|');
dat.pep.Vero.GN(dat.S2.Vero) = extractAfter(dat.pep.Vero.GN(dat.S2.Vero),'|');

dat.pep.A549.GN(dat.S2.A549) = extractBefore(dat.pep.A549.GN(dat.S2.A549),'|');
dat.pep.Vero.GN(dat.S2.Vero) = extractBefore(dat.pep.Vero.GN(dat.S2.Vero),'|');

% Now: Vero fastas need annotating from file: Gene names are contained
% within the string for A549
dat.pep.A549.GN(~dat.S2.A549) = extractAfter(dat.pep.A549.GN(~dat.S2.A549),'|');
dat.pep.A549.GN(~dat.S2.A549) = extractAfter(dat.pep.A549.GN(~dat.S2.A549),'|');
dat.pep.A549.GN(~dat.S2.A549) = extractBefore(dat.pep.A549.GN(~dat.S2.A549),'_');

%% Match Vero to GN
dat.pep.Vero.GN(~dat.S2.Vero) = extractAfter(dat.pep.Vero.GN(~dat.S2.Vero),'|');
dat.pep.Vero.GN(~dat.S2.Vero) = extractBefore(dat.pep.Vero.GN(~dat.S2.Vero),'|');

% Import Vero annotations
dat.VeroGN = readtable([path , 'vero_acc_to_gn_signaltransit.csv']);

% Match up accessions
[lia locb] = ismember(dat.pep.Vero.GN , dat.VeroGN.Entry);

% Replace accessions with Gene Names
dat.pep.Vero.GN(lia) = dat.VeroGN.GeneNames(locb(lia));

% Note PUR6 = PAICS, so replacing this manually for consistency throughout
% the manuscript
dat.pep.A549.GN{contains(dat.pep.A549.GN,'PUR6')} = 'PAICS'; 

%% Plot all data

figure
subplot(2,6,7:9) % Row 2, 1:2 (2,4,4)
scatter(log2(dat.mean.A549(~dat.hits.A549)) , -log10(dat.q2.A549(~dat.hits.A549)),'filled','k','MarkerFaceAlpha',0.2)
ylim([0 -pmax])
xlabel('Mean Log_2 fold-change (Infected over Mock)');
ylabel('-Log_1_0 Q-value');
hold on
% Testing in progress
scatter(log2(dat.mean.A549(dat.hits.A549 + dat.viral.A549 == 1)) , -log10(dat.q2.A549(dat.hits.A549 + dat.viral.A549 == 1)),'filled','r','MarkerFaceAlpha',0.4)

scatter(log2(dat.mean.A549(dat.viral.A549)) , -log10(dat.q2.A549(dat.viral.A549)),'filled','b','MarkerFaceAlpha',0.4)

%

line([-fccut -fccut],[0 -pmax],'Color','k','LineStyle','--')
line([fccut fccut],[0 -pmax],'Color','k','LineStyle','--')
line([-xcut xcut],[-log10(qcut) -log10(qcut)],'Color','k','LineStyle','--')
text(log2(dat.mean.A549(dat.viral.A549)) , -log10(dat.q2.A549(dat.viral.A549)),dat.pep.A549.GN(dat.viral.A549),'FontSize',10);
text(log2(dat.mean.A549(log2(dat.mean.A549) >= 2)) , -log10(dat.q2.A549(log2(dat.mean.A549) >= 2)),dat.pep.A549.GN(log2(dat.mean.A549) >= 2),'FontSize',10);
%set(gca,'Fontsize',14)
legend('N-termini','Cellular','Viral','location','northwest')
legend('boxoff')
title('neo-N-termini: A549-Ace2')
set(gca,'Fontsize',14)
hold off
xlim([-3,6])

writematrix([dat.mean.A549,dat.q2.A549],'Fig1E.csv');

subplot(2,6,1:2) % Row 1: 1 (2,4,5)
% Uncomment section to see individual datapoints
% for ii = 1:3
%     scatter(timecourseTimes , tc_r_A(ii,:),'MarkerFaceColor', 'k'); 
%     hold on; 
% end
% 
% for ii = 1:3
%     scatter(timecourseTimes , tc_r_V(ii,:),'MarkerFaceColor', 'r'); 
%     hold on; 
% end


errorbar(timecourseTimes, mean(tc_r_A),std(tc_r_A),'-ok','LineWidth',1,'MarkerFaceColor','k');
hold on
errorbar(timecourseTimes, mean(tc_r_V),std(tc_r_V),'-or','LineWidth',1,'MarkerFaceColor','r');

writematrix([tc_r_A;tc_r_V],'Fig1B.csv')

title('Viral RNA levels')
xlabel('Hours post-infection')
ylabel('RNA copies (PFU equivalent)/mL')
legend({'A549-Ace2','Vero E6'},'Location','southeast')
legend('boxoff')
xticks(tcTicks)
xticklabels(tcTicklabels)
xlim([-1,25])
set(gca,'yscale','log')
ylim([1e1,1e7])
set(gca,'Fontsize',14)
hold off

subplot(2,6,3:4) % Row 1, 2 (2,4,6)
errorbar(timecourseTimesP, mean(datP.Rs.A549,2),std(datP.Rs.A549'),'-ok','LineWidth',1,'MarkerFaceColor','k');
hold on
errorbar(timecourseTimesP, mean(datP.Rs.Vero,2),std(datP.Rs.Vero'),'-or','LineWidth',1,'MarkerFaceColor','r');

errorbar(timecourseTimesP, mean(datP.Rnsp8.A549,2),std(datP.Rnsp8.A549'),'--sk','LineWidth',1,'MarkerFaceColor','k');
errorbar(timecourseTimesP, mean(datP.Rnsp8.Vero,2),std(datP.Rnsp8.Vero'),'--sr','LineWidth',1,'MarkerFaceColor','r');

writematrix([datP.Rs.A549';datP.Rs.Vero';datP.Rnsp8.A549';datP.Rnsp8.Vero'],'Fig1C.csv')

title('Viral Protein levels')
xlabel('Hours post-infection')
ylabel('Log_2 Protein Intensity')
legend({'Spike: A549-Ace2','Spike: Vero E6','nsp8: A549-Ace2','nsp8: Vero E6'},'Location','southeast')
legend('boxoff')
xticks(tcTicks)
xticklabels(tcTicklabels)
xlim([-1,25])
%set(gca,'yscale','log2')
ylim([20,36])
set(gca,'Fontsize',14)
hold off

subplot(2,6,5:6) % Row 1: 3 (2,4,7)
% Uncomment section to see individual datapoints
% for ii = 1:3
%     scatter(timecourseTimes , tc_p_A(ii,:),'MarkerFaceColor', 'k'); 
%     hold on; 
% end
% 
% for ii = 1:3
%     scatter(timecourseTimes , tc_p_V(ii,:),'MarkerFaceColor', 'r'); 
%     hold on; 
% end

errorbar(timecourseTimes, mean(tc_p_A),std(tc_p_A),'-ok','LineWidth',1,'MarkerFaceColor','k');
hold on
errorbar(timecourseTimes, mean(tc_p_V),std(tc_p_V),'-or','LineWidth',1,'MarkerFaceColor','r');

writematrix([tc_p_A;tc_p_V],'Fig1D.csv')

title('Viral Titres')
xlabel('Hours post-infection')
ylabel('PFU/mL')
legend({'A549-Ace2','Vero E6'},'Location','southeast')
legend('boxoff')
xticks(tcTicks)
xticklabels(tcTicklabels)
xlim([-1,25])
set(gca,'yscale','log')
ylim([1e1,1e5])
set(gca,'Fontsize',14)
hold off
% Add Vero and A549 as colored text

subplot(2,6,10:12) % Row 2: 3:4 (2,4,8)
scatter(log2(dat.mean.Vero(~dat.hits.Vero)) , -log10(dat.q2.Vero(~dat.hits.Vero)),'filled','k','MarkerFaceAlpha',0.2)
ylim([0 -pmax])
xlabel('Mean Log_2 fold-change (Infected over Mock)');
ylabel('-Log_1_0 Q-value');
hold on
scatter(log2(dat.mean.Vero(dat.hits.Vero + dat.viral.Vero == 1)) , -log10(dat.q2.Vero(dat.hits.Vero + dat.viral.Vero == 1)),'filled','r','MarkerFaceAlpha',0.4)

scatter(log2(dat.mean.Vero(dat.viral.Vero)) , -log10(dat.q2.Vero(dat.viral.Vero)),'filled','b','MarkerFaceAlpha',0.4)

line([-fccut -fccut],[0 -pmax],'Color','k','LineStyle','--')
line([fccut fccut],[0 -pmax],'Color','k','LineStyle','--')
line([-xcut xcut],[-log10(qcut) -log10(qcut)],'Color','k','LineStyle','--')
text(log2(dat.mean.Vero(dat.viral.Vero)) , -log10(dat.q2.Vero(dat.viral.Vero)),dat.pep.Vero.GN(dat.viral.Vero),'FontSize',10);

text(log2(dat.mean.Vero(log2(dat.mean.Vero) >= 2)) , -log10(dat.q2.Vero(log2(dat.mean.Vero) >= 2)),dat.pep.Vero.GN(log2(dat.mean.Vero) >= 2),'FontSize',10);
set(gca,'Fontsize',14)
legend('N-termini','Cellular','Viral','location','northwest')
legend('boxoff')
title('neo-N-termini: Vero E6')
xlim([-3,6])
hold off


writematrix([dat.mean.Vero,dat.q2.Vero],'Fig1F.csv');


% Save figure
print([path , '/Figures/Fig_Timecourse.pdf'],'-dpdf');




