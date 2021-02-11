% SARS-CoV-2 protease inhibitor experiment

clear 
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
dat = struct();

% Load the data
dat = struct();

% Samples:
% 1) Mock
% 2) 24h Infected
% 3) 24h Infected + Camostat (inhibits TMPRSS)
% 4) 24h Infected + Calpeptin (Inhibits Calpains, Cathepsins)

% Drugs were added at 12hpi. Samples havested 24hpi.

% Note: 27C,28N,29C,30N omitted (though included in search).
% Means channels 3,4,7,8 removed and omitted from rearrangement
% calculations

% Enriched
dat.evi = readtable([path , 'tSARS2Inhib/evidence.txt']);
% evi
dat.pep = readtable([path , 'tSARS2Inhib/peptides.txt']);
% pep

% Enriched Var Search
% evi
dat.var = readtable([path , 'tSARS2Inhib_var/evidence.txt']);

% Unenriched
% protein
dat.prot = readtable([path , 'tSARS2Inhib_unenriched/proteinGroups.txt']);

% TMT labeling strategy
dat.tmt = readmatrix([path , 'SARS2proteaseinhibitors_TMT_RandomisationStrategy.csv']);

% Rearrange TMT channels, omitting unused channels and accounting for randomisation
% First remove unwanted channels (3,4,7,8)
% Note data was searched as 16plex, though only 12 channels were used. The
% 4 unused channels are removed a this time.
dat.evi(:,[54,55,58,59]) = [];
dat.pep(:,[46,47,50,51]) = [];
dat.prot(:,[21,22,25,26]) = [];

% Then reorder columns by the labeling strategy
dat.evi(:,[52:63]) = dat.evi(:,dat.tmt + 51);
dat.pep(:,[44:55]) = dat.pep(:,dat.tmt + 43);
dat.prot(:,[19:30]) = dat.prot(:,dat.tmt + 18);

%% Filter results
% remove Rev, Con and low PEP/FDR IDs.

dat.evi.PotentialContaminant = categorical(dat.evi.PotentialContaminant);
dat.evi.Reverse = categorical(dat.evi.Reverse);

dat.evi = dat.evi(dat.evi.PotentialContaminant ~= '+',:);
dat.evi = dat.evi(dat.evi.Reverse ~= '+',:);

dat.pep.PotentialContaminant = categorical(dat.pep.PotentialContaminant);
dat.pep.Reverse = categorical(dat.pep.Reverse);

dat.pep = dat.pep(dat.pep.PotentialContaminant ~= '+',:);
dat.pep = dat.pep(dat.pep.Reverse ~= '+',:);

dat.prot.PotentialContaminant = categorical(dat.prot.PotentialContaminant);
dat.prot.Reverse = categorical(dat.prot.Reverse);

dat.prot = dat.prot(dat.prot.PotentialContaminant ~= '+',:);
dat.prot = dat.prot(dat.prot.Reverse ~= '+',:);

%% First up: analysis of the protein-level data
% Goal: examine total viral protein levels, both as a direct result, but
% also for normalising cleavage site abundance to.

% 0 to NaN
dat.protM = table2array(dat.prot(:,[19:30]));
dat.protM(dat.protM == 0) = NaN;

% Median normalise columns
dat.protM = dat.protM ./ nanmedian(dat.protM , 1);

% Remove rows with mainly missing data (>= 6)
logNan = sum(isnan(dat.protM) , 2) >= 6;

dat.protM = dat.protM(~logNan , :);
dat.prot  = dat.prot(~logNan , :); 

% Knn impute missing data
dat.protM = knnimpute(dat.protM);

% Row normalise by dividing by the mean of the mock channels
dat.protM = dat.protM ./ mean(dat.protM(:,[1:3]) , 2);

% Generate new cell array containing viral gene names to annotate table
% with.
s2names = {'M','N','ORF9B','nsp3','ORF3A','ORF7A','pp1ab','S'};
    
%% Generate heatmap of the results
[b,idx] = sort(max(dat.protM([1:8],:),[],2),'ascend')

figure
imagesc(log2(dat.protM([idx],:)));
colormap(redblue(65));
c = colorbar;
c.Label.String = 'log_2 fold-change relative to Mock (24h post-infection)';
yticks(1:8);
yticklabels(s2names(idx));
xticks([2,5,8,11]);
xticklabels({'Mock','Infected','Infected + Camostat','Infected + Calpeptin'});
xtickangle(45)
title('Viral protein levels 24h post-infection');

% Stats on protein expression differences
% Note: Q values obtained given the small sample size are smaller than the
% p-values. As such the more conservative p-values will be used.
% Infected vs. Infected + Camostat
[~ , dat.p.InfVsInfCam , ~ , ~] = ttest2(log2(dat.protM(idx,[4:6])) , log2(dat.protM(idx,[7:9])),'Dim',2);

% Infected vs. Infected + Calpeptin
[~ , dat.p.InfVsInfCal , ~ , ~] = ttest2(log2(dat.protM(idx,[4:6])) , log2(dat.protM(idx,[10:12])),'Dim',2);

set(gca,'FontSize',12)

print([path , '/Figures/PI_Viral_proteinData.pdf'],'-dpdf');
%% Analysis part 2: Examine if viral neo-N-termini show regulation by the protease inhibitors.

% Cleanup was performed in part earlier: Filter on PEP 0.02
dat.pep = dat.pep(dat.pep.PEP <= 0.02,:);

% Extract TMT RI matrix
dat.pepM = table2array(dat.pep(:,[44:55]));

% 0 to NaN
dat.pepM(dat.pepM == 0) = NaN;

% Column normalise by median
dat.pepM = dat.pepM ./ nanmedian(dat.pepM , 1);

% Remove rows with >=50% NaN;
logpNaN = sum(isnan(dat.pepM) , 2) >= 6;

dat.pepM(logpNaN , :) = [];
dat.pep(logpNaN , :) = [];

% Knnimpute missing data
dat.pepM = knnimpute(dat.pepM);

% Row normalise by dividing by the mean of the mock channels
dat.pepM = dat.pepM ./ mean(dat.pepM(:,[1:3]) , 2);

% Select viral cleavage sites
dat.vpepM = dat.pepM(contains(dat.pep.Proteins , 'SARS2'),:);
dat.vpep = dat.pep(contains(dat.pep.Proteins , 'SARS2'),:);
%%
% Reorder by: max, then by viral protein
% Reorder by max
% Reorder by start position
[b , idx] = sort(dat.vpep.StartPosition,'ascend');

dat.vpepM = dat.vpepM(idx,:);
dat.vpep  = dat.vpep(idx,:);

clear b idx

% Reorder by viral protein
[b, idx] = sort(dat.vpep.Proteins);

dat.vpepM = dat.vpepM(idx,:);
dat.vpep  = dat.vpep(idx,:);

clear b idx

%% Visualise

figure
imagesc(log2(dat.vpepM)); colormap(redblue(65)); 
c = colorbar;
c.Label.String = 'Log_2 fold-change over Mock';

xticks([2,5,8,11]);
xticklabels({'Mock','Infected','Infected + Camostat','Infected + Calpeptin'});
xtickangle(45)

ylb = {'N','N','N','N','N','N','N','N','N','N','N','N','N','ORF7A','ORF9B','S','S','S'};

for ii = 1:numel(ylb)
    ylb{ii} = [ylb{ii},' (',num2str(dat.vpep.StartPosition(ii)),')'];
end

set(gca, 'YTick', [1:18], 'YTickLabel', ylb) 

title('Abundance of viral neo-N-termini following inhibitor treatment');
ylabel('Viral protein and neo-N-terminus start position')

set(gca, 'FontSize',12);

print([path , '/Figures/PI_Viral_Unormalised.pdf'],'-dpdf');
%% STATS on unnormalised
[~ , dat.p.PEPIvCam , ~ , ~] = ttest2(log2(dat.vpepM(:,[4:6])) , log2(dat.vpepM(:,[7:9])),'Dim',2);
dat.sig.PEPIvCam = dat.p.PEPIvCam <= 0.05

% WT vs KO
[~ , dat.p.PEPIvCal , ~ , ~] = ttest2(log2(dat.vpepM(:,[4:6])) , log2(dat.vpepM(:,[10:12])),'Dim',2);
dat.sig.PEPIvCal = dat.p.PEPIvCal <= 0.05

%% Now normalise to viral protein levels (per sample normalisation);
% N is row 7, ORF7A row 4, orf9b row 5 and S row 6 (using sorted data).
[b,idx] = sort(max(dat.protM([1:8],:),[],2),'ascend');

dat.normprotM = dat.protM(idx,:);

% Normalise N data (1 to 13)
dat.normvpepM = dat.vpepM;
dat.normvpepM(1:13,:) = dat.normvpepM(1:13,:) ./ dat.normprotM(7,:);

% Normalise ORF7A data (14)
dat.normvpepM(14,:) = dat.normvpepM(14,:) ./ dat.normprotM(4,:);

% Normalise ORF9B data (15)
dat.normvpepM(15,:) = dat.normvpepM(15,:) ./ dat.normprotM(5,:);

% Normalise S data (6)
dat.normvpepM(16:18,:) = dat.normvpepM(16:18,:) ./ dat.normprotM(6,:);

% Row normalise again
dat.normvpepM = dat.normvpepM(:,4:12) ./ mean(dat.normvpepM(:,4:12),2);

%% visualise
figure
imagesc(log2(dat.normvpepM)); colormap(redblue(65)); 
c = colorbar;
c.Label.String = 'Log_2 fold-change over Mock (normalised)';

xticks([2,5,8]);

xticklabels({'Infected','Infected + Camostat','Infected + Calpeptin'});
xtickangle(45)


ylb = {'N','N','N','N','N','N','N','N','N','N','N','N','N','ORF7A','ORF9B','S','S','S'};

for ii = 1:numel(ylb)
    ylb{ii} = [ylb{ii},' (',num2str(dat.vpep.StartPosition(ii)),')'];
end

caxis([-2.5 2.5])

set(gca, 'YTick', [1:18], 'YTickLabel', ylb) 

title('Abundance of viral neo-N-termini following inhibitor treatment');
ylabel('Viral protein (neo-N-terminus start position)')

set(gca,'FontSize',12)

print([path , '/Figures/PI_Viral_Normalised.pdf'],'-dpdf');
%% STATS on normalised

[~ , dat.p.nPEPIvCam , ~ , ~] = ttest2(log2(dat.normvpepM(:,[1:3])) , log2(dat.normvpepM(:,[4:6])),'Dim',2);
dat.sig.nPEPIvCam = dat.p.nPEPIvCam <= 0.05

% WT vs KO
[~ , dat.p.nPEPIvCal , ~ , ~] = ttest2(log2(dat.normvpepM(:,[1:3])) , log2(dat.normvpepM(:,[7:9])),'Dim',2);
dat.sig.nPEPIvCal = dat.p.nPEPIvCal <= 0.05
