% QC data plots

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Goals:
% 1. Show enrichment % in Nterm-enriched fractions
% 2. Show breakdown of N-acetyl, pyroGly and TMT
% 3. Principal Component analysis showing clustering of replicates and
% separation of Mock/Infected

% The data is save as two separate figures: in /Figures/
% 1. NtermQC.pdf
% 2. PCA.pdf


% Notes: MaxQuant treats TMT modifications as fixed modifications, the
% unmodified peptides identified in the Var search file will not be picked
% up (or will have poor PEP and be filtered out).

dat = struct();
% Load the peptide files - Variable search
dat.var.A549 = readtable([path , 'tA549_EnrichVar/evidence.txt']);
dat.var.Vero = readtable([path , 'tVero_EnrichVar/evidence.txt']);

% Load the peptide files - Quant search
dat.pep.A549 = readtable([path , 'tA549_Enrich/peptides.txt']);
dat.pep.Vero = readtable([path , 'tVero_Enrich/peptides.txt']);

% Load in peptide data for unenriched - Quant searched
dat.upep.A549 = readtable([path , 'tA549_TotalFrac/peptides.txt']);
dat.upep.Vero = readtable([path , 'tVero_TotalFrac/peptides.txt']);

% Load in the TMTpro design matrices containing the randomised layouts
% A549-ACE2
dat.tmt.A549 = readmatrix([path , 'SARS2a549tmtlabelling_20200507.csv']);
% VeroE6
dat.tmt.Vero = readmatrix([path , 'Verotmtlabelling_20200511.csv']);



% Rearrange RI channels
% Enriched
dat.pep.A549(:,44:59) = dat.pep.A549(: , dat.tmt.A549(2,:) + 43);
dat.pep.Vero(:,44:59) = dat.pep.Vero(: , dat.tmt.Vero(2,:) + 43);
% Unenriched
dat.upep.A549(:,44:59) = dat.upep.A549(: , dat.tmt.A549(2,:) + 43);
dat.upep.Vero(:,44:59) = dat.upep.Vero(: , dat.tmt.Vero(2,:) + 43);

%% Data cleanup

samples = {'A549','Vero'};
% Filter on PEP, Rev, Con
% Peptide PEP <= 0.02
% Reverse database hits and potential contaminants (from MaxQuant
% contaminants list) removed
for ii = 1:numel(samples)
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).PEP <= 0.02 , :);
    dat.pep.(samples{ii}).Reverse = categorical(dat.pep.(samples{ii}).Reverse);
    dat.pep.(samples{ii}).PotentialContaminant = categorical(dat.pep.(samples{ii}).PotentialContaminant);
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).Reverse ~= '+' , :);
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).PotentialContaminant ~= '+' , :);
    
    dat.upep.(samples{ii}) = dat.upep.(samples{ii})(dat.upep.(samples{ii}).PEP <= 0.02 , :);
    dat.upep.(samples{ii}).Reverse = categorical(dat.upep.(samples{ii}).Reverse);
    dat.upep.(samples{ii}).PotentialContaminant = categorical(dat.upep.(samples{ii}).PotentialContaminant);
    dat.upep.(samples{ii}) = dat.upep.(samples{ii})(dat.upep.(samples{ii}).Reverse ~= '+' , :);
    dat.upep.(samples{ii}) = dat.upep.(samples{ii})(dat.upep.(samples{ii}).PotentialContaminant ~= '+' , :);
    
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).PEP <= 0.02 , :);
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).PIF >= 0.7  , :);
    dat.var.(samples{ii}).Reverse = categorical(dat.var.(samples{ii}).Reverse);
    dat.var.(samples{ii}).PotentialContaminant = categorical(dat.var.(samples{ii}).PotentialContaminant);
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).Reverse ~= '+' , :);
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).PotentialContaminant ~= '+' , :);
end

%
%% From variable search work out fraction of blocked N-termini

blocked = struct();
% First identify unique peptide sequences (with modifications)
for ii = 1:numel(samples)
   [blocked.unq.(samples{ii}), blocked.ia.(samples{ii})] = unique(dat.var.(samples{ii}).ModifiedSequence); 
end

% All modified N-termini begin: '_('
% Identify all sequences with a modified N-terminus (As due to our search
% criteria, such sequences can only be TMT labeled, Acetylated, or
% Pyroglutamate.
for ii = 1:numel(samples)
   blocked.logMod.(samples{ii}) = contains(blocked.unq.(samples{ii}) , '_(');
end

% Calculate fraction of total peptides that are blocked
for ii = 1:numel(samples)
   blocked.frac.(samples{ii}) = sum(blocked.logMod.(samples{ii})) / numel(blocked.logMod.(samples{ii})); 
end

% Plot
figure
subplot(2,3,1)
bar([blocked.frac.A549*100, blocked.frac.Vero*100],'k')
xticklabels({'A549-Ace2','Vero E6'})
xtickangle(45)
ylabel('% Blocked N-termini')
title('Blocked N-terminus enrichment')
xlabel('Sample')

% Work out relative intensity of modified/unmodified
% If enrichment has worked, while unenriched will likely be identified due
% to the bRP fractionation to enhance coverage, they should be less
% intense.
subplot(2,3,2)
scatter(ones(numel(dat.var.A549.Intensity(blocked.ia.A549(blocked.logMod.A549 == 0))),1),...
    log10(dat.var.A549.Intensity(blocked.ia.A549(blocked.logMod.A549 == 0))),'k','filled','MarkerFaceAlpha',0.05,'jitter','on','jitteramount',0.25)
hold on
scatter(ones(numel(dat.var.A549.Intensity(blocked.ia.A549(blocked.logMod.A549 == 1))),1)+1,...
    log10(dat.var.A549.Intensity(blocked.ia.A549(blocked.logMod.A549 == 1))),'k','filled','MarkerFaceAlpha',0.05,'jitter','on','jitteramount',0.25)

h = boxplot(log10(dat.var.A549.Intensity(blocked.ia.A549)),blocked.logMod.A549);
set(h,{'linew'},{3})
xticklabels({'Unblocked','Blocked'})
xtickangle(45)
ylabel('Log_1_0 peptide intensity')
title('A549-Ace2: N-termini')
xlabel('N-terminus')
ylim([6 11])

subplot(2,3,3)
scatter(ones(numel(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 0))),1),...
    log10(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 0))),'k','filled','MarkerFaceAlpha',0.05,'jitter','on','jitteramount',0.25)
hold on
scatter(ones(numel(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 1))),1)+1,...
    log10(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 1))),'k','filled','MarkerFaceAlpha',0.05,'jitter','on','jitteramount',0.25)

h = boxplot(log10(dat.var.Vero.Intensity(blocked.ia.Vero)),blocked.logMod.Vero);
set(h,{'linew'},{3})
xticklabels({'Unblocked','Blocked'})
xtickangle(45)
ylabel('Log_1_0 peptide intensity')
title('Vero E6: N-termini')
xlabel('N-terminus')
ylim([6 11])

% Now want distribution of N-terminal Mods.

% So N-terminal mods will all begin: '_(A , _(G or _(T for N-terminal Ac,
% Pyroglu and TMT respectively

% First select modified only
for ii = 1:numel(samples)
   blocked.unqM.(samples{ii}) = blocked.unq.(samples{ii})(blocked.logMod.(samples{ii}) == 1); 
end

% Then grab numbers of each
for ii = 1:numel(samples)
    blocked.mods.(samples{ii}).N = sum(contains(blocked.unqM.(samples{ii}) , '_(A')) / numel(blocked.unqM.(samples{ii}));
    blocked.mods.(samples{ii}).P = sum(contains(blocked.unqM.(samples{ii}) , '_(G')) / numel(blocked.unqM.(samples{ii}));
    blocked.mods.(samples{ii}).T = sum(contains(blocked.unqM.(samples{ii}) , '_(T')) / numel(blocked.unqM.(samples{ii}));
end

%figure
subplot(2,3,4)
bar([blocked.mods.A549.N*100 , blocked.mods.A549.P*100 , blocked.mods.A549.T*100], 'k');
xticklabels({'Acetylation','Pyroglutamine','TMTpro'})
xtickangle(45)
ylabel('% of blocked N-termini')
title('A549-Ace2: N-termini')
xlabel('N-terminal modification')

subplot(2,3,5)
bar([blocked.mods.Vero.N*100 , blocked.mods.Vero.P*100 , blocked.mods.Vero.T*100], 'k');
xticklabels({'Acetylation','Pyroglutamine','TMTpro'})
xtickangle(45)
ylabel('% of blocked N-termini')
title('Vero E6: N-termini')
xlabel('N-terminal modification')

subplot(2,3,6)
bar([numel(blocked.unqM.A549) * blocked.mods.A549.T,...
    numel(blocked.unqM.Vero) * blocked.mods.Vero.T],'k')
xticklabels({'A549-Ace2','Vero E6'});
xtickangle(45)
ylabel('# TMTpro-labelled N-termini')
title('# TMTpro-labelled N-termini')
xlabel('Sample')

% Save figure as .pdf
print([path , '/Figures/NtermQC.pdf'],'-dpdf');

%% Quantitative comparisons: PCA

% Convert table to array
for ii = 1:numel(samples)
    dat.mat.(samples{ii})  = table2array(dat.pep.(samples{ii})(: , 44:59));
    dat.umat.(samples{ii}) = table2array(dat.upep.(samples{ii})(:, 44:59));
end

% Normalise, Impute, Remove rows with all missing data

for ii = 1:numel(samples)
    % Convert 0 to NaN
    dat.mat.(samples{ii})(dat.mat.(samples{ii})   == 0) = NaN;
    dat.umat.(samples{ii})(dat.umat.(samples{ii}) == 0) = NaN;
    % Normalise by column median
   dat.mat.(samples{ii})  = dat.mat.(samples{ii}) ./ nanmedian(dat.mat.(samples{ii})); 
   dat.umat.(samples{ii}) = dat.umat.(samples{ii}) ./ nanmedian(dat.umat.(samples{ii}));
   % Renove rows with >= 50% missing data
   dat.mat.(samples{ii}) = dat.mat.(samples{ii})(sum(isnan(dat.mat.(samples{ii})) , 2) <= 8 , :);
   dat.umat.(samples{ii}) = dat.umat.(samples{ii})(sum(isnan(dat.umat.(samples{ii})) , 2) <= 8 , :);
   % KNN impute (k = 3)
   dat.mat.(samples{ii}) = knnimpute(dat.mat.(samples{ii}) , 3);
   dat.umat.(samples{ii}) = knnimpute(dat.umat.(samples{ii}) , 3);
   % Row normalise by row mean
   dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ mean(dat.mat.(samples{ii}) , 2);
   dat.umat.(samples{ii}) = dat.umat.(samples{ii}) ./ mean(dat.umat.(samples{ii}) , 2);
end

% Now prepare data for PCA:
% Select 10% most variable peptides.
pcntV = 10;

for ii = 1:numel(samples)
    %Select calculate variance for each peptide
    dat.matvar.(samples{ii})  = var(dat.mat.(samples{ii}) , 0 , 2);
    dat.umatvar.(samples{ii}) = var(dat.umat.(samples{ii}) , 0 , 2);
    % Sort on variance: return index
    [~ , dat.mvIA.(samples{ii})]  = sort(dat.matvar.(samples{ii}) , 'descend');
    [~ , dat.umvIA.(samples{ii})] = sort(dat.umatvar.(samples{ii}) , 'descend');
    % Grab index for top X% values (defined by pcntV)
    dat.varIa.(samples{ii})  = dat.mvIA.(samples{ii})(1:round(numel(dat.mvIA.(samples{ii})) / pcntV),:);
    dat.uvarIa.(samples{ii}) = dat.umvIA.(samples{ii})(1:round(numel(dat.umvIA.(samples{ii})) / pcntV),:);
end

% PCA plots
% enriched
[dat.pca.A549,~,~,~,dat.pcaV.A549,~] = pca(log2(dat.mat.A549(dat.varIa.A549 , :)))
[dat.pca.Vero,~,~,~,dat.pcaV.Vero,~] = pca(log2(dat.mat.Vero(dat.varIa.Vero , :)))
% unenriched
[dat.upca.A549,~,~,~,dat.upcaV.A549,~] = pca(log2(dat.umat.A549(dat.uvarIa.A549 , :)))
[dat.upca.Vero,~,~,~,dat.upcaV.Vero,~] = pca(log2(dat.umat.Vero(dat.uvarIa.Vero , :)))

grpVar = categorical({'0M','0','0','0','6','6','6','12','12','12','24','24','24','24M','24M','24M'});

figure
subplot(1,4,1) % PCA on A549, enriched data
pcA = 1;
pcB = 2;
hgA = gscatter(dat.pca.A549(:,pcA) , dat.pca.A549(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.pcaV.A549(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.pcaV.A549(pcB))),'%)']);
title('A549-Ace2 - Enriched')
h6 = hgA(6)
h6.Color = 'k';

subplot(1,4,2) % PCA on Vero, enriched data
pcA = 1;
pcB = 2;
hgB = gscatter(dat.pca.Vero(:,pcA) , dat.pca.Vero(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.pcaV.Vero(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.pcaV.Vero(pcB))),'%)']);
title('Vero E6 - Enriched')
h6 = hgB(6)
h6.Color = 'k';

subplot(1,4,3) % PCA on A549, unenriched data
pcA = 3; % Note: PCA 3 in this specific case gave better separation than PCA 1.
pcB = 2;
hgC = gscatter(dat.upca.A549(:,pcA) , dat.upca.A549(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.upcaV.A549(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.upcaV.A549(pcB))),'%)']);
title('A549-Ace2 - Unenriched')
h6 = hgC(6)
h6.Color = 'k';

subplot(1,4,4) % PCA on Vero, unenriched data
pcA = 1;
pcB = 2;
hgD = gscatter(dat.upca.Vero(:,pcA) , dat.upca.Vero(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.upcaV.Vero(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.upcaV.Vero(pcB))),'%)']);
title('Vero E6 - Unenriched')
h6 = hgD(6)
h6.Color = 'k';

% Save figure
print([path , '/Figures/PCA.pdf'],'-dpdf');
