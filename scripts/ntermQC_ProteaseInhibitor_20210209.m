% QC data plots - Vero E6 infections following protease inhibitor treatment

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/'

% Goals:
% 1. Show enrichment % in Nterm-enriched fractions
% 2. Show breakdown of N-acetyl, pyroGly and TMT
% 3. Principal Component analysis showing clustering of replicates and
% separation of Mock/Infected

% The data is save as two separate figures: in /Figures/
% 1. PI_NtermQC.pdf
% 2. PI_PCA.pdf


% Notes: MaxQuant treats TMT modifications as fixed modifications, the
% unmodified peptides identified in the Var search file will not be picked
% up (or will have poor PEP and be filtered out).

dat = struct();
% Load the peptide files - Variable search
dat.var.Vero = readtable([path , 'tSARS2Inhib_var/evidence.txt']);

% Load the peptide files - Quant search
dat.pep.Vero = readtable([path , 'tSARS2Inhib/peptides.txt']);

% Load in peptide data for unenriched - Quant searched
dat.upep.Vero = readtable([path , 'tSARS2Inhib_Unenriched/peptides.txt']);

% Load in the TMTpro design matrices containing the randomised layouts

% VeroE6
dat.tmt = readmatrix([path , 'SARS2proteaseinhibitors_TMT_RandomisationStrategy.csv']);

%% Rearrange RI channels
% First remove unwanted channels
% Note data was searched as 16plex, though only 12 channels were used. The
% 4 unused channels are removed a this time.
dat.pep.Vero(:,[46,47,50,51]) = [];
dat.upep.Vero(:,[46,47,50,51]) = [];

% Rearrange remaining channels as per randomisation strategy
dat.pep.Vero(:,[44:55]) = dat.pep.Vero(:,dat.tmt + 43);
dat.upep.Vero(:,[44:55]) = dat.upep.Vero(:,dat.tmt + 43);

%% Data cleanup

samples = {'Vero'};
% Filter on PEP, Rev, Con
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

%% From variable search work out fraction of blocked N-termini

blocked = struct();

for ii = 1:numel(samples)
   [blocked.unq.(samples{ii}), blocked.ia.(samples{ii})] = unique(dat.var.(samples{ii}).ModifiedSequence); 
end

% All modified N-termini begin: '_('

for ii = 1:numel(samples)
   blocked.logMod.(samples{ii}) = contains(blocked.unq.(samples{ii}) , '_(');
end

% Calculate fraction

for ii = 1:numel(samples)
   blocked.frac.(samples{ii}) = sum(blocked.logMod.(samples{ii})) / numel(blocked.logMod.(samples{ii})); 
end

% Plot
figure
subplot(1,4,1)
%bar([blocked.frac.A549*100, blocked.frac.Vero*100],'k')
bar([blocked.frac.Vero*100],'k')
xticklabels({'Vero E6'})
xtickangle(45)
ylabel('% Blocked N-termini')
title('Blocked N-terminus enrichment')
ylim([0 100])
xlabel('Sample')

% Work out relative intensity of modified/unmodified
subplot(1,4,2)
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

writematrix(log10(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 0))),'Fig_S5d_unblocked.csv')
writematrix(log10(dat.var.Vero.Intensity(blocked.ia.Vero(blocked.logMod.Vero == 1))),'Fig_S5d_blocked.csv')


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

subplot(1,4,3)
bar([blocked.mods.Vero.N*100 , blocked.mods.Vero.P*100 , blocked.mods.Vero.T*100], 'k');
xticklabels({'Acetylation','Pyroglutamine','TMTpro'})
xtickangle(45)
ylabel('% of blocked N-termini')
title('Vero E6: N-termini')
xlabel('N-terminal modification')

subplot(1,4,4)
bar([...%numel(blocked.unqM.A549) * blocked.mods.A549.T,...
    numel(blocked.unqM.Vero) * blocked.mods.Vero.T],'k')
xticklabels({'A549-Ace2','Vero E6'});
xtickangle(45)
ylabel('# TMTpro-labelled N-termini')
title('# TMTpro-labelled N-termini')
xlabel('Sample')


print([path , '/Figures/PI_NtermQC.pdf'],'-dpdf');

%% Quantitative comparisons: PCA

% Convert table to array
for ii = 1:numel(samples)
    dat.mat.(samples{ii})  = table2array(dat.pep.(samples{ii})(: , 44:55));
    dat.umat.(samples{ii}) = table2array(dat.upep.(samples{ii})(:, 44:55));
end

% Normalise, Impute, Remove rows with all missing data
%
for ii = 1:numel(samples)
    % Convert 0 to NaN
    dat.mat.(samples{ii})(dat.mat.(samples{ii})   == 0) = NaN;
    dat.umat.(samples{ii})(dat.umat.(samples{ii}) == 0) = NaN;
    % Normalise by column median
   dat.mat.(samples{ii})  = dat.mat.(samples{ii}) ./ nanmedian(dat.mat.(samples{ii})); 
   dat.umat.(samples{ii}) = dat.umat.(samples{ii}) ./ nanmedian(dat.umat.(samples{ii}));
   % Renove rows with >= 50% missing data
   dat.mat.(samples{ii}) = dat.mat.(samples{ii})(sum(isnan(dat.mat.(samples{ii})) , 2) <= 6 , :);
   dat.umat.(samples{ii}) = dat.umat.(samples{ii})(sum(isnan(dat.umat.(samples{ii})) , 2) <= 6 , :);
   % KNN impute (k = 3)
   dat.mat.(samples{ii}) = knnimpute(dat.mat.(samples{ii}) , 3);
   dat.umat.(samples{ii}) = knnimpute(dat.umat.(samples{ii}) , 3);
   % Row normalise by row mean
   dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ mean(dat.mat.(samples{ii}) , 2);
   dat.umat.(samples{ii}) = dat.umat.(samples{ii}) ./ mean(dat.umat.(samples{ii}) , 2);
end

% Now prepare data for PCA
% Perform PCA on 10% most variable genes.
pcntV = 10;

for ii = 1:numel(samples)
    %Select % most variable genes
    dat.matvar.(samples{ii})  = var(dat.mat.(samples{ii}) , 0 , 2);
    dat.umatvar.(samples{ii}) = var(dat.umat.(samples{ii}) , 0 , 2);
    % Sort: return index
    [~ , dat.mvIA.(samples{ii})]  = sort(dat.matvar.(samples{ii}) , 'descend');
    [~ , dat.umvIA.(samples{ii})] = sort(dat.umatvar.(samples{ii}) , 'descend');
    % Grab index for top X% values
    dat.varIa.(samples{ii})  = dat.mvIA.(samples{ii})(1:round(numel(dat.mvIA.(samples{ii})) / pcntV),:);
    dat.uvarIa.(samples{ii}) = dat.umvIA.(samples{ii})(1:round(numel(dat.umvIA.(samples{ii})) / pcntV),:);
end

% PCA plots
% enriched
[dat.pca.Vero,~,~,~,dat.pcaV.Vero,~] = pca(log2(dat.mat.Vero(dat.varIa.Vero , :)))
% unenriched
[dat.upca.Vero,~,~,~,dat.upcaV.Vero,~] = pca(log2(dat.umat.Vero(dat.uvarIa.Vero , :)))

grpVar = categorical({'Mock','Mock','Mock','Infected','Infected','Infected'...
        'Infected + Camostat','Infected + Camostat','Infected + Camostat',...
        'Infected + Calpeptin','Infected + Calpeptin','Infected + Calpeptin'});


figure
subplot(1,2,1) % PCA on Vero, enriched data
pcA = 1;
pcB = 2;
hgB = gscatter(dat.pca.Vero(:,pcA) , dat.pca.Vero(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.pcaV.Vero(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.pcaV.Vero(pcB))),'%)']);
title('Vero E6 - Enriched')

writematrix([dat.pca.Vero(:,1) , dat.pca.Vero(:,2)],'Fig_S5a.csv')


subplot(1,2,2) % PCA on Vero, unenriched data
pcA = 1;
pcB = 2;
hgD = gscatter(dat.upca.Vero(:,pcA) , dat.upca.Vero(:,pcB) , grpVar);
xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.upcaV.Vero(pcA))),'%)']);
ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.upcaV.Vero(pcB))),'%)']);
title('Vero E6 - Unenriched')

writematrix([dat.upca.Vero(:,1) , dat.upca.Vero(:,2)],'Fig_S5b.csv')


print([path , '/Figures/PI_PCA.pdf'],'-dpdf');
