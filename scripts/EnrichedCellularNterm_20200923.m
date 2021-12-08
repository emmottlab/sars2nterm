%% Figure 3 and STables - Cellular factors modulated by infection


clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
dat = struct();


% Load the Quant search evidence and peptide files
dat.evi.A549 = readtable([path , 'tA549_Enrich/evidence.txt']);1
dat.evi.Vero = readtable([path , 'tVero_Enrich/evidence.txt']);2

dat.pep.A549 = readtable([path , 'tA549_Enrich/peptides.txt']);3
dat.pep.Vero = readtable([path , 'tVero_Enrich/peptides.txt']);4

% Load the Var search evidence files
dat.var.A549 = readtable([path , 'tA549_EnrichVar/evidence.txt']);5
dat.var.Vero = readtable([path , 'tVero_EnrichVar/evidence.txt']);6

samples = {'A549','Vero'};

%% Import gene name/signal peptide data from uniprot
gn = struct();

gn.human = readtable([path , 'human_acc_to_gn_signaltransit.csv']);
gn.vero  = readtable([path , 'vero_acc_to_gn_signaltransit.csv']);

% Load in the TMTpro design matrices containing the randomised layouts
% A549-ACE2
dat.tmt.A549 = readmatrix([path , 'SARS2a549tmtlabelling_20200507.csv']);
% VeroE6
dat.tmt.Vero = readmatrix([path , 'Verotmtlabelling_20200511.csv']);

% Rearrange the Reporter channels in dat.evi to obtain the correct order of
% Channels:
% Channel layout from 1-16 after reordering is:
% 0h  Mock
% 0h  A
% 0h  B
% 0h  C
% 6h  A
% 6h  B
% 6h  C
% 12h A
% 12h B
% 12h C
% 24h A
% 24h B
% 24h C
% 24h Mock A
% 24h Mock B
% 24h Mock C

% Reorder the corrected RI intensity channels

dat.pep.A549(:,44:59) = dat.pep.A549(: , dat.tmt.A549(2,:) + 43);
dat.pep.Vero(:,44:59) = dat.pep.Vero(: , dat.tmt.Vero(2,:) + 43);


% Backup temp
backup = dat;
%% Data cleanup
dat = backup;

% Remove Reverse hits, Contaminants and low PEP identifications
for ii = 1:numel(samples)
    % Remove PEP <= 0.02
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).PEP <= 0.02 , :);
    % Remove Reverse hits
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(categorical(dat.pep.(samples{ii}).Reverse) ~= '+' , :);
    % Remove potential contaminants
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(categorical(dat.pep.(samples{ii}).PotentialContaminant) ~= '+' , :);
end
%
% Now Convert the TMT RI data to matrix form

for ii = 1:numel(samples)
   % Convert table columns to a matrix
    dat.mat.(samples{ii}) = table2array(dat.pep.(samples{ii})(: , [44:59]));
   
   % Convert 0 to NaN;
   dat.mat.(samples{ii})(dat.mat.(samples{ii}) == 0) = NaN; 
end

% Remove unquantified hits, column normalise, impute and row normalise.
for ii = 1:numel(samples)
   % identify rows with >75% NaN
   allNaN = sum(isnan(dat.mat.(samples{ii})) , 2) >13;
   dat.mat.(samples{ii}) = dat.mat.(samples{ii})(~allNaN , :);
   dat.pep.(samples{ii}) = dat.pep.(samples{ii})(~allNaN , :);
   
   clear allNaN
   
   % Median normalise
    dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ nanmedian(dat.mat.(samples{ii}));

    % Knn impute missing data
    dat.mat.(samples{ii}) = knnimpute(dat.mat.(samples{ii}));

    % Row normalise by mean (whole dataset) 
    dat.mat.(samples{ii}) = dat.mat.(samples{ii}) ./ mean(dat.mat.(samples{ii}) , 2);
 
    % Generate matrix with 24h data only
    dat.h24.(samples{ii}) = dat.mat.(samples{ii})(:,[end-5: end]);
   
end

backup2 = dat;

%%
dat = backup2;

 for ii = 1:numel(samples)
    dat.h24.(samples{ii}) = dat.h24.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
    dat.mat.(samples{ii}) = dat.mat.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
    dat.pep.(samples{ii}) = dat.pep.(samples{ii})(dat.pep.(samples{ii}).StartPosition > 2 , :);
 end
 
 
%% KS testing of proteins matching LxGG or (STVP)xLQ, DEVD?
 
 % Extract the P4 to P1 sequence from the Nterm cleavage window
 for ii = 1:numel(samples) 
     
     for jj = 1:numel(dat.pep.(samples{ii}).N_termCleavageWindow)
     
         dat.pep.(samples{ii}).p4p1{jj} = dat.pep.(samples{ii}).N_termCleavageWindow{jj}(12:15);
     
     end
     
 end
 
 % Regexp
 dat.reg.Nsp5 = '[A|S|T|V|P].LQ'; % Note added 'A'
 dat.reg.Nsp3 = 'L.GG';
 
for ii = 1:numel(samples)
    regNsp3  = regexp(dat.pep.(samples{ii}).p4p1 , dat.reg.Nsp3);
    regNsp5  = regexp(dat.pep.(samples{ii}).p4p1 , dat.reg.Nsp5);
    
    dat.viralProMatch.nsp3.(samples{ii}) = ~cellfun(@isempty , regNsp3);
    dat.viralProMatch.nsp5.(samples{ii}) = ~cellfun(@isempty , regNsp5);
    
    clear regNsp3 regNsp5
end
 
 %% Plot distribution of results

yl = 3.5; % set y axis limit

figure
t = tiledlayout(2,2,'TileSpacing','compact');

% A549 nsp5
 dat.vmat.A549 = [dat.h24.A549(:,1);dat.h24.A549(:,2);dat.h24.A549(:,3)];
[B, I] = sort(dat.vmat.A549, 'descend');
 
nexttile
 b = bar(log2(B),'FaceColor','k');
 b.FaceColor = 'flat';
 b.FaceAlpha = 0.3;

 dat.vpv.nsp5.A549 = [dat.viralProMatch.nsp5.A549;dat.viralProMatch.nsp5.A549;dat.viralProMatch.nsp5.A549];
 
  test = 1:numel(dat.vmat.A549);

  test = test(dat.vpv.nsp5.A549(I));
hold on
 
  scatter(test,log2(B(dat.vpv.nsp5.A549(I))),'filled','MarkerFaceAlpha',0.5); % Add scatter plot indicating peptides matching protease consensus
  

 [h,p,k] =  kstest2(test, 1:numel(dat.vmat.A549)) % KS test (two-tailed)
 hold on
 text(numel(dat.vpv.nsp5.A549)/3 , 1.5,['p = ',num2str(p)],'FontSize',14) % Add KS p-value to plot
 title('A549-Ace2: (A|P|S|T|V)xLQ enrichment')
 ylim([-yl yl])
 hold off
 
 writematrix(log2(B),'FigS13_Aall.csv')
 writematrix([test',log2(B(dat.vpv.nsp5.A549(I)))],'FigS13_Asig.csv')
 
 % A549 nsp3
 dat.vmat.A549 = [dat.h24.A549(:,1);dat.h24.A549(:,2);dat.h24.A549(:,3)];
[B, I] = sort(dat.vmat.A549, 'descend');
 
nexttile
 b = bar(log2(B),'FaceColor','k');
 b.FaceColor = 'flat';
 b.FaceAlpha = 0.3;

 dat.vpv.nsp3.A549 = [dat.viralProMatch.nsp3.A549;dat.viralProMatch.nsp3.A549;dat.viralProMatch.nsp3.A549];
 
  test = 1:numel(dat.vmat.A549);
  test = test(dat.vpv.nsp3.A549(I));
hold on
 
scatter(test,log2(B(dat.vpv.nsp3.A549(I))),'filled','MarkerFaceAlpha',0.5); % Add scatter plot indicating peptides matching protease consensus
  
 [h,p,k] =  kstest2(test, 1:numel(dat.vmat.A549)) % KS test (two-tailed)
 hold on
 text(numel(dat.vpv.nsp3.A549)/3 , 1.5,['p = ',num2str(p)],'FontSize',14) % Add KS p-value to plot
 title('A549-Ace2: LxGG enrichment')
 ylim([-yl yl])
 hold off
 
  writematrix(log2(B),'FigS13_Ball.csv')
 writematrix([test',log2(B(dat.vpv.nsp3.A549(I)))],'FigS13_Bsig.csv')
 %
 dat.vmat.Vero = [dat.h24.Vero(:,1);dat.h24.Vero(:,2);dat.h24.Vero(:,3)];
[B, I] = sort(dat.vmat.Vero, 'descend');
 
nexttile
 b = bar(log2(B),'FaceColor','k');
 b.FaceColor = 'flat';
 b.FaceAlpha = 0.3;

 dat.vpv.nsp5.Vero = [dat.viralProMatch.nsp5.Vero;dat.viralProMatch.nsp5.Vero;dat.viralProMatch.nsp5.Vero];
 
  test = 1:numel(dat.vmat.Vero);
  test = test(dat.vpv.nsp5.Vero(I));
  
hold on
   scatter(test,log2(B(dat.vpv.nsp5.Vero(I))),'filled','MarkerFaceAlpha',0.5); % Add scatter plot indicating peptides matching protease consensus
  
 [h,p,k] =  kstest2(test, 1:numel(dat.vmat.Vero)) % KS test (two-tailed)
 hold on
 text(numel(dat.vpv.nsp5.Vero)/3 , 1.5,['p = ',num2str(p)],'FontSize',14) % Add KS p-value to plot
 title('Vero E6: (A|P|S|T|V)xLQ enrichment')
 ylim([-yl yl])
 hold off
  
   writematrix(log2(B),'FigS13_Call.csv')
 writematrix([test',log2(B(dat.vpv.nsp5.Vero(I)))],'FigS13_Csig.csv')
 
 % Vero nsp3
 dat.vmat.Vero = [dat.h24.Vero(:,1);dat.h24.Vero(:,2);dat.h24.Vero(:,3)];
[B, I] = sort(dat.vmat.Vero, 'descend');
 
nexttile
 b = bar(log2(B),'FaceColor','k');
 b.FaceColor = 'flat';
 b.FaceAlpha = 0.3;

 dat.vpv.nsp3.Vero = [dat.viralProMatch.nsp3.Vero;dat.viralProMatch.nsp3.Vero;dat.viralProMatch.nsp3.Vero];
 
  test = 1:numel(dat.vmat.Vero);
  test = test(dat.vpv.nsp3.Vero(I));
hold on
  scatter(test,log2(B(dat.vpv.nsp3.Vero(I))),'filled','MarkerFaceAlpha',0.5); % Add scatter plot indicating peptides matching protease consensus
 
 [h,p,k] =  kstest2(test, 1:numel(dat.vmat.Vero)); % KS test (two-tailed)
 hold on
 text(numel(dat.vpv.nsp3.Vero)/3 , 1.5,['p = ',num2str(p)],'FontSize',14) % Add KS p-value to plot

 title('Vero E6: LxGG enrichment')
 ylim([-yl yl])
 hold off

 xlabel(t,'Neo-N-terminal peptides');
 ylabel(t,'Log_2 24h Infected / Mock')
 
   writematrix(log2(B),'FigS13_Dall.csv')
 writematrix([test',log2(B(dat.vpv.nsp3.Vero(I)))],'FigS13_Dsig.csv')
 
 print([path , '/Figures/KStest_ConsensusNterm.pdf'],'-dpdf');
 
%% Heat maps - Vero and A549

% Note peptide hits manually curated from those showing t-test significance
% at 24h and matching or close to matching nsp3/5 protease consensus
% sequences and the subset showing significant fold-change)


a549seqs = {'ASQDENFGNTTPR',... % NUP107
    'ASSAASSASPVSR',... % XRCC1
    'FNSSDTVTSPQR',... % SRC - LxGG
    'SKDQITAGNAAR',... % PAICS
    'SSVVATSKER',... % PNN
    'VTTTANKVGR'}; % WNK1
    %'TQGPPDYPR',... % MRPL49  - removed as matches degrabase


veroseqs = {'ALEVLPVAPPPEPR',... % ATAD2
    'APAPWHGEGTSPQLR',... % BST1 - AEGG
    'AQVECSHSSQQR',... % GOLGA3
    'ASSAASSASPVSR',... % XRCC1
    'FQESDDADEDYGR',... % NUCKS1
    'GAAGAGGGGSGAGGGSGGSGGR',... % KLHDC10 - RRGG
    'SFGTEEPAYSTR',... KAT7
    'TSPSPKAGAATGR',... % ATP51B
    'TTSSSITLR'}; % MYLK
    
% T test and multiple hypothesis correction
% Perform t-test for significance
[~ , dat.p.A549 , ~ , ~] = ttest2(log2(dat.h24.A549(:,[1:3])) , log2(dat.h24.A549(:,[4:6])),'Dim',2);
[~ , dat.p.Vero , ~ , ~] = ttest2(log2(dat.h24.Vero(:,[1:3])) , log2(dat.h24.Vero(:,[4:6])),'Dim',2);

%Correct for multiple hypothesis testing
[dat.fdr.A549 , dat.q.A549] = mafdr(dat.p.A549);
[dat.fdr.Vero , dat.q.Vero] = mafdr(dat.p.Vero);

%% Select those samples which made the q-value cutoff
qcut = 0.05; % Q value cutoff

% Select based on qvalue cutoff
dat.matq.A549 = dat.mat.A549(dat.q.A549 <= qcut , :);
dat.matq.Vero = dat.mat.Vero(dat.q.Vero <= qcut , :);

dat.pepq.A549 = dat.pep.A549(dat.q.A549 <= qcut , :);
dat.pepq.Vero = dat.pep.Vero(dat.q.Vero <= qcut , :);

%% Now need to select based on sequences of interest

% Ismember gives only exact, not partial matches
 dat.heatmap.A549 = dat.matq.A549(ismember(dat.pepq.A549.Sequence,a549seqs),:);
 dat.heatmap.Vero = dat.matq.Vero(ismember(dat.pepq.Vero.Sequence,veroseqs),:);

 dat.pepH.A549 = dat.pepq.A549(ismember(dat.pepq.A549.Sequence,a549seqs),:);
 dat.pepH.Vero = dat.pepq.Vero(ismember(dat.pepq.Vero.Sequence,veroseqs),:);



%% Annotate pepH with Gene Names
dat.pepH.A549.Acc = extractBetween(dat.pepH.A549.LeadingRazorProtein,4,9);
dat.pepH.Vero.Acc = extractBetween(dat.pepH.Vero.LeadingRazorProtein,4,13);

[lia, locb] = ismember(dat.pepH.A549.Acc , gn.human.Entry);
dat.pepH.A549.GN = gn.human.GeneNames(locb);
clear lia locb

[lia, locb] = ismember(dat.pepH.Vero.Acc , gn.vero.Entry);
dat.pepH.Vero.GN = gn.vero.GeneNames(locb);
clear lia locb

% Note: GN contains multiple gene names. Keep only the first
% Need to do in a loop so that it only does for those which contain ' '
lia = contains(dat.pepH.A549.GN , ' ');
dat.pepH.A549.GN(lia) = extractBefore(dat.pepH.A549.GN(lia) , ' ');
clear lia
lia = contains(dat.pepH.Vero.GN , ' ');
dat.pepH.Vero.GN(lia) = extractBefore(dat.pepH.Vero.GN(lia) , ' ');
clear lia

% Also want to sort on variance.
dat.hmV.A549 = mean(dat.heatmap.A549(:,11:13),2);
dat.hmV.Vero = mean(dat.heatmap.Vero(:,11:13),2);

[~ , dat.vI.A549] = sort(dat.hmV.A549,'descend');
[~ , dat.vI.Vero] = sort(dat.hmV.Vero,'descend');

%% Want to Add indication of consensus sequence and match to it.
a549yticks = dat.pepH.A549.GN;
veroyticks = dat.pepH.Vero.GN;

% Note: these will be incorrect if the identify/number of selected hits
% changes (shouldn't happen due to manual selection)

% Blue = nsp5 consensus match, red = nsp3 consensus match
a549consensus = {': \color{magenta}V\color{black}l\color{magenta}LQ';...
                 ': \color{magenta}A\color{black}t\color{magenta}LQ';...
                 ': \color{darkGreen}L\color{black}f\color{darkGreen}GG';...
                 ': \color{magenta}V\color{black}l\color{magenta}LQ';...
                 ': \color{magenta}P\color{black}a\color{magenta}LQ';...
                 ': \color{black}QrF\color{magenta}Q'};

veroconsensus = {': \color{black}Ev\color{magenta}LQ';...
                 ': \color{black}Ae\color{darkGreen}GG';...
                 ': \color{magenta}T\color{black}k\color{magenta}LQ';...
                 ': \color{magenta}A\color{black}t\color{magenta}LQ';...
                 ': \color{black}DyS\color{magenta}Q';...
                 ': \color{black}Rr\color{darkGreen}GG';...
                 ': \color{black}Rn\color{magenta}LQ';...
                 ': \color{black}YaA\color{magenta}Q';...
                 ': \color{magenta}P\color{black}v\color{magenta}LQ'};
             
% Now concatenate the Ytick labels so contain gene name, position, and P4 to P1 colored by match to consensus             
for ii = 1:numel(a549yticks)
   a549yticks{ii} = [a549yticks{ii} ,' (',dat.pepH.A549.FirstAminoAcid{ii},num2str(dat.pepH.A549.StartPosition(ii)),')', a549consensus{ii}]; 
end

for ii = 1:numel(veroyticks)
   veroyticks{ii} = [veroyticks{ii} ,' (',dat.pepH.Vero.FirstAminoAcid{ii},num2str(dat.pepH.Vero.StartPosition(ii)),')', veroconsensus{ii}]; 
end
           

figure
t = tiledlayout(2,1)
nexttile
imagesc(log2(dat.heatmap.A549(dat.vI.A549,:))); colormap( redbluecmap( 64 ) ); %rb;
c = colorbar('Ticks',[-2 , -1 , 0 , 1 , 2],'TickLabels',{'\leq-2','-1','0','1','\geq2'})
c.Label.String = 'Log_2 fold-change';
caxis([-2 2])
yticks([1:1:numel(veroyticks)])
yticklabels(a549yticks(dat.vI.A549));
xticks([])
title('A549-Ace2 cellular neo-N-termini')

nexttile
imagesc(log2(dat.heatmap.Vero(dat.vI.Vero,:))); colormap( redbluecmap( 64 ) );%rb;
c = colorbar('Ticks',[-2 , -1 , 0 , 1 , 2],'TickLabels',{'\leq-2','-1','0','1','\geq2'})
c.Label.String = 'Log_2 fold-change';
caxis([-2 2])
yticks([1:1:numel(veroyticks)])
yticklabels(veroyticks(dat.vI.Vero));
xticks([1,3,6,9,12,15])
xticklabels({'0M','0h','6h','12h','24h','24M'})
ylabel(t,'Gene name of cleaved protein')
xlabel(t,'Hours post-infection')
title('Vero E6 cellular neo-N-termini')

% Save figure
print([path , '/Figures/Fig_heatmaps.pdf'],'-dpdf');

%% Export full tables of N-termimi Quantification

% Get the accession number
dat.pep.A549.Acc = extractBetween(dat.pep.A549.LeadingRazorProtein,4,9);
dat.pep.Vero.Acc = extractBetween(dat.pep.Vero.LeadingRazorProtein,4,13);

% Match accession to gene name
[lia, locb] = ismember(dat.pep.A549.Acc , gn.human.Entry);
for ii = 1:numel(dat.pep.A549(:,1))
    dat.pep.A549.GN(ii) = {''};
end
dat.pep.A549.GN(lia) = gn.human.GeneNames(locb(lia));
clear lia locb

[lia, locb] = ismember(dat.pep.Vero.Acc , gn.vero.Entry);
for ii = 1:numel(dat.pep.Vero(:,1))
    dat.pep.Vero.GN(ii) = {''};
end
dat.pep.Vero.GN(lia) = gn.vero.GeneNames(locb(lia));
clear lia locb

% Note: GN contains multiple gene names. Keep only the first
% Need to do in a loop so that it only does for those which contain ' '
lia = contains(dat.pep.A549.GN , ' ');
dat.pep.A549.GN(lia) = extractBefore(dat.pep.A549.GN(lia) , ' ');
clear lia
lia = contains(dat.pep.Vero.GN , ' ');
dat.pep.Vero.GN(lia) = extractBefore(dat.pep.Vero.GN(lia) , ' ');
clear lia

% Generate new tables to hold the output.
TabA = table();
TabV = table();

% Export only the indicated columns
TabA = dat.pep.A549(:,[1,2,103,35,36,105,37,38,42,43,92,99,100]);

TabV = dat.pep.Vero(:,[1,2,103,35,36,105,37,38,42,43,92,99,100]);

TabA = [TabA , array2table(log2(dat.mat.A549))];
TabV = [TabV , array2table(log2(dat.mat.Vero))];

% Rename the TMT columns
timepoints = {'TMT_0M','TMT_0hA','TMT_0hB','TMT_0hC','TMT_6hA','TMT_6hB','TMT_6hC','TMT_12hA','TMT_12hB','TMT_12hC','TMT_24hA','TMT_24hB','TMT_24hC','TMT_24MA','TMT_24MB','TMT_24MC'};
TabA.Properties.VariableNames(14:29) = timepoints; %13:28
TabV.Properties.VariableNames(14:29) = timepoints; %13:28

% Write the output tables
writetable(TabA,[path , 'Tables/A549_Nterm_Quant.csv']);
writetable(TabV,[path , 'Tables/Vero_Nterm_Quant.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align Vero and A549 data

[c , ia, iv] = intersect(TabA.Sequence , TabV.Sequence);


% c = common peptides

% Extract 24h data
cAmat = 2.^table2array(TabA(ia,24:29)); % 23:28
cVmat = 2.^table2array(TabV(iv,24:29)); % 23:28


%%
% calculate mean 24h infected over mock
cAmatRatio = mean(cAmat(:,1:3),2) ./ mean(cAmat(:,4:6),2);
cVmatRatio = mean(cVmat(:,1:3),2) ./ mean(cVmat(:,4:6),2);


% 
% cAmat = (mean(cAmat(:,1:3),2) ./ mean(cAmat(:,4:6),2))  .* TabA.Intensity(ia);
% cVmat = (mean(cVmat(:,1:3),2) ./ mean(cVmat(:,4:6),2))  .* TabV.Intensity(iv);

cAmat = mean(cAmat(:,1:3),2) .* TabA.Intensity(ia);
cVmat = mean(cVmat(:,1:3),2) .* TabV.Intensity(iv);

% get hits in separate table and sort everything equivalently
TabAcutdown = TabA(ia,:);

[~,idx] = sort(TabAcutdown.TMT_24hA,'descend');
TabAcutdown = TabAcutdown(idx,:);

cAmat = cAmat(idx);
cVmat = cVmat(idx);

cAmatRatio = cAmatRatio(idx);
cVmatRatio = cVmatRatio(idx);

% Hits are:
% PNN (1)
% Viral (2-9)
% Nup107 (10)
% XRCC1 (15)
hits = [1:10,15];

%% Provide text annotation
TabAcutdown.GN(2:9) = TabAcutdown.LeadingRazorProtein(2:9);
TabAcutdown.GN(2:9) = extractAfter(TabAcutdown.GN(2:9),'|');
TabAcutdown.GN(2:9) = extractBefore(TabAcutdown.GN(2:9),'|');

%%
figure
subplot(1,2,1)
scatter(log2(cAmat), log2(cVmat),'filled','k','MarkerFaceAlpha',0.3);
hold on
scatter(log2(cAmat(hits)), log2(cVmat(hits)),'filled','r');

%title('Correlation of neo-N-termini intensity shared between the A549-Ace2 and Vero E6 datasets')
xlabel('Log_2 24h peptide intensity in A549-Ace2')
ylabel('Log_2 24h peptide intensity in Vero E6')
xlim([20,40])
ylim([20,40])

line([20,40],[20,40],'Color','k','LineStyle','--');
text(log2(cAmat(hits)), log2(cVmat(hits)),TabAcutdown.GN(hits),'FontSize',14);


% Calculate R2
R2 = corrcoef(log2(cAmat), log2(cVmat))
text(22,38,['Pearsons \rho = ',num2str(R2(1,2))],'FontSize',14);
set(gca,'FontSize',14);

writematrix([log2(cAmat), log2(cVmat)],'FigS2_A.csv')

subplot(1,2,2)

scatter(log2(cAmatRatio), log2(cVmatRatio),'filled','k','MarkerFaceAlpha',0.3);
hold on
scatter(log2(cAmatRatio(hits)), log2(cVmatRatio(hits)),'filled','r');
%title('Correlation of neo-N-termini intensity shared between the A549-Ace2 and Vero E6 datasets')
xlabel('Mean Log_2 24h Infected/24h Mock in A549-Ace2')
ylabel('Mean Log_2 24h Infected/24h Mock in Vero E6')
xlim([-2,7])
ylim([-2,7])

line([-2,7],[-2,7],'Color','k','LineStyle','--');
text(log2(cAmatRatio(hits)), log2(cVmatRatio(hits)),TabAcutdown.GN(hits),'FontSize',14);

% Calculate R2
R2 = corrcoef(log2(cAmatRatio), log2(cVmatRatio))
text(-1,6.1,['Pearsons \rho = ',num2str(R2(1,2))],'FontSize',14);

set(gca,'FontSize',14);


writematrix([log2(cAmatRatio), log2(cVmatRatio)],'FigS2_B.csv')

print([path , '/Figures/SFig_CorrelationScatter.pdf'],'-dpdf');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save all variables for analysis in separate downstream HIquant analysis
% script.


%% Generate files for downstream HIquant analysis of cellular proteins

% Need to: 
% Prelim: subset the sequences of interest
HIquant = struct();
HIquant.neoNterm.A549 = TabA(contains(TabA.Sequence,a549seqs),:);
HIquant.neoNterm.Vero = TabV(contains(TabV.Sequence,veroseqs),:);

% Note A549 includes an extra non-LxGG SRC cleavage site, possibly from
% cellular proteases as per previous published descriptions. This is removed manually
HIquant.neoNterm.A549(4,:) = [];


% 1. import unenriched peptide data
dat.Unfrpep.A549 = readtable([path , 'tA549_TotalFrac/peptides.txt']);
dat.Unfrpep.Vero = readtable([path , 'tVero_TotalFrac/peptides.txt']);

% 2. normalise as for the enriched

% Reorder the corrected RI intensity channels
dat.Unfrpep.A549(:,44:59) = dat.Unfrpep.A549(: , dat.tmt.A549(2,:) + 43);
dat.Unfrpep.Vero(:,44:59) = dat.Unfrpep.Vero(: , dat.tmt.Vero(2,:) + 43);
% 
% 
% Backup temp
backup = dat;

%% Data cleanup
dat = backup;

% Remove Reverse hits, Contaminants and low PEP identifications
for ii = 1:numel(samples)
    % Remove PEP <= 0.02
    dat.Unfrpep.(samples{ii}) = dat.Unfrpep.(samples{ii})(dat.Unfrpep.(samples{ii}).PEP <= 0.02 , :);
    % Remove Reverse hits
    dat.Unfrpep.(samples{ii}) = dat.Unfrpep.(samples{ii})(categorical(dat.Unfrpep.(samples{ii}).Reverse) ~= '+' , :);
    % Remove potential contaminants
    dat.Unfrpep.(samples{ii}) = dat.Unfrpep.(samples{ii})(categorical(dat.Unfrpep.(samples{ii}).PotentialContaminant) ~= '+' , :);
    % Filter on PIF

end

% Now Convert the TMT RI data to matrix form

for ii = 1:numel(samples)
   % Convert table columns to a matrix
    dat.Umat.(samples{ii}) = table2array(dat.Unfrpep.(samples{ii})(: , [44:59]));
   
   % Convert 0 to NaN;
   dat.Umat.(samples{ii})(dat.Umat.(samples{ii}) == 0) = NaN; 
end

% Remove unquantified hits, column normalise, impute and row normalise.
for ii = 1:numel(samples)
   % identify rows with >75% NaN
   allNaN = sum(isnan(dat.Umat.(samples{ii})) , 2) >13;
   dat.Umat.(samples{ii}) = dat.Umat.(samples{ii})(~allNaN , :);
   dat.Unfrpep.(samples{ii}) = dat.Unfrpep.(samples{ii})(~allNaN , :);
   
   clear allNaN
   
   % Median normalise
    dat.Umat.(samples{ii}) = dat.Umat.(samples{ii}) ./ nanmedian(dat.Umat.(samples{ii}));

    % Knn impute missing data
    dat.Umat.(samples{ii}) = knnimpute(dat.Umat.(samples{ii}));

    % Row normalise by mean (whole dataset) 
    dat.Umat.(samples{ii}) = dat.Umat.(samples{ii}) ./ mean(dat.Umat.(samples{ii}) , 2);
 
    % Generate matrix with 24h data only
    dat.Uh24.(samples{ii}) = dat.Umat.(samples{ii})(:,[end-5: end]);
   
end

backup2 = dat;

%%
% 3. for those select candidates: 
% a) Subset the data on cellular proteins:
for ii = 1:numel(samples)
    prots = HIquant.neoNterm.(samples{ii}).GN;
   for jj = 1:numel(HIquant.neoNterm.(samples{ii}).GN)
      HIquant.cell.(samples{ii}).(prots{jj,1}) = dat.Unfrpep.(samples{ii})(ismember(dat.Unfrpep.(samples{ii}).LeadingRazorProtein,HIquant.neoNterm.(samples{ii}).LeadingRazorProtein(jj)),:);
      HIquant.cellMat.(samples{ii}).(prots{jj,1}) = dat.Umat.(samples{ii})(ismember(dat.Unfrpep.(samples{ii}).LeadingRazorProtein,HIquant.neoNterm.(samples{ii}).LeadingRazorProtein(jj)),:);  
   end
end

% b) Loop through each protein and make sure there is no overlap with
% peptides and cleavage site: this makes analysis simple, but potentially
% loses some power to infer proteoforms.

% So needs to either start after end position, or end before start position
for ii = 1:numel(samples)
    prots = HIquant.neoNterm.(samples{ii}).GN;
    for jj = 1:numel(HIquant.neoNterm.(samples{ii}).GN)
    
        HIquant.cell.(samples{ii}).(prots{jj,1}) = HIquant.cell.(samples{ii}).(prots{jj,1})(HIquant.cell.(samples{ii}).(prots{jj}).StartPosition > HIquant.neoNterm.(samples{ii}).EndPosition(jj) + HIquant.cell.(samples{ii}).(prots{jj}).EndPosition < HIquant.neoNterm.(samples{ii}).StartPosition(1) == 1 , :);
        HIquant.cellMat.(samples{ii}).(prots{jj,1}) = HIquant.cellMat.(samples{ii}).(prots{jj,1})(HIquant.cell.(samples{ii}).(prots{jj}).StartPosition > HIquant.neoNterm.(samples{ii}).EndPosition(jj) + HIquant.cell.(samples{ii}).(prots{jj}).EndPosition < HIquant.neoNterm.(samples{ii}).StartPosition(1) == 1 , :);

    end
end

%% Tweak to add clustering filter: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c) group on a per-protein/cleavage site basis.

% HIquant requires input tables, without column headers, with:
% Col 1: colon-separated protein IDs
% Col 2: modified/peptide sequence
% Col 3+ level in conditions 1-n

for ii = 1:numel(samples)
       prots = HIquant.neoNterm.(samples{ii}).GN;
       HIquant.exportFinal.(samples{ii}) = table();
    % For each gene name in the neo-N-termini dataset (one cleavage site per GN)   
    for jj = 1:numel(HIquant.neoNterm.(samples{ii}).GN)
        HIquant.export.(samples{ii}).(prots{jj}).Seq = HIquant.cell.(samples{ii}).(prots{jj}).Sequence;
        HIquant.export.(samples{ii}).(prots{jj}).Mat = HIquant.cellMat.(samples{ii}).(prots{jj}); % Note: non-Neo-N-terminal peptides not log2 yet so doing here. neo-N-termini were immediately prior to exporting tables above.
        for kk = 1:numel(HIquant.export.(samples{ii}).(prots{jj}).Seq)
           HIquant.export.(samples{ii}).(prots{jj}).ProteinID{kk,1} = [prots{jj},';',prots{jj},'_',num2str(HIquant.neoNterm.(samples{ii}).StartPosition(jj)),'cl'];%prots{jj}; 
        end
        
            % Next section only performed if peptides from the protein
            % exist in the unenriched dataset.
            if isempty(HIquant.export.(samples{ii}).(prots{jj}).Mat) == 0 
                % Need to add in the cleaved peptide
                HIquant.export.(samples{ii}).(prots{jj}).Seq{kk + 1,1} = HIquant.neoNterm.(samples{ii}).Sequence{jj};
                HIquant.export.(samples{ii}).(prots{jj}).Mat(kk + 1,:) = 2.^table2array(HIquant.neoNterm.(samples{ii})(jj,[13:28]));
                HIquant.export.(samples{ii}).(prots{jj}).ProteinID{kk + 1,1} = [prots{jj},'_',num2str(HIquant.neoNterm.(samples{ii}).StartPosition(jj)),'cl'];

                % Convert output to table and rearrange into HIquant input
                % format.
                HIquant.exportTab.(samples{ii}).(prots{jj}) = struct2table(HIquant.export.(samples{ii}).(prots{jj}));
                HIquant.exportTab.(samples{ii}).(prots{jj}) = HIquant.exportTab.(samples{ii}).(prots{jj})(:,[3,1,2]);
                HIquant.exportTab.(samples{ii}).(prots{jj});%
                % Concatonate within each sample
                %HIquant.exportFinal.(samples{ii}) = [HIquant.exportFinal.(samples{ii}); HIquant.exportTab.(samples{ii}).(prots{jj})];
            end
    
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%% End clustering filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Write all these outputs to put through HIquant
    if ~exist('/Users/ed/Documents/HiquantTest/', 'dir')
       mkdir('/Users/ed/Documents/HiquantTest')
    end

writetable(HIquant.exportTab.A549.PAICS,'/Users/ed/Documents/HiquantTest/Input.txt','Delimiter','tab','WriteVariableNames',0);

% Absolutely CRUCIAL - test HIquant from terminal first to be sure of
% correct installation and all library requirements are in place. THEN you
% MUST start matlab from the command line.
tic
runs = [];

   runs = length(fieldnames(HIquant.exportTab.A549)) + length(fieldnames(HIquant.exportTab.Vero)); 

kk = 1;
for ii = 1:numel(samples)
    fnames = fieldnames(HIquant.exportTab.(samples{ii}));
    for jj = 1:length(fieldnames(HIquant.exportTab.(samples{ii})));
        
        writetable(HIquant.exportTab.(samples{ii}).(fnames{jj}),['/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Input.txt'],'Delimiter','tab','WriteVariableNames',0);
        
        commd = ['python3 ~/Documents/GitHub/HIquant/HIquant_run.py ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Input.txt'' ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Output_1.csv'' ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Output_2.csv'''];
        %commdV = ['python3 ~/Documents/GitHub/HIquant/HIquant_run.py ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_InputEE.txt'' ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Output_1EE.csv'' ''/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Output_2EE.csv'''];
      
        system(commd);

        disp([num2str(kk),' of ',num2str(runs)])
        kk = kk + 1;
    end
end
disp('Matlab HIquant loop completed')
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate CV's for cellular peptides

for ii = 1:numel(samples)
    
    fnames = fieldnames(HIquant.exportTab.(samples{ii}));
    
        for jj = 1:length(fnames)
            
            mat = table2array(HIquant.exportTab.(samples{ii}).(fnames{jj})(:,3));
            
            if size(mat,1) > 2
                HIquant.QC.(samples{ii}).Prots{jj,1} = fnames{jj};
                HIquant.QC.(samples{ii}).(fnames{jj}).CVs =   std(mat(1:end-1,:)) ./ mean(mat(1:end-1,:));
                %HIquant.QC.(samples{ii}).MedSD(jj,1) = median(std(mat(1:end-1,:)));
                HIquant.QC.(samples{ii}).MedCV(jj,1) = median(std(mat(1:end-1,:)) ./ mean(mat(1:end-1,:)));
                
                % Identify number of clusters in the unenriched peptide
                % data
                eva = evalclusters(HIquant.cellMat.(samples{ii}).(fnames{jj,1}),'linkage','gap','KList',[1:6])
                HIquant.QC.(samples{ii}).Clusters(jj,1) = eva.OptimalK;
                clear mat
            else
                
                HIquant.QC.(samples{ii}).Prots{jj,1} = fnames{jj};
                HIquant.QC.(samples{ii}).(fnames{jj}).CVs =   NaN;
                HIquant.QC.(samples{ii}).Clusters(jj,1) = NaN;
                %HIquant.QC.(samples{ii}).MedSD(jj,1) = NaN;
                HIquant.QC.(samples{ii}).MedCV(jj,1) = NaN;
                clear mat eva
                
            end
     
        end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Note SRC is identified in the unenriched fraction by a single peptide,
% but is added back into the analysis as it forms a central part of the
% story. 
HIquant.QC.A549.Clusters(3) = 1;

%% Import and display analysis

% Need to:

HIqImport = struct();

% Loop created files into HIqImport
for ii = 1:numel(samples)
    fnames = fieldnames(HIquant.exportTab.(samples{ii}));
    fnames = fnames(HIquant.QC.(samples{ii}).Clusters == 1); % Test;
        for jj = 1:length(fnames) %(fieldnames(HIquant.exportTab.(samples{ii})))
           % Import Inferred proteoform stochiometry
            HIqImport.(samples{ii}).stoich.(fnames{jj}) = readtable(['/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'_Output_1.csv']);
           % Remove first row
           HIqImport.(samples{ii}).stoich.(fnames{jj})(1,:) = [];
           % Import errors/confidence limits
           HIqImport.(samples{ii}).errors.(fnames{jj}) = readtable(['/Users/ed/Documents/HiquantTest/',samples{ii},fnames{jj},'Output_2.csv']);
           % Remove first row
           HIqImport.(samples{ii}).errors.(fnames{jj})(1,:) = [];
        end
end

% Generate tables and visualisations

% First up: immediate interest is fraction at each timepoint, as well as
% the total


%
% Separate the timepoints/replicates along the X axis
x = [0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5];
x2 = [0,1,2,3,4,5];

%HIqImport.A549.allerrors = table();
%HIqImport.Vero.allerrors = table();

%HIqImport.A549.allerrors(:,1) = HIqImport.A549.errors.NUP107.Var1;
%HIqImport.Vero.allerrors(:,1) = HIqImport.Vero.errors.ATAD2.Var1;

%%
% Now plot
% TODO: Add in total abundance:
fig = figure
t = tiledlayout(2,3)
set(fig,'defaultAxesColorOrder',[0 0 0;0 0 1]);
for ii = 1:numel(samples)
   fnames = fieldnames(HIqImport.(samples{ii}).stoich);
   %fnames = fnames(HIquant.QC.(samples{ii}).Clusters == 1); % test
   

   
   
   
   for jj = 1:numel(fnames)
       nexttile
       % Total abundance
       %yyaxis right
       scatter(x,...
           sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))) ./ (mean(sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))))),...
           'filled','b','MarkerFaceAlpha',0.3)
           
       writematrix(sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))) ./ (mean(sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))))),['FigS14_total',fnames{jj},'.csv'])
       
       hold on
       % Uncleaved protein abundance
       yyaxis left
        scatter(x,...
            table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))),...
            'filled','k','MarkerFaceAlpha',0.3)

        writematrix(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))),['FigS14_uncleaved',fnames{jj},'.csv'])
        
        %cleaved protein abundance
        scatter(x,...
            table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(2,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))),...
            'filled','r','MarkerFaceAlpha',0.3)
        
        writematrix( table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(2,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))),['FigS14_cleaved',fnames{jj},'.csv'])
        
        text(-0.5,1.4,['R_2: ',num2str(HIqImport.(samples{ii}).errors.(fnames{jj}){1,2},'%.2f'),'; Overall Error: ',...
            num2str(HIqImport.(samples{ii}).errors.(fnames{jj}){2,2}, '%.2f'),'; Ratio Error: ',...
            num2str(HIqImport.(samples{ii}).errors.(fnames{jj}){3,2}, '%.2f')],'FontSize',11)
  
        
                % Notes: for the purposes of this plotting, replicate the 0h Mock sample to
% give 18 samples and simplify analysis. THis will not alter interpretation
% as this sample will not show error bars as it will have SD of 0, so this
% data should be purely used for graphing purposes.
        
        % uncleaved
        aa = table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17)));
        
        aa = [aa(1), aa(1), aa];
        aa = reshape(aa,3,[]);
        aa2 = mean(aa);
        aaErr = std(aa);
        
        errorbar(x2,aa2,aaErr,'k','LineStyle','none');
        errorbar(x2(2:5),aa2(2:5),aaErr(2:5),'k');

        % cleaved
        bb = table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(2,2:17)) ./ ...
            sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17)));
        
        bb = [bb(1), bb(1), bb];
        bb = reshape(bb,3,[]);
        bb2 = mean(bb);
        bbErr = std(bb);
        
        errorbar(x2,bb2,bbErr,'r','LineStyle','none');
        errorbar(x2(2:5),bb2(2:5),bbErr(2:5),'r');
%         ylabel('Fraction of whole')
%         ylim([-0.2 1.5])
        
        % total
        %yyaxis right
        cc = sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))) ./ mean(sum(table2array(HIqImport.(samples{ii}).stoich.(fnames{jj})(1:2,2:17))));
        
        cc = [cc(1), cc(1), cc];
        cc = reshape(cc,3,[]);
        cc2 = mean(cc);
        ccErr = std(cc);
        errorbar(x2,cc2,ccErr,'b','LineStyle','none'); % No line wanted to connect the mocks
        errorbar(x2(2:5),cc2(2:5),ccErr(2:5),'b');
        
        % Adding annotation to separate the mocks from the samples
        line([0.5,0.5],[-0.2,1.5],'Color','k','LineStyle',':');
        line([4.5,4.5],[-0.2,1.5],'Color','k','LineStyle',':');
        % End annotation
        
        
        title([samples{ii},' ',fnames{jj}])
        xlabel('Time post-infection (h)')
        ylabel('Fraction of whole, per timepoint')
        xlim([-1 6])
        xticks([0:1:5])
        xticklabels({'0h Mock','0h','6h','12h','24h','24h Mock'})
        ylim([-0.2 1.5])
        legend({'Total','Uncleaved','Cleaved'},'Location','eastoutside')
        yyaxis right
                ylabel('Total abundance/mean')
        ylim([-0.2 1.5])
        
        hold off

        
   end
    
end



print([path , '/Figures/Fig_cleavagestoichiometry.pdf'],'-dpdf');
