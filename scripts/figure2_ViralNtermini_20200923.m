%% Identifying viral N-termini and Neo-N-termini

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Goal is SData tables, as well as 5 panels for Fig. 2.
% - Panels 1-3 will show S, ORF3a and N
% - Panels 4-5 will show the S cleavage site
% - Panel 4 will show S cleavage site abundance in Vero/A549.
% - Panel 5 will show this normalised to total S.


% Load data
dat = struct();


% Load the Var search evidence files
dat.var.A549 = readtable([path , 'tA549_EnrichVar/evidence.txt']);
dat.var.Vero = readtable([path , 'tVero_EnrichVar/evidence.txt']);

% Load the peptide files
dat.pep.A549 = readtable([path , 'tA549_EnrichVar/peptides.txt']);
dat.pep.Vero = readtable([path , 'tVero_EnrichVar/peptides.txt']);

samples = {'A549','Vero'};

%% Import gene name/signal peptide data from uniprot

viralseq = fastaread('/Users/ed/Dropbox/Liverpool/TempAnalysis/SARS2_Custom_20200518.fasta');

% Backup temp
backup = dat;
%% Data cleanup
dat = backup;

% Remove Reverse hits, Contaminants and low PEP identifications
for ii = 1:numel(samples)
    % Remove PEP <= 0.02
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).PEP <= 0.02 , :);
    % Remove Reverse hits
    dat.var.(samples{ii}) = dat.var.(samples{ii})(categorical(dat.var.(samples{ii}).Reverse) ~= '+' , :);
    % Remove potential contaminants
    dat.var.(samples{ii}) = dat.var.(samples{ii})(categorical(dat.var.(samples{ii}).PotentialContaminant) ~= '+' , :);
    % Filter on PIF >= 0.7
    dat.var.(samples{ii}) = dat.var.(samples{ii})(dat.var.(samples{ii}).PIF >= 0.7 , :);
end



for ii = 1:numel(samples)
   lia = contains(dat.var.(samples{ii}).Proteins , 'SARS2');
   dat.var.(samples{ii}) = dat.var.(samples{ii})(lia , :); 
   clear lia
end
    

%% Want to now sort on: has TMT/as N-terminal, has pyroglu? preceeded by?

% All modified N-termini begin '_(' (Acetyl, TMT, Pyroglu)
for ii = 1:numel(samples)
   lia = contains(dat.var.(samples{ii}).ModifiedSequence , '_('); 
   dat.var.(samples{ii}) = dat.var.(samples{ii})(lia , :);
   clear lia
end

% Want to work out:
% A) Which are top-scoring - sort on descending score
for ii = 1:numel(samples)

    [~ , idx] = sort(dat.var.(samples{ii}).Score , 'descend');
    dat.var.(samples{ii}) = dat.var.(samples{ii})(idx , :);
end

% B) Which are unique (sort on unique sequence as all are N-terminally
% modified) - as sorted by score this will take the first (i.e. highest
% scoring)
for ii = 1:numel(samples)
    [c , ia, ic] = unique(dat.var.(samples{ii}).Sequence);
    dat.var.(samples{ii}) = dat.var.(samples{ii})(ia , :);
end

%%
% C) Where is the cleavage site (i.e. first AA in peptide) - grab from the peptides file based on
% sequence
for ii = 1:numel(samples)
   [lia , locb] = ismember(dat.var.(samples{ii}).Sequence , dat.pep.(samples{ii}).Sequence);
   dat.var.(samples{ii}).StartPosition   = dat.pep.(samples{ii}).StartPosition(locb);
   % Useful for next step - P1 AA N-terminal of peptide.
   dat.var.(samples{ii}).AminoAcidBefore = dat.pep.(samples{ii}).AminoAcidBefore(locb);
end

% D) Pyrogly - if preceeded by R or K is lower confidence.
%%
for ii = 1:numel(samples)
   % Identify pyroglu peptides
   dat.pg.(samples{ii}) = contains(dat.var.(samples{ii}).ModifiedSequence , '_(G');
   dat.r.(samples{ii}) = contains(dat.var.(samples{ii}).AminoAcidBefore , 'R');
   dat.k.(samples{ii}) = contains(dat.var.(samples{ii}).AminoAcidBefore , 'K');
   % Logical sort based on pyrogly and R/K
   dat.pgrk.(samples{ii}) = [dat.pg.(samples{ii}) , dat.r.(samples{ii}) , dat.k.(samples{ii})]
   dat.logpgrk.(samples{ii}) = sum(dat.pgrk.(samples{ii}) , 2) >= 2;
end

% Generate two separate lists: 
% 1. High confidence (- pyroGlu preceeded by R or K)
% 2. Lower confidence (just those excluded)
%%
for ii = 1:numel(samples)
  dat.pgvar.(samples{ii}) = dat.var.(samples{ii})(dat.logpgrk.(samples{ii}) == 1 , :);
  dat.var.(samples{ii})   = dat.var.(samples{ii})(dat.logpgrk.(samples{ii}) == 0 , :);
end
%% Want to work out how many unique viral proteins I identified Neo-N-termini for.
dat.viralProt = unique([dat.var.A549.LeadingRazorProtein ; dat.var.Vero.LeadingRazorProtein]);

% Note: comes out as 8, as NiORF1 and 9B are identical.
%% Calculate sequence length for all fasta entries
for ii = 1:numel(viralseq)
    viralseq(ii).length = length(viralseq(ii).Sequence);
end

%% Plotting
%%%%%%%%%%%%%%%%%%% Nucleocapsid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define settings
 colors = {'k','r'}
 pos = [2 , 3 ,4 , 5] + 4;
 lnwdth = 1;


 ia = struct();
 % Now add N-termini
 % Obtain logical vectors for N peptides
 for ii = 1:numel(samples)
    ia.nuc.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|N|Nucleoprotein'); 
 end
 
 %
 figs = struct();
 
 % This sectionn notes the p1' location (N-terminus)
 for ii = 1:numel(samples)
    figs.nuc.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.nuc.(samples{ii}))
 end

% section grabs peptide intensity
int = struct;
for ii = 1:numel(samples)
    int.nuc.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.nuc.(samples{ii}));
end


figure
subplot(4,1,1)
line([0,419],[...
    min(log2([int.nuc.(samples{1}) ; int.nuc.(samples{2})]))-1 ,...
    min(log2([int.nuc.(samples{1}) ; int.nuc.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.nuc.(samples{ii}))
        line([figs.nuc.(samples{ii})(jj) , figs.nuc.(samples{ii})(jj)],...
            [min(log2([int.nuc.(samples{1}) ; int.nuc.(samples{2})]))-1 ,
            log2(int.nuc.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
aa1 =       scatter(figs.nuc.(samples{ii}),log2(int.nuc.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
aa2 =       scatter(figs.nuc.(samples{ii}),log2(int.nuc.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 419]);
ylim([min(log2([int.nuc.(samples{1}) ; int.nuc.(samples{2})]))-1,... 
    max(log2([int.nuc.(samples{1}) ; int.nuc.(samples{2})]))+3]);

title('Nucleocapsid (N) neo-N-termini')
xlabel('Amino Acid Position')
hold off
legend([aa1,aa2],{'A549-Ace2','Vero E6'})
legend('boxoff')
set(gca,'Fontsize',14)

 %%%%%%%%%%%%%%%%%%%%%%%% ORF 3a  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  % Obtain logical vectors for 3a peptides
 for ii = 1:numel(samples)
    ia.o3a.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|ORF3A|ORF'); 
 end
 
  for ii = 1:numel(samples)
    figs.o3a.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.o3a.(samples{ii}))
  end
 

for ii = 1:numel(samples)
    int.o3a.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.o3a.(samples{ii}));
end

 subplot(4,1,2)
line([0,276],[...
    min(log2([int.o3a.(samples{1}) ; int.o3a.(samples{2})]))-1 ,...
    min(log2([int.o3a.(samples{1}) ; int.o3a.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.o3a.(samples{ii}))
        line([figs.o3a.(samples{ii})(jj) , figs.o3a.(samples{ii})(jj)],...
            [min(log2([int.o3a.(samples{1}) ; int.o3a.(samples{2})]))-1 ,
            log2(int.o3a.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.o3a.(samples{ii}),log2(int.o3a.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.o3a.(samples{ii}),log2(int.o3a.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 276]);
ylim([min(log2([int.o3a.(samples{1}) ; int.o3a.(samples{2})]))-1,... 
    max(log2([int.o3a.(samples{1}) ; int.o3a.(samples{2})]))+1]);

title('ORF3a neo-N-termini')
xlabel('Amino Acid Position')
hold off
legend([ab1,ab2],{'A549-Ace2','Vero E6'})
legend('boxoff')
set(gca,'Fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Spike %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Obtain logical vectors for S peptides
 for ii = 1:numel(samples)
    ia.spk.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|S|Spike'); 
 end
 
  for ii = 1:numel(samples)
    figs.spk.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.spk.(samples{ii}))
  end
 
for ii = 1:numel(samples)
    int.spk.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.spk.(samples{ii}));
end

subplot(4,1,3)
line([0,1273],[...
    min(log2([int.spk.(samples{1}) ; int.spk.(samples{2})]))-1 ,...
    min(log2([int.spk.(samples{1}) ; int.spk.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.spk.(samples{ii}))
        line([figs.spk.(samples{ii})(jj) , figs.spk.(samples{ii})(jj)],...
            [min(log2([int.spk.(samples{1}) ; int.spk.(samples{2})]))-1 ,
            log2(int.spk.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ac1 =       scatter(figs.spk.(samples{ii}),log2(int.spk.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ac2 =       scatter(figs.spk.(samples{ii}),log2(int.spk.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 1273]);
ylim([min(log2([int.spk.(samples{1}) ; int.spk.(samples{2})]))-1,... 
    max(log2([int.spk.(samples{1}) ; int.spk.(samples{2})]))+1]);

title('Spike (S) neo-N-termini')
xlabel('Amino Acid Position')
hold off
legend([ac1,ac2],{'A549-Ace2','Vero E6'})
legend('boxoff')
set(gca,'Fontsize',14)

% Now do for remaining identified viral proteins to prepare supplemental tables 

%  %%%%%%%%%%%%%%%%%%%%%%%% Membrane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Obtain logical vectors for M peptides
 for ii = 1:numel(samples)
    ia.mem.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|M|Membrane'); 
 end
 
  for ii = 1:numel(samples)
    figs.mem.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.mem.(samples{ii}))
  end

%int = struct;
for ii = 1:numel(samples)
    int.mem.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.mem.(samples{ii}));
end
 

  %%%%%%%%%%%%%%%%%%%%%%%% ORF 7a iORF1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Obtain logical vectors for 3a peptides
 for ii = 1:numel(samples)
    ia.o7ai.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|ORF7aiORF1|ORF'); 
 end
 
  for ii = 1:numel(samples)
    figs.o7ai.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.o7ai.(samples{ii}))
  end
 
  
% int = struct;
for ii = 1:numel(samples)
    int.o7ai.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.o7ai.(samples{ii}));
end
 
 %%%%%%%%%%%%%%%%%%%%%%%% ORF 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  % Obtain logical vectors for ORF8 peptides
 for ii = 1:numel(samples)
    ia.o8.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|ORF8|ORF'); 
 end
 
  for ii = 1:numel(samples)
    figs.o8.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.o8.(samples{ii}))
  end
 
for ii = 1:numel(samples)
    int.o8.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.o8.(samples{ii}));
end

%%%%%%%%%%%%%%%%%%%%%%%% ORF 9B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Obtain logical vectors for ORF9B peptides
 for ii = 1:numel(samples)
    ia.o9b.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , {'SARS2|ORF9B|ORF', 'SARS2|NiORF1|Nucleoprotein'}); 
 end
 
  for ii = 1:numel(samples)
    figs.o9b.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.o9b.(samples{ii}))
  end
 
for ii = 1:numel(samples)
    int.o9b.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.o9b.(samples{ii}));
end

 %%%%%%%%%%%%%%%%%%%%%%%% pp1ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain logical vectors for pp1ab peptides
 for ii = 1:numel(samples)
    ia.pp1ab.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|pp1ab|Replicase');
 end
    
  for ii = 1:numel(samples)
    figs.pp1ab.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.pp1ab.(samples{ii}))
  end
  
for ii = 1:numel(samples)
    int.pp1ab.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.pp1ab.(samples{ii}));
end

 % Now import Quant data for panel looking at Spike abundace over time.
dat2 = struct
dat2.pep.A549 = readtable([path , 'tA549_Enrich/peptides.txt']);
dat2.pep.Vero = readtable([path , 'tVero_Enrich/peptides.txt']);

% TMT labelling strategies
% A549-ACE2
dat2.tmt.A549 = readmatrix([path , 'SARS2a549tmtlabelling_20200507.csv']);
% VeroE6
dat2.tmt.Vero = readmatrix([path , 'Verotmtlabelling_20200511.csv']);

% Reorder the corrected RI intensity channels to account for the labelling
% strategy
dat2.pep.A549(:,44:59) = dat2.pep.A549(: , dat2.tmt.A549(2,:) + 43);
dat2.pep.Vero(:,44:59) = dat2.pep.Vero(: , dat2.tmt.Vero(2,:) + 43);

 % Clean up data
 % Remove Reverse hits, Contaminants and low PEP identifications
for ii = 1:numel(samples)
    % Remove PEP <= 0.02
    dat2.pep.(samples{ii}) = dat2.pep.(samples{ii})(dat2.pep.(samples{ii}).PEP <= 0.02 , :);
    % Remove Reverse hits
    dat2.pep.(samples{ii}) = dat2.pep.(samples{ii})(categorical(dat2.pep.(samples{ii}).Reverse) ~= '+' , :);
    % Remove potential contaminants
    dat2.pep.(samples{ii}) = dat2.pep.(samples{ii})(categorical(dat2.pep.(samples{ii}).PotentialContaminant) ~= '+' , :);
end

% Now Convert the TMT RI data to matrix form

for ii = 1:numel(samples)
   % Convert table columns to a matrix
    dat2.mat.(samples{ii}) = table2array(dat2.pep.(samples{ii})(: , [44:59]));
   
   % Convert 0 to NaN;
   dat2.mat.(samples{ii})(dat2.mat.(samples{ii}) == 0) = NaN; 
end


% Remove unquantified hits, column normalise, impute and row normalise.
for ii = 1:numel(samples)
   % identify rows with >75% NaN
   allNaN = sum(isnan(dat2.mat.(samples{ii})) , 2) > 13;
   dat2.mat.(samples{ii}) = dat2.mat.(samples{ii})(~allNaN , :);
   dat2.pep.(samples{ii}) = dat2.pep.(samples{ii})(~allNaN , :);
   
   clear allNaN
   
   % Median normalise
    dat2.mat.(samples{ii}) = dat2.mat.(samples{ii}) ./ nanmedian(dat2.mat.(samples{ii}));

    % KNNimpute missing data
    dat2.mat.(samples{ii}) = knnimpute(dat2.mat.(samples{ii}));

    % Row normalise by mean (whole dataset) - not log2 data
    dat2.mat.(samples{ii}) = dat2.mat.(samples{ii}) ./ mean(dat2.mat.(samples{ii}) , 2);
   
end
 %
% Now need to identify my Spike peptide of interest in both samples
% STGSNVFQTR

for ii = 1:numel(samples)
    dat2.spikePeplog.(samples{ii}) = contains(dat2.pep.(samples{ii}).Sequence, 'STGSNVFQTR');
    dat2.spikeMat.(samples{ii})    = dat2.mat.(samples{ii})(dat2.spikePeplog.(samples{ii}),:);
    dat2.spikeInt.(samples{ii})    = dat2.pep.(samples{ii}).Intensity(dat2.spikePeplog.(samples{ii}));

end

dat2.spikeReformat.A549 = dat2.spikeMat.A549([2,3,4;5,6,7;8,9,10;11,12,13]);
dat2.spikeReformat.Vero = dat2.spikeMat.Vero([2,3,4;5,6,7;8,9,10;11,12,13]);

   dat2.spikeRfInt.A549 = log2(dat2.spikeInt.A549 .* (dat2.spikeReformat.A549 ./ sum(sum(dat2.spikeReformat.A549))));
   dat2.spikeRfInt.Vero = log2(dat2.spikeInt.Vero .* (dat2.spikeReformat.Vero ./ sum(sum(dat2.spikeReformat.Vero))));

timecourseTimesP = [0,6,12,24];
tcTicks = [0,3,6,9,12,15,18,21,24];
tcTicklabels = {'0','3','6','9','12','','18','','24'};
%
subplot(4,1,4)
errorbar(timecourseTimesP, mean(dat2.spikeRfInt.A549'),std(dat2.spikeRfInt.A549'),'-ok','LineWidth',1,'MarkerFaceColor','k');
hold on
errorbar(timecourseTimesP, mean(dat2.spikeRfInt.Vero'),std(dat2.spikeRfInt.Vero'),'-or','LineWidth',1,'MarkerFaceColor','r');

title('Spike S637 neo-N-terminus')
xlabel('Hours post-infection')
legend({'A549-Ace2','Vero E6'},'Location','southeast')
legend('boxoff')

xticks(tcTicks)
xticklabels(tcTicklabels)
xlim([-1,25])
ylim([17,27])
set(gca,'Fontsize',14)
hold off

[ax , h] = suplabel('Log2 Peptide Intensity','y')
set(h,'FontSize',16)    
    
% Save figure
print([path , '/Figures/Fig_ViralNterm.pdf'],'-dpdf');


 %% Write Supplementary Data Tables Covering confident viral neo-N-termini
 % Remove excess columns from tables
 % Remove 35-42,54,56-65
 dat.tbl.A549 = dat.var.A549(:,[1:34,43:53,55,66:67]);
 dat.tbl.Vero = dat.var.Vero(:,[1:34,43:53,55,66:67]);
 
 % Re-sort data by: a) start position, then b) protein
 [~ , ia] = sort(dat.tbl.A549.StartPosition,'ascend');
 dat.tbl.A549 = dat.tbl.A549(ia , :);
 [~ , ia] = sort(dat.tbl.A549.LeadingRazorProtein);
 dat.tbl.A549 = dat.tbl.A549(ia , :);
 
 clear ia
  [~ , ia] = sort(dat.tbl.Vero.StartPosition,'ascend');
 dat.tbl.Vero = dat.tbl.Vero(ia , :);
 [~ , ia] = sort(dat.tbl.Vero.LeadingRazorProtein);
 dat.tbl.Vero = dat.tbl.Vero(ia , :);
 clear ia
 
 % Write tables
 writetable(dat.tbl.A549 , [path , '/Tables/viralNeoN-termini_A549.csv']);
 
 writetable(dat.tbl.Vero , [path , '/Tables/viralNeoN-termini_VeroE6.csv']);