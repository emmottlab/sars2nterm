%% Identifying viral N-termini and Neo-N-termini
% Supplementary Figure (Further viral neo-N-termini)

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';
 
% - Panels 1-5 will show M, 7a, 8, 9B and pp1ab
% Troubleshooting - the tiledlayout feature used in this script was
% introduced in Matlab R2019b, it will not work with earlier versions.
% However, to make this script work in earlier versions, the tiledlayout
% and nexttile commands will need replacing with subplot.

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

viralseq = fastaread([path ,'SARS2_Custom_20200518.fasta']);

% Backup temp
backup = dat;
%% Data cleanup
dat = backup;

% Remove Reverse hits, Contaminants and low PEP and low PIF identifications
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

% Select for SARS2 proteins ONLY.
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
% A) Which are top-scoring - first sort on descending score
for ii = 1:numel(samples)

    [~ , idx] = sort(dat.var.(samples{ii}).Score , 'descend');
    dat.var.(samples{ii}) = dat.var.(samples{ii})(idx , :);
end

% B) Which are unique (sort on unique sequence as all are N-terminally
% modified) - as sorted by score this will take the first (i.e. highest
% scoring) for each peptide
for ii = 1:numel(samples)
    [c , ia, ic] = unique(dat.var.(samples{ii}).Sequence);
    dat.var.(samples{ii}) = dat.var.(samples{ii})(ia , :);
end


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
   % Identify pyroglu peptides and whether these follow R/K which could
   % indicate tryptic cleavage
   dat.pg.(samples{ii}) = contains(dat.var.(samples{ii}).ModifiedSequence , '_(G');
   dat.r.(samples{ii}) = contains(dat.var.(samples{ii}).AminoAcidBefore , 'R');
   dat.k.(samples{ii}) = contains(dat.var.(samples{ii}).AminoAcidBefore , 'K');
   % Logical sort based on pyrogly and R/K
   dat.pgrk.(samples{ii}) = [dat.pg.(samples{ii}) , dat.r.(samples{ii}) , dat.k.(samples{ii})];
   dat.logpgrk.(samples{ii}) = sum(dat.pgrk.(samples{ii}) , 2) >= 2;
end

% Generate two separate lists: 
% 1. High confidence (Potentially pyroGlu, but not preceeded by R/K as these could be experimental artifacts generated duriing trypsin digestion)
% 2. Lower confidence (just those excluded - pyroGlu, and preceeded by R/K)

% Separate out the high and low confidence lists.
for ii = 1:numel(samples)
  dat.pgvar.(samples{ii}) = dat.var.(samples{ii})(dat.logpgrk.(samples{ii}) == 1 , :);
  dat.var.(samples{ii})   = dat.var.(samples{ii})(dat.logpgrk.(samples{ii}) == 0 , :);
end
% Want to work out how many unique viral proteins I identified Neo-N-termini for.
dat.viralProt = unique([dat.var.A549.LeadingRazorProtein ; dat.var.Vero.LeadingRazorProtein]);

% Note: comes out as 8, as NiORF1 and 9B are identical.

% Calculate sequence length for all fasta entries
for ii = 1:numel(viralseq)
    viralseq(ii).length = length(viralseq(ii).Sequence);
end

%% Plotting
% Define settings
 colors = {'k','r'}
 pos = [2 , 3 ,4 , 5] + 4;
 lnwdth = 1;


 ia = struct();

 figs = struct();
 
int = struct;

%  %%%%%%%%%%%%%%%%%%%%%%%% Membrane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Obtain logical vectors for M peptides
 for ii = 1:numel(samples)
    ia.mem.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|M|Membrane'); 
 end
 
 % Identify the start position of the peptides within the protein
  for ii = 1:numel(samples)
    figs.mem.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.mem.(samples{ii}))
  end

%int = struct;
% Grab the peptide intensity for the peptides of interest
for ii = 1:numel(samples)
    int.mem.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.mem.(samples{ii}));
end

% Plot data
figure
t = tiledlayout(5,1); 
nexttile
line([0,222],[...
    min(log2([int.mem.(samples{1}) ; int.mem.(samples{2})]))-1 ,...
    min(log2([int.mem.(samples{1}) ; int.mem.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.mem.(samples{ii}))
        line([figs.mem.(samples{ii})(jj) , figs.mem.(samples{ii})(jj)],...
            [min(log2([int.mem.(samples{1}) ; int.mem.(samples{2})]))-1 ,
            log2(int.mem.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.mem.(samples{ii}),log2(int.mem.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.mem.(samples{ii}),log2(int.mem.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 222]);
ylim([min(log2([int.mem.(samples{1}) ; int.mem.(samples{2})]))-1,... 
    max(log2([int.mem.(samples{1}) ; int.mem.(samples{2})]))+1]);

title('Membrane (M) neo-N-termini')

hold off
legend([ab1,ab2],{'A549-Ace2','Vero E6'},'Location','west')
legend('boxoff')
set(gca,'Fontsize',14)

writematrix([figs.mem.A549,log2(int.mem.A549)],'FigS3_memA.csv')
writematrix([figs.mem.Vero,log2(int.mem.Vero)],'FigS3_memV.csv')


 %%%%%%%%%%%%%%%%%%%%%%%% ORF 7a iORF1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Note: the plotting and analysis steps for subsequent proteins are as for M.
 % Note is 7a, but begins at Amino acid 3. (119 vs 121 AA)


% Obtain logical vectors for 3a peptides
 for ii = 1:numel(samples)
    ia.o7a.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein , 'SARS2|ORF7aiORF1|ORF'); 
 end
 
  for ii = 1:numel(samples)
    figs.o7a.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.o7a.(samples{ii}))
  end

%int = struct;
for ii = 1:numel(samples)
    int.o7a.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.o7a.(samples{ii}));
end

nexttile

line([0,119],[...
    min(log2([int.o7a.(samples{1}) ; int.o7a.(samples{2})]))-1 ,...
    min(log2([int.o7a.(samples{1}) ; int.o7a.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.o7a.(samples{ii}))
        line([figs.o7a.(samples{ii})(jj) , figs.o7a.(samples{ii})(jj)],...
            [min(log2([int.o7a.(samples{1}) ; int.o7a.(samples{2})]))-1 ,
            log2(int.o7a.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.o7a.(samples{ii}),log2(int.o7a.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.o7a.(samples{ii}),log2(int.o7a.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 119]);
ylim([min(log2([int.o7a.(samples{1}) ; int.o7a.(samples{2})]))-1,... 
    max(log2([int.o7a.(samples{1}) ; int.o7a.(samples{2})]))+1]);

title({'ORF7a iORF1 neo-N-termini','(in-frame, begins ORF7a I3)'})
hold off
set(gca,'Fontsize',14)


writematrix([figs.o7a.A549,log2(int.o7a.A549)],'FigS3_o7aA.csv')
writematrix([figs.o7a.Vero,log2(int.o7a.Vero)],'FigS3_o7aV.csv')

 
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

nexttile
%subplot(5,1,3)
line([0,121],[...
    min(log2([int.o8.(samples{1}) ; int.o8.(samples{2})]))-1 ,...
    min(log2([int.o8.(samples{1}) ; int.o8.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.o8.(samples{ii}))
        line([figs.o8.(samples{ii})(jj) , figs.o8.(samples{ii})(jj)],...
            [min(log2([int.o8.(samples{1}) ; int.o8.(samples{2})]))-1 ,
            log2(int.o8.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.o8.(samples{ii}),log2(int.o8.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.o8.(samples{ii}),log2(int.o8.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 121]);
ylim([min(log2([int.o8.(samples{1}) ; int.o8.(samples{2})]))-1,... 
    max(log2([int.o8.(samples{1}) ; int.o8.(samples{2})]))+1]);

title({'ORF8 neo-N-termini'})
hold off

set(gca,'Fontsize',14)


writematrix([figs.o8.A549,log2(int.o8.A549)],'FigS3_o8A.csv')
writematrix([figs.o8.Vero,log2(int.o8.Vero)],'FigS3_o8V.csv')


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

nexttile
line([0,97],[...
    min(log2([int.o9b.(samples{1}) ; int.o9b.(samples{2})]))-1 ,...
    min(log2([int.o9b.(samples{1}) ; int.o9b.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.o9b.(samples{ii}))
        line([figs.o9b.(samples{ii})(jj) , figs.o9b.(samples{ii})(jj)],...
            [min(log2([int.o9b.(samples{1}) ; int.o9b.(samples{2})]))-1 ,
            log2(int.o9b.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.o9b.(samples{ii}),log2(int.o9b.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.o9b.(samples{ii}),log2(int.o9b.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 97]);
ylim([min(log2([int.o9b.(samples{1}) ; int.o9b.(samples{2})]))-1,... 
    max(log2([int.o9b.(samples{1}) ; int.o9b.(samples{2})]))+1]);

title({'ORF9B N-terminus'})
hold off

set(gca,'Fontsize',14)


writematrix([figs.o9b.A549,log2(int.o9b.A549)],'FigS3_o9bA.csv')
writematrix([figs.o9b.Vero,log2(int.o9b.Vero)],'FigS3_o9bV.csv')


 %%%%%%%%%%%%%%%%%%%%%%%% pp1ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain logical vectors for ORF9B peptides
 for ii = 1:numel(samples)
    ia.pp.(samples{ii}) = ismember(dat.var.(samples{ii}).LeadingRazorProtein ,'SARS2|pp1ab|Replicase'); 
 end
 
  for ii = 1:numel(samples)
    figs.pp.(samples{ii}) = dat.var.(samples{ii}).StartPosition(ia.pp.(samples{ii}))
  end

for ii = 1:numel(samples)
    int.pp.(samples{ii}) = dat.var.(samples{ii}).Intensity(ia.pp.(samples{ii}));
end

nexttile
line([0,7096],[...
    min(log2([int.pp.(samples{1}) ; int.pp.(samples{2})]))-1 ,...
    min(log2([int.pp.(samples{1}) ; int.pp.(samples{2})]))-1],'LineWidth',lnwdth,'Color','k')
 hold on
 
 for ii = 1:2
    for jj = 1:numel(int.pp.(samples{ii}))
        line([figs.pp.(samples{ii})(jj) , figs.pp.(samples{ii})(jj)],...
            [min(log2([int.pp.(samples{1}) ; int.pp.(samples{2})]))-1 ,
            log2(int.pp.(samples{ii})(jj))],'Color',[17,17,17]/255)
        
    end
 end
      ii = 1;
ab1 =       scatter(figs.pp.(samples{ii}),log2(int.pp.(samples{ii})),[],colors{ii},'filled'); 
      ii = 2;
ab2 =       scatter(figs.pp.(samples{ii}),log2(int.pp.(samples{ii})),[],colors{ii},'filled'); 

clear ii
xlim([0, 7096]);
ylim([min(log2([int.pp.(samples{1}) ; int.pp.(samples{2})]))-1,... 
    max(log2([int.pp.(samples{1}) ; int.pp.(samples{2})]))+1]);

title({'pp1ab neo-N-termini'})
xlabel('Amino Acid Position')
hold off

ylabel(t,'Log_2 Peptide Intensity','FontSize',14)
set(gca,'Fontsize',14)

writematrix([figs.pp.A549,log2(int.pp.A549)],'FigS3_ppA.csv')
writematrix([figs.pp.Vero,log2(int.pp.Vero)],'FigS3_ppV.csv')

% Save results
print([path , '/Figures/SFig_SupViralNtermini.pdf'],'-dpdf');
