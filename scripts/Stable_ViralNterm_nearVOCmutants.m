% SARS-CoV-2 N-terminomics: identification of viral cleavage sites lying
% proximal to mutations in currently circulating variants of interest.

% Characteristic VoI/VoC mutations identified from covariants.org, data obtained on 15th March 2021.

% Stains examined:
% B.1.1.1.7 aka 20I/501Y.V1
% B.1.351   aka 20H/501Y.V2
% P.1       aka 20J/501Y.V3
% B.1.427/9 aka 20C/S.452R
% B.1.526   aka 20C/S.484K
% B.1.525   aka 20A/S.484K


clear
clc

%% Import data
% Import stain/mutation table generated data available from CoVariants.org
dat = struct()

%dat.variants =
%readtable('/Users/ed/Documents/GitHub/sars2nterm/data/covariants_20210315.csv');
% Older list, depreciated.
dat.variants = readtable('/Users/ed/Documents/GitHub/sars2nterm/data/covariants_20210524.csv');

% Import viral cleavage data from tables S3/4.
dat.a549 = readtable('/Users/ed/Documents/GitHub/sars2nterm/data/Tables/viralNeoN-termini_A549.csv');
dat.vero = readtable('/Users/ed/Documents/GitHub/sars2nterm/data/Tables/viralNeoN-termini_VeroE6.csv');
dat.total = readtable('/Users/ed/Documents/GitHub/sars2nterm/data/Tables/ViralPeptides.csv');
%% Reannotate protein identification in dat.a549/vero to match that in dat.variants.

% Note that complete reannotation is unnecessary as only a subset of viral
% proteins are identified

% These are: S, ORF1a, N, ORF8, E, ORF1b,ORF3a.
% Note ORF1b is not accepted coronavirus nomenclature as it fails to
% account for the frameshift and N-terminal shared region with ORF1a.

% ORF1a: 1-4405
% ORF1b: 4393-7096

% Note: for the purposes of this script, the 4393-4405 region is considered the beginning of pp1b.

dat.a549.varStartPos = dat.a549.StartPosition;

dat.a549.varProt(contains(dat.a549.LeadingRazorProtein , 'SARS2|N|Nucleoprotein')) = {'N'};
dat.a549.varProt(contains(dat.a549.LeadingRazorProtein , 'SARS2|S|Spike')) = {'S'};
dat.a549.varProt(contains(dat.a549.LeadingRazorProtein , 'SARS2|ORF3A|ORF')) = {'ORF3a'};

% If is pp1ab & is >= 4393 start position: is pp1b. if <4393 is pp1a.
dat.a549.varProt(contains(dat.a549.LeadingRazorProtein , 'SARS2|pp1ab|Replicase')) = {'ORF1a'};

% The SARS2 proteins longer than 4393 are pp1b.
dat.a549.varProt(dat.a549.StartPosition >= 4393) = {'ORF1b'};
dat.a549.StartPosition(dat.a549.StartPosition >= 4393) = dat.a549.StartPosition(dat.a549.StartPosition >= 4393) - 4392;

% Now do for Vero
dat.vero.varStartPos = dat.vero.StartPosition;

dat.vero.varProt(contains(dat.vero.LeadingRazorProtein , 'SARS2|N|Nucleoprotein')) = {'N'};
dat.vero.varProt(contains(dat.vero.LeadingRazorProtein , 'SARS2|S|Spike')) = {'S'};
dat.vero.varProt(contains(dat.vero.LeadingRazorProtein , 'SARS2|ORF3A|ORF')) = {'ORF3a'};
dat.vero.varProt(contains(dat.vero.LeadingRazorProtein , 'SARS2|ORF8|ORF')) = {'ORF8'};

% If is pp1ab & is >= 4393 start position: is pp1b. if <4393 is pp1a.
dat.vero.varProt(contains(dat.vero.LeadingRazorProtein , 'SARS2|pp1ab|Replicase')) = {'ORF1a'};

% The SARS2 proteins longer than 4393 are pp1b.
dat.vero.varProt(dat.vero.StartPosition >= 4393) = {'ORF1b'};
dat.vero.StartPosition(dat.vero.StartPosition >= 4393) = dat.vero.StartPosition(dat.vero.StartPosition >= 4393) - 4392;

% Now do for total
dat.total.varStartPos = dat.total.StartPosition;

dat.total.varProt(contains(dat.total.Proteins , 'SARS2|N|Nucleoprotein')) = {'N'};
dat.total.varProt(contains(dat.total.Proteins , 'SARS2|S|Spike')) = {'S'};
dat.total.varProt(contains(dat.total.Proteins , 'SARS2|ORF3A|ORF')) = {'ORF3a'};
dat.total.varProt(contains(dat.total.Proteins , 'SARS2|ORF8|ORF')) = {'ORF8'};

% If is pp1ab & is >= 4393 start position: is pp1b. if <4393 is pp1a.
dat.total.varProt(contains(dat.vero.Proteins , 'SARS2|pp1ab|Replicase')) = {'ORF1a'};

% The SARS2 proteins longer than 4393 are pp1b.
dat.total.varProt(dat.total.StartPosition >= 4393) = {'ORF1b'};
dat.total.StartPosition(dat.total.StartPosition >= 4393) = dat.total.StartPosition(dat.total.StartPosition >= 4393) - 4392;

%%

% Now complete the remaining proteins so that matlab plays nice
dat.a549.varProt(find(cellfun(@isempty,dat.a549.varProt))) = {'_'};
dat.vero.varProt(find(cellfun(@isempty,dat.vero.varProt))) = {'_'};
dat.total.varProt(find(cellfun(@isempty,dat.total.varProt))) = {'_'};

%% Now match cleavage sites with variants data and identify how far each cleavage site is from its closest variant
% Concept: 
% 1. select those matching each variant protein:
% 2. Identify the mutation closest to each cleavage site
% 3. Identify the distance
% 4. Highlight if <5 and/or <10AA away

% Loop through every line in dat.variants
for ii = 1:numel(dat.variants.Protein)
    
    % Ace2-A549
    [lia, locb] = ismember(dat.a549.varProt , dat.variants.Protein(ii));
    
    % If there is at least one match to the protein
    if sum(lia) >= 1
        posVar = dat.variants.Position_AA(ii);
        NtermPos = dat.a549.StartPosition(lia);
        
        [~, minA] = min(abs(NtermPos-posVar));
        dat.variants.closestNterm_A549(ii) = NtermPos(minA);
        
        dat.variants.distFromClosest_A549(ii) = dat.variants.Position_AA(ii) - dat.variants.closestNterm_A549(ii);
      clear posVar NtermPos minA
    end
    
    dat.variants.within_10_AA_A549 = abs(dat.variants.distFromClosest_A549) <= 10;
    dat.variants.within_5_AA_A549 = abs(dat.variants.distFromClosest_A549) <= 5;
    
    
    
    % Vero E6
    [lia, locb] = ismember(dat.vero.varProt , dat.variants.Protein(ii));
    
    % If there is at least one match to the protein
    if sum(lia) >= 1
        posVar = dat.variants.Position_AA(ii);
        NtermPos = dat.vero.StartPosition(lia);
        
        [~, minA] = min(abs(NtermPos-posVar));
        dat.variants.closestNterm_Vero(ii) = NtermPos(minA);
        
        dat.variants.distFromClosest_Vero(ii) = dat.variants.Position_AA(ii) - dat.variants.closestNterm_Vero(ii);
      clear posVar NtermPos minA
    end
    
        dat.variants.within_10_AA_Vero = abs(dat.variants.distFromClosest_Vero) <= 10;
    dat.variants.within_5_AA_Vero = abs(dat.variants.distFromClosest_Vero) <= 5;
    
    mat = table2array(dat.variants(:,[5,6,9,10]));
    mat(mat == 0) = NaN;
    
    dat.variants(:,[5,6,9,10]) = array2table(mat);
    

end

% This is complete and ready for export as S. Table
writetable(dat.variants,'/Users/ed/Documents/GitHub/sars2nterm/data/Tables/Stable_NtermNearVOCmutations.csv');

%% So really want a reverse of the previous 

% Sort on proteins included in variant mutations
logNVOC = ismember(dat.total.varProt,'_');
dat.totalS = dat.total(~logNVOC,:);

% Loop through all entries in total
for ii = 1:numel(dat.totalS.Sequence)
      [lia, locb] = ismember(dat.variants.Protein, dat.totalS.varProt(ii));
 
      if sum(lia) >= 1
            posVar = dat.totalS.StartPosition(ii);
            NtermPos = dat.variants.Position_AA(lia);
            
          [~, minA] = min(abs(NtermPos-posVar));
         dat.totalS.closestNterm(ii) = NtermPos(minA);
         
         dat.totalS.distFromClosest(ii) = abs(dat.totalS.varStartPos(ii) - dat.totalS.closestNterm(ii));
       clear posVar NtermPos minA           
      end
      
end


%%

% Remove duplicates
[c, ia, ic] = unique(dat.totalS.Sequence);
dat.totalS = dat.totalS(ia,:);

% Now sort and do stats
[~ , Idx] = sort(dat.totalS.distFromClosest , 'ascend');

dat.totalS = dat.totalS(Idx,:);

% Indicate enriched
logEnr = contains(dat.totalS.RawFile , 'Enriched');

% Now do KS test
num = 1:numel(logEnr);
%%
figure
scatter(num(~logEnr),log10(dat.totalS.distFromClosest(~logEnr)),'filled','k');
hold on
scatter(num(logEnr),log10(dat.totalS.distFromClosest(logEnr)),'filled','r');
xlabel('Peptides ranked by distance from VOC mutation');
ylabel('Log_1_0 amino acid distance from VOC mutation');
legend({'Viral Peptides','Viral neo-N-termini'},'Location','northwest')
title('neo-N-terminus proximity to VOC mutations')


[h,p,k] = kstest2(num(logEnr),num);
text(40,0.5,['KS test: p =',num2str(p)],'FontSize',14);

set(gca,'FontSize',14);

%%
% Subset most interesting
temp = dat.variants.within_5_AA_A549 + dat.variants.within_5_AA_Vero >= 1;
dat.subsvariants = dat.variants(temp,:);

temp = dat.variants.within_10_AA_A549 + dat.variants.within_10_AA_Vero >= 1;
dat.subsvariants2 = dat.variants(temp,:);

% Remove redundancy
for ii = 1:numel(dat.subsvariants.Protein)
    dat.subsvariants.ProtPos{ii} = [dat.subsvariants.Protein{ii},'_', num2str(dat.subsvariants.Position_AA(ii))];
end
size(unique(dat.subsvariants.ProtPos),1) % How many unique sites make cutoff

for ii = 1:numel(dat.subsvariants2.Protein)
    dat.subsvariants2.ProtPos{ii} = [dat.subsvariants2.Protein{ii},'_', num2str(dat.subsvariants2.Position_AA(ii))];
end
size(unique(dat.subsvariants2.ProtPos),1) % How many unique sites make cutoff

% Number of VOC/VOIs 
categories(categorical(dat.variants.Strain))

% Total sites/mutations
for ii = 1:numel(dat.variants.Protein)
   dat.variants.ProtMut{ii} = [dat.variants.Protein{ii},'_',dat.variants.Mutation{ii}]; 
end

size(unique(dat.variants.ProtMut),1)

%%  Write tables with the subset closest
writetable(dat.subsvariants, '/Users/ed/Documents/GitHub/sars2nterm/data/Tables/VOC_neoNterm_within5aa.csv');
writetable(dat.subsvariants2,'/Users/ed/Documents/GitHub/sars2nterm/data/Tables/VOC_neoNterm_within10aa.csv');
