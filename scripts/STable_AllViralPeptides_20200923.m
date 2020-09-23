% All viral peptides script for table S1

% Script draws together all viral peptide identifications from the
% unenriched and enriched samples.

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data
% Need to import both evidence (for modified peptide sequence) and peptides
% (for position in protein)

% Going from both enriched and unenriched samples

dat = struct();

% Load Evidence files
dat.evi.A549tot = readtable([path , 'tA549_TotalFrac/evidence.txt']); 1
dat.evi.Verotot = readtable([path , 'tVero_TotalFrac/evidence.txt']); 2
dat.evi.A549enr = readtable([path , 'tA549_EnrichVar/evidence.txt']); 3
dat.evi.Veroenr = readtable([path , 'tVero_EnrichVar/evidence.txt']); 4

% Load Peptide files
dat.pep.A549tot = readtable([path , 'tA549_TotalFrac/peptides.txt']); 5
dat.pep.Verotot = readtable([path , 'tVero_TotalFrac/peptides.txt']); 6
dat.pep.A549enr = readtable([path , 'tA549_EnrichVar/peptides.txt']); 7
dat.pep.Veroenr = readtable([path , 'tVero_EnrichVar/peptides.txt']); 8

% Sample names
samples = {'A549tot','Verotot','A549enr','Veroenr'};

% Data cleanup:
% Filters applied to evidence files only as peptides.txt is only used to
% source protein start/end positions on filtered peptides
% Remove contaminants, reverse hits, filter data on PEP, PIF

for ii = 1:numel(samples)
  % Filter on PEP
  dat.evi.(samples{ii}) = dat.evi.(samples{ii})(dat.evi.(samples{ii}).PEP <= 0.02 , :);
   
  % Filer on PIF 
  dat.evi.(samples{ii}) = dat.evi.(samples{ii})(dat.evi.(samples{ii}).PIF >= 0.7 , :);
  
  % Filter on Reverse database hits
  dat.evi.(samples{ii}) = dat.evi.(samples{ii})(categorical(dat.evi.(samples{ii}).Reverse) ~= '+' , :);  

  % Filter on Common contaminants.
  dat.evi.(samples{ii}) = dat.evi.(samples{ii})(categorical(dat.evi.(samples{ii}).PotentialContaminant) ~= '+' , :);

  % Filter data on viral peptides only (All contain 'SARS2')
  dat.evi.(samples{ii}) = dat.evi.(samples{ii})(contains(dat.evi.(samples{ii}).Proteins , 'SARS2') == 1 , :);
  
end

% Format table for output
%%

% Combine ACE2-A549 and Vero data.
combTab = table();

for ii = 1:numel(samples)
   
    tempTab = table();
    
   tempTab.Sequence         = dat.evi.(samples{ii}).Sequence;
   tempTab.ModifiedSequence = dat.evi.(samples{ii}).ModifiedSequence;
   tempTab.Z                = dat.evi.(samples{ii}).Charge;
   tempTab.RT               = dat.evi.(samples{ii}).RetentionTime;
   tempTab.RawFile          = dat.evi.(samples{ii}).RawFile;
   tempTab.ScanNumber       = dat.evi.(samples{ii}).MS_MSScanNumber;
   tempTab.Proteins         = dat.evi.(samples{ii}).Proteins;
   tempTab.Score            = dat.evi.(samples{ii}).Score;
   tempTab.PEP              = dat.evi.(samples{ii}).PEP;
   
   combTab = [combTab ; tempTab];
   
   clear tempTab;
    
end

% Sort on score, then select unique
[~ , idx] = sort(combTab.Score , 'descend');
combTab = combTab(idx,:); 

% Select unique peptides based on the modified sequece. Matlab 'unique' 
% returns the first occurance of each which will now be the highest scoring.
[~ , ia , ~] = unique(combTab.ModifiedSequence);
combTab = combTab(ia , :);

% Fill in the Observed in ACE2A549 and Vero columns
A549Peps = [dat.evi.A549enr.ModifiedSequence ; dat.evi.A549tot.ModifiedSequence];
VeroPeps = [dat.evi.Veroenr.ModifiedSequence ; dat.evi.Verotot.ModifiedSequence];

combTab.ObsInACE2A549(ismember(combTab.ModifiedSequence , A549Peps) == 1) = '+';
combTab.ObsInVeroE6  (ismember(combTab.ModifiedSequence , VeroPeps) == 1) = '+';

% Add Enzyme specificity (Trypsin, Arg-C)
combTab.EnzymeSpecificitySearched(contains(combTab.RawFile , 'Enriched') == 1) = 'A';
combTab.EnzymeSpecificitySearched(contains(combTab.RawFile , 'Total'   ) == 1) = 'T';

% Add start and end positions of peptides using the peptides file
% first generate a combined table containing modified sequences, start
% and end position 
pepTab = table();

for ii = 1:numel(samples)
    tempTab = table();
    
    tempTab.Seq   = dat.pep.(samples{ii}).Sequence;
    tempTab.Start = dat.pep.(samples{ii}).StartPosition;
    tempTab.End   = dat.pep.(samples{ii}).EndPosition;
   
    pepTab = [pepTab ; tempTab];
    clear tempTab;
end

% Map start and end positions to peptides using the peptides file
[~ , locb] = ismember(combTab.Sequence , pepTab.Seq);

combTab.StartPosition = pepTab.Start(locb);
combTab.EndPosition   = pepTab.End(locb);

% Tidy by sorting based on amino acid sequence
[~ , idx] = sort(combTab.Sequence);
combTab = combTab(idx , :);


% Write table
writetable(combTab , [path , '/Tables/ViralPeptides.csv']);



clear A549Peps combTab dat ia idx ii locb pepTab samples VeroPeps