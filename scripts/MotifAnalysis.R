
## Motif analysis script

library(protr)
library(stringr)
library(seqinr)
library(dagLogo)



load('DiffExprAnalysis.RData')

# Using dagLogo to produce overrepresented AAs

## Making background sets of proteins detected in the experiment
Prot_detected_A549 = Prot_A549[which(names(Prot_A549) %in% pep_A549$Leading.razor.protein)]
Prot_detected_vero = Prot_vero[which(names(Prot_vero) %in% pep_vero$Leading.razor.protein)]

prtm_A549 = prepareProteome(fasta = 'Prot_A549_trimmed.fasta')
prtm_vero = prepareProteome(fasta = 'Prot_vero_trimmed.fasta')

## Peptide sets

pep_set_A549 = do.call(rbind, lapply(1:nrow(nntpep_A549f), function(i) {
    c(protein = nntpep_A549f$Leading.razor.protein[i], 
      peptide_seq = nntpep_A549f$Sequence[i],
      #cleav_region = paste0(substr(nntpep_A549f$N.term.cleavage.window[i], 25, 30), 
      #                      substr(nntpep_A549f$C.term.cleavage.window[i], 1, 5)
      cleave_region = substr(Prot_A549[[nntpep_A549f$Leading.razor.protein[i]]],
                             nntpep_A549f$Start.position[i]-5, nntpep_A549f$Start.position[i]+4))
}))
pep_set_A549 = as.data.frame(pep_set_A549)


pep_set_vero = do.call(rbind, lapply(1:nrow(nntpep_verof), function(i) {
    c(protein = nntpep_verof$Leading.razor.protein[i], 
      peptide_seq = nntpep_verof$Sequence[i],
      #cleav_region = paste0(substr(nntpep_verof$N.term.cleavage.window[i], 25, 30), 
      #                      substr(nntpep_verof$C.term.cleavage.window[i], 1, 5)))
      cleave_region = substr(Prot_vero[[nntpep_verof$Leading.razor.protein[i]]],
                             nntpep_verof$Start.position[i]-5, nntpep_verof$Start.position[i]+4))
}))
pep_set_vero = as.data.frame(pep_set_vero)

# remove remaining SARS proteins from the dataset
pep_set_A549 = pep_set_A549[substr(pep_set_A549$protein, 1,2) == 'sp',]
pep_set_vero = pep_set_vero[substr(pep_set_vero$protein, 1,2) == 'tr',]


# Some cleavage sites don't have 5 AAs before so padded by "_"
pep_set_A549_padded = pep_set_A549
pep_set_A549_padded$cleave_region = sapply(pep_set_A549$cleave_region, 
                                           function(x) {
                                               if(nchar(x)<10) paste0(paste0(rep('_', 10-nchar(x)), collapse=''),x)
                                               else x
                                               })

pep_set_vero_padded = pep_set_vero
pep_set_vero_padded$cleave_region = sapply(pep_set_vero$cleave_region, 
                                           function(x) {
                                               if(nchar(x)<10) paste0(paste0(rep('_', 10-nchar(x)), collapse=''),x)
                                               else x
                                               })


seq_obj_A549 = formatSequence(seq = pep_set_A549_padded$cleave_region, proteome = prtm_A549, upstreamOffset = 5, downstreamOffset = 5)
seq_obj_vero = formatSequence(seq = pep_set_vero_padded$cleave_region, proteome = prtm_vero, upstreamOffset = 5, downstreamOffset = 5)

# bakground models
bg_fisher_A549 = buildBackgroundModel(seq_obj_A549, background = "wholeProteome", proteome = prtm_A549, testType = "fisher")
bg_ztest_A549 = buildBackgroundModel(seq_obj_A549, background = "wholeProteome", proteome = prtm_A549, testType = "ztest")

bg_fisher_vero = buildBackgroundModel(seq_obj_vero, background = "wholeProteome", proteome = prtm_vero, testType = "fisher")
bg_ztest_vero = buildBackgroundModel(seq_obj_vero, background = "wholeProteome", proteome = prtm_vero, testType = "ztest")

# Statistical tests
t0_fisher_A549 = testDAU(seq_obj_A549, dagBackground = bg_fisher_A549)
t0_ztest_A549 = testDAU(seq_obj_A549, dagBackground = bg_ztest_A549)

t0_fisher_vero = testDAU(seq_obj_vero, dagBackground = bg_fisher_vero)
t0_ztest_vero = testDAU(seq_obj_vero, dagBackground = bg_ztest_vero)

# Visualisations

dagHeatmap(t0_fisher_A549)
dagHeatmap(t0_ztest_A549) 

dagHeatmap(t0_fisher_vero) 
dagHeatmap(t0_ztest_vero) 

dagLogo(t0_fisher_A549) 
dagLogo(t0_ztest_A549) 

dagLogo(t0_fisher_vero) 
dagLogo(t0_ztest_vero) 


# Repeating with significant peptides only
## Differential expression analysis (DiffExpression.R) should be run here!
## At this point the NeoNT peptides are assigned with more stringent rules that are used
## in the second half of differential expression analysis.


pep_set_A549sig_padded = pep_set_A549_padded[pep_set_A549_padded$peptide_seq %in% sig_pep_A549$Sequence,]
pep_set_verosig_padded = pep_set_vero_padded[pep_set_vero_padded$peptide_seq %in% sig_pep_vero$Sequence,]


seq_obj_A549sig = formatSequence(seq = pep_set_A549sig_padded$cleave_region, proteome = prtm_A549, upstreamOffset = 5, downstreamOffset = 5)
seq_obj_verosig = formatSequence(seq = pep_set_verosig_padded$cleave_region, proteome = prtm_vero, upstreamOffset = 5, downstreamOffset = 5)

# bakground models
bg_fisher_A549sig = buildBackgroundModel(seq_obj_A549sig, background = "wholeProteome", proteome = prtm_A549, testType = "fisher")
bg_ztest_A549sig = buildBackgroundModel(seq_obj_A549sig, background = "wholeProteome", proteome = prtm_A549, testType = "ztest")

bg_fisher_verosig = buildBackgroundModel(seq_obj_verosig, background = "wholeProteome", proteome = prtm_vero, testType = "fisher")
bg_ztest_verosig = buildBackgroundModel(seq_obj_verosig, background = "wholeProteome", proteome = prtm_vero, testType = "ztest")

# Statistical tests
t0_fisher_A549sig = testDAU(seq_obj_A549sig, dagBackground = bg_fisher_A549)
t0_ztest_A549sig = testDAU(seq_obj_A549sig, dagBackground = bg_ztest_A549)

t0_fisher_verosig = testDAU(seq_obj_verosig, dagBackground = bg_fisher_vero)
t0_ztest_verosig = testDAU(seq_obj_verosig, dagBackground = bg_ztest_vero)

# Visualisations

dagHeatmap(t0_fisher_A549sig) 
dagHeatmap(t0_ztest_A549sig) 

dagHeatmap(t0_fisher_verosig) 
dagHeatmap(t0_ztest_verosig) 

dagLogo(t0_fisher_A549sig) 
dagLogo(t0_ztest_A549sig) 

dagLogo(t0_fisher_verosig) 
dagLogo(t0_ztest_verosig) 


# Repeating the motif analysis without duplicated proteins

pep_set_A549sig_padded_nodupes = pep_set_A549_padded[pep_set_A549_padded$peptide_seq %in% sig_pep_A549$Sequence,]
pep_set_A549sig_padded_nodupes = pep_set_A549sig_padded_nodupes[-duplicated(pep_set_A549sig_padded_nodupes$protein),]

seq_obj_A549sig_nodupes = formatSequence(seq = pep_set_A549sig_padded_nodupes$cleave_region[-1], 
                                 proteome = prtm_A549, 
                                 upstreamOffset = 5, 
                                 downstreamOffset = 5)

bg_fisher_A549sig_noDupes = buildBackgroundModel(seq_obj_A549sig_nodupes, 
                                                 background = "wholeProteome", 
                                                 proteome = prtm_A549, testType = "fisher")
bg_ztest_A549sig_noDupes = buildBackgroundModel(seq_obj_A549sig_nodupes, 
                                                background = "wholeProteome", 
                                                proteome = prtm_A549, testType = "ztest")

# Statistical tests
t0_fisher_A549sig_noDupes = testDAU(seq_obj_A549sig_nodupes, dagBackground = bg_fisher_A549)
t0_ztest_A549sig_noDupes = testDAU(seq_obj_A549sig_nodupes, dagBackground = bg_ztest_A549)


# Visualisations

dagHeatmap(t0_fisher_A549sig_noDupes) 
dagHeatmap(t0_ztest_A549sig_noDupes) 

dagLogo(t0_fisher_A549sig_noDupes) 
dagLogo(t0_ztest_A549sig_noDupes) 

save.image('MotifAnalysis.RData')
# =================================================== 
# ==================== Functions ==================== 
# =================================================== 

annotate_peptide = function(pep, prot, prots, subseq_length=10, NeoNT='Original'){
    offset = round(subseq_length/2)
    
    idx = grep(prot, names(prots), fixed = T)
    if(length(idx)<1) return(c(peptide = pep, 
                               protein = prot, 
                               position = 0, 
                               subseq = '',
                               NeoNT = NeoNT
                               ))
    if (length(idx) >1) {
        print(length(idx))
        print(pep)
        print(prot)
        return(c(peptide = pep, 
                 protein = prot, 
                 position = 0, 
                 subseq = '',
                 NeoNT = NeoNT
        ))
        }
    do.call(rbind, lapply(idx, function(ii){
        prot_seq = prots[[ii]]
        loc = str_locate(prot_seq, pep)
        
        c(peptide = pep, 
          protein = prot, 
          position = loc[1], 
          subseq = substr(prot_seq, loc[1]-offset, loc[1]+offset-1),
          NeoNT = NeoNT
        )
    }))
    
}
