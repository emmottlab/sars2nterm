# Processing protein/peptide data for motif analysis

load('ProcAnnot.Rdata')

Prot_A549 = readFASTA('../data/HomoSapiens_spONLY_20350_dl20200508.fasta')
Prot_vero = readFASTA('../data/GreenMonkey_19223_dl20200516_UP000029965.fasta')
Prot_SARS = readFASTA('../data/SARS2_Custom_20200518.fasta')

table(dat$A549$pep$Amino.acid.before, dat$A549$pep$NeoNT)
table(dat$A549$pep$First.amino.acid, dat$A549$pep$NeoNT)
table(dat$A549$pep$Last.amino.acid, dat$A549$pep$NeoNT)
table(dat$A549$pep$Amino.acid.after, dat$A549$pep$NeoNT)

table(dat$Vero$pep$Amino.acid.before, dat$Vero$pep$NeoNT)
table(dat$Vero$pep$First.amino.acid, dat$Vero$pep$NeoNT)
table(dat$Vero$pep$Last.amino.acid, dat$Vero$pep$NeoNT)
table(dat$Vero$pep$Amino.acid.after, dat$Vero$pep$NeoNT)

# Further filtering
nntpep_A549 = dat$A549$pep[dat$A549$pep$NeoNT == 'NeoNT' & is.na(dat$A549$annot$signalpep),]
nntpep_vero = dat$Vero$pep[dat$Vero$pep$NeoNT == 'NeoNT' & is.na(dat$Vero$annot$signalpep),]

## Removing sequences that have R/K as AA one before the start
nntpep_A549 = nntpep_A549[!(nntpep_A549$Amino.acid.before %in% c('R','K')),]
nntpep_vero = nntpep_vero[!(nntpep_vero$Amino.acid.before %in% c('R','K')),]

# Predicting signal peptides
## Trimming the fasta files to only include proteins of remaining peptides

Prot_A549_trimmed = Prot_A549[names(Prot_A549) %in% nntpep_A549$Leading.razor.protein]
Prot_vero_trimmed = Prot_vero[names(Prot_vero) %in% nntpep_vero$Leading.razor.protein]

write.fasta(Prot_A549_trimmed, names=names(Prot_A549_trimmed), file.out = 'Prot_A549_trimmed.fasta')
write.fasta(Prot_vero_trimmed, names=names(Prot_vero_trimmed), file.out = 'Prot_vero_trimmed.fasta')

## Using http://www.cbs.dtu.dk/services/SignalP/ for prediction of signal peptides

SigPep_A549 = read.table('output_protein_type_A549.txt', sep='\t')
SigPep_vero = read.table('output_protein_type_vero.txt', sep='\t')

SigPep_A549 = SigPep_A549[SigPep_A549$V2 != 'OTHER',]
SigPep_A549 = SigPep_A549[order(SigPep_A549$V3, decreasing = T),]
SigPep_A549$pos = as.numeric(substr(SigPep_A549$V5, 12,13))

SigPep_vero = SigPep_vero[SigPep_vero$V2 != 'OTHER',]
SigPep_vero = SigPep_vero[order(SigPep_vero$V3, decreasing = T),]
SigPep_vero$pos = as.numeric(substr(SigPep_vero$V5, 12,13))

#Matching up the predicted signal peptides to the peptide table
nntpep_A549$signalPepPos = SigPep_A549$pos[match(gsub('[|]', '_', nntpep_A549$Leading.razor.protein), SigPep_A549$V1)]
nntpep_A549$signalPep = nntpep_A549$signalPepPos - nntpep_A549$Start.position

nntpep_A549$SigPosStr = SigPep_A549$V5[match(gsub('[|]', '_', nntpep_A549$Leading.razor.protein), SigPep_A549$V1)]
nntpep_A549$SigConf = SigPep_A549$V3[match(gsub('[|]', '_', nntpep_A549$Leading.razor.protein), SigPep_A549$V1)]
nntpep_A549$PosConf = substr(nntpep_A549$SigPosStr,27,34)


nntpep_vero$signalPepPos = SigPep_vero$pos[match(gsub('[|]', '_', nntpep_vero$Leading.razor.protein), SigPep_vero$V1)]
nntpep_vero$signalPep = nntpep_vero$signalPepPos - nntpep_vero$Start.position

nntpep_vero$SigPosStr = SigPep_vero$V5[match(gsub('[|]', '_', nntpep_vero$Leading.razor.protein), SigPep_vero$V1)]
nntpep_vero$SigConf = SigPep_vero$V3[match(gsub('[|]', '_', nntpep_vero$Leading.razor.protein), SigPep_vero$V1)]
nntpep_vero$PosConf = substr(nntpep_vero$SigPosStr,27,34)

nntpep_A549$NeoNT[abs(nntpep_A549$signalPep) <= 5 & nntpep_A549$SigConf > 0.9 & nntpep_A549$PosConf < 90] = 'Signal'
nntpep_vero$NeoNT[abs(nntpep_vero$signalPep) <= 5 & nntpep_vero$SigConf > 0.9 & nntpep_vero$PosConf < 90] = 'Signal'

table(nntpep_A549$NeoNT)
table(nntpep_vero$NeoNT)

nntpep_A549f = nntpep_A549[nntpep_A549$NeoNT=='NeoNT',]
nntpep_verof = nntpep_vero[nntpep_vero$NeoNT=='NeoNT',]

save.image('ProtFiltered.RData')