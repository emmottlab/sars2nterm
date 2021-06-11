
## A script for generating figures and inputs for TopFind

library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(grImport2)

load('Final.RData')

## 1. TopFind

A549_sig_topfind = datF_ff$A549$pep[datF_ff$A549$annot$NeoNT=='NeoNT',][res_A549$q.values <= 0.05,]
Vero_sig_topfind = datF_ff$Vero$pep[datF_ff$Vero$annot$NeoNT=='NeoNT',][res_vero$q.values <= 0.05,]

A549_topfind_input = A549_sig_topfind[,c('up_id','Sequence')]
write.csv(A549_topfind_input, 'A549_topfind_input.csv')

Vero_topfind_input = Vero_sig_topfind[,c('up_id','Sequence')]
write.csv(Vero_topfind_input, 'Vero_topfind_input.csv')


A549_sig_topfind_up = datF_ff$A549$pep[datF_ff$A549$annot$NeoNT=='NeoNT',][res_A549$q.values <= 0.05 & res_A549$logFC > 0,]
A549_sig_topfind_down = datF_ff$A549$pep[datF_ff$A549$annot$NeoNT=='NeoNT',][res_A549$q.values <= 0.05 & res_A549$logFC < 0,]

A549_topfind_input_up = A549_sig_topfind_up[,c('up_id','Sequence')]
write.csv(A549_topfind_input_up, 'A549_topfind_input_up.csv')

A549_topfind_input_down = A549_sig_topfind_down[,c('up_id','Sequence')]
write.csv(A549_topfind_input_down, 'A549_topfind_input_down.csv')


## 2. Can I have the motif analysis as a .svg as it is 
# (two panels side by side, with larger font sizes for the X/Y ticks.

svglite::svglite('Figures/Motifs_A549.svg')
dagLogo(t0_fisher_A549sig, fontsize = 10, title='A549')
dev.off()

svglite::svglite('Figures/Motifs_Vero.svg')
dagLogo(t0_fisher_verosig, fontsize = 10, title='Vero') 
dev.off()


## 3. Can you do me the gene ratio/p-adjust plots as a 3x2 (or 4x2) set of panels.
# I note that Vero there wasn’t a biological process slide unlike the A549 
# - if this was because it wasn’t significant, please just leave an empty spot where the panel would go.
# Have A549 as the left set, and Vero as the right hand set.



p1_1 = dotplot(enr1$A549$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, showCategory=20)+ 
    ggtitle('Transcirption factor targets')+
    theme(plot.title = element_text(hjust = 0.5))

p1_2 = dotplot(enr1$A549$neopeps_24h$c5.bp.v7.1.entrez.gmt, showCategory=20)+
    ggtitle('GO: Biological process')+
    theme(plot.title = element_text(hjust = 0.5))

p1_3 = dotplot(enr1$A549$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)+
    ggtitle('GO: Celular compartment')+
    theme(plot.title = element_text(hjust = 0.5))

tmp = enr1$A549$neopeps_24h$c5.mf.v7.1.entrez.gmt
tmp@result$p.adjust = round(tmp@result$p.adjust,3)
p1_4 = dotplot(tmp, showCategory=20)+ 
    ggtitle('GO: Molecular function')+
    theme(plot.title = element_text(hjust = 0.5))


p2_1 = dotplot(enr1$Vero$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, showCategory=20)+
    ggtitle('Transcirption factor targets')+
    theme(plot.title = element_text(hjust = 0.5))

p2_2 = ggplot(NULL)+theme_minimal()+ 
    ggtitle('GO: Biological process (No significant terms)')+
    theme(plot.title = element_text(hjust = 0.5))
    
p2_3 = dotplot(enr1$Vero$neopeps_24h$c5.cc.v7.1.entrez.gmt, showCategory=20)+ 
    ggtitle('GO: Celular compartment')+
    theme(plot.title = element_text(hjust = 0.5))

p2_4 = dotplot(enr1$Vero$neopeps_24h$c5.mf.v7.1.entrez.gmt, showCategory=20)+ 
    ggtitle('GO: Molecular function')+
    theme(plot.title = element_text(hjust = 0.5))

svglite::svglite('Figures/EnrichmentPways.svg', width = 20, height=30)
ggarrange(ggarrange(p1_1, p1_2, p1_3, p1_4, nrow=4, ncol=1),
          ggarrange(p2_1, p2_2, p2_3, p2_4, nrow=4, ncol=1),
          ncol=2, nrow=1, labels=c('A549', 'Vero'))
dev.off()

## 4. Can you please do me a similarly laid-out figure with the spider diagram panels
# that are paired with the above plots. These are currently numbered. 
# If the individual points could be labeled with ‘GENENAME_XXX’ where xxx is the position
# of the first amino acid of the Neo-N-terminus.


# Replace entrez_ids with gene names + peptide start pos

replace_genePepPos = function(enr_obj, datF, is.vero=F){
    idx = which(enr_obj@result$qvalue <=0.05)
    genes = strsplit(enr_obj@result$geneID[idx], '/')
    
    if(is.vero) entrez_ids = datF$annot$entrezid.hs
    else entrez_ids = datF$annot$entrezid
    
    replaced = lapply(genes, function(lst){
        idxs = match(lst, entrez_ids)
        paste(datF$annot$rep.gene[idxs],
              '-',
              datF$pep$Start.position[idxs],
              collapse='/', sep='')
    })
    reps = unlist(replaced)
    enr_obj@result$geneID[idx] = reps
    enr_obj
}



# Make plots
q1_1 = cnetplot(replace_genePepPos(enr1$A549$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, datF_ff$A549), showCategory=20)+ 
    ggtitle('Transcirption factor targets')+
    theme(plot.title = element_text(hjust = 0.5))

q1_2 = cnetplot(replace_genePepPos(enr1$A549$neopeps_24h$c5.bp.v7.1.entrez.gmt, datF_ff$A549), showCategory=20)+
    ggtitle('GO: Biological process')+
    theme(plot.title = element_text(hjust = 0.5))

q1_3 = cnetplot(replace_genePepPos(enr1$A549$neopeps_24h$c5.cc.v7.1.entrez.gmt, datF_ff$A549), showCategory=20)+
    ggtitle('GO: Celular compartment')+
    theme(plot.title = element_text(hjust = 0.5))

q1_4 = cnetplot(replace_genePepPos(enr1$A549$neopeps_24h$c5.mf.v7.1.entrez.gmt, datF_ff$A549), showCategory=20)+ 
    ggtitle('GO: Molecular function')+
    theme(plot.title = element_text(hjust = 0.5))


q2_1 = cnetplot(replace_genePepPos(enr1$Vero$neopeps_24h$c3.tft.gtrd.v7.1.entrez.gmt, datF_ff$Vero, is.vero=T), showCategory=20)+
    ggtitle('Transcirption factor targets')+
    theme(plot.title = element_text(hjust = 0.5))

q2_2 = ggplot(NULL)+theme_minimal()+ 
    ggtitle('GO: Biological process (No significant terms)')+
    theme(plot.title = element_text(hjust = 0.5))
    
q2_3 = cnetplot(replace_genePepPos(enr1$Vero$neopeps_24h$c5.cc.v7.1.entrez.gmt, datF_ff$Vero, is.vero=T), showCategory=20)+ 
    ggtitle('GO: Celular compartment')+
    theme(plot.title = element_text(hjust = 0.5))

q2_4 = cnetplot(replace_genePepPos(enr1$Vero$neopeps_24h$c5.mf.v7.1.entrez.gmt, datF_ff$Vero, is.vero=T), showCategory=20)+ 
    ggtitle('GO: Molecular function')+
    theme(plot.title = element_text(hjust = 0.5))


svglite::svglite('Figures/EnrichmentNetworks.svg', width = 20, height=30)
ggarrange(ggarrange(q1_1, q1_2, q1_3, q1_4, nrow=4, ncol=1),
          ggarrange(q2_1, q2_2, q2_3, q2_4, nrow=4, ncol=1),
          ncol=2, nrow=1, labels=c('A549', 'Vero'))
dev.off()





