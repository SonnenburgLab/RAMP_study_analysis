library('phangorn')
library('dada2')
library('DECIPHER')

output_path = "/home/mmcarter/user_data/Projects/RAMP/16S/03.output"
save_path = "/home/mmcarter/user_data/Projects/RAMP/16S/03.output"

seqtab.nochim = readRDS(file.path(output_path, "210831_ramp_seqtab.rds"))

# seqs <- getSequences(seqtab.nochim)
# names(seqs) <- seqs # This propagates to the tip labels of the tree
# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# saveRDS(alignment, file = file.path(save_path, "210831_ramp_alignment.rds"))

alignment = readRDS(file = file.path(save_path, "210831_ramp_alignment.rds"))

# print('alignment saved')

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# dm <- dist.ml(phang.align)
# saveRDS(dm, file = file.path(save_path, "210831_ramp_phang.align.dm.rds"))

# print('dm saved')

dm = readRDS(file = file.path(save_path, "210831_ramp_phang.align.dm.rds"))

n <- NJ(dm) # Note, tip order != sequence order

fit = pml(n, data=phang.align)
saveRDS(fit, file = file.path(save_path, "210831_ramp_fit.phang.align.dm.rds"))
## negative edges length changed to 0!
# 

print('fit saved')

fitGTR <- update(fit, k=4, inv=0.2)

saveRDS(fitGTR, file = file.path(save_path, "210831_ketomed_fitGTR.rds"))

print('fitGTR saved')

fitGTR_optim <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR_optim, file = file.path(save_path, "210831_ramp_fitGTR_optim.rds"))

