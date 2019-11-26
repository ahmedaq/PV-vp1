function E = calcSeqEnergy(seq,H)

E = seq*triu(-H)*seq';