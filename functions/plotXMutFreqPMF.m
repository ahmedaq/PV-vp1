
%%
function [freqCountNumMutPerSeq numSeq] = plotXMutFreqPMF(mutMtxIn, protLen)

numMutPerSeq = sum(mutMtxIn, 2);

numSeq = size(mutMtxIn,1);

freqCountNumMutPerSeq = zeros(1,protLen);

for i = 1:protLen+1
   freqCountNumMutPerSeq(i) = sum(numMutPerSeq == i-1) ;
    
end

return;