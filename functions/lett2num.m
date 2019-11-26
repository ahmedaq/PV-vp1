function msa_num=lett2num(msa_lett)
% usage: msa_num=lett2num(msa_lett)
%
% Translates an alignment from a representation where the 20 natural amino
% acids are represented by letters to a representation where they are
% represented by the numbers 1,...,20. Any symbol not corresponding to an
% amino acid is represented by B = 21.
%
% Author: Olivier Rivoire
% 2/2010
%
%**************************************************************************

code='ACDEFGHIKLMNPQRSTVWYB-';

msa_num=zeros(size(msa_lett));

for s=1:size(msa_lett,1)
    for a=1:size(msa_lett,2)
        lett=msa_lett(s,a);
        for index=1:length(code)
            if lett==code(index)
                msa_num(s,a)=index;
            end
        end
    end
end

