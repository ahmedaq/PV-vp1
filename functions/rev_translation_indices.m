function Out_sector = rev_translation_indices(In_sector,V)

Out_sector = [];

for mm = 1:length(In_sector)
    for kk = 1:length(V)
        if V(kk)==In_sector(mm)
            Out_sector(mm)=kk;
        end
    end
end
    


