%This code analyzes the sequences downloaded from ncbi with search details
%"poliovirus"[Organism]

function preprocessing_seqs_PV(inputfile)

run startup.m

[header,seqs] = fastaread(inputfile);

%% Filtering the sequences 

indices_seqs_poliovirus = [];
indices_seqs_human_poliovirus = [];
indices_seqs_human_poliovirus1 = [];
indices_seqs_human_poliovirus2 = [];
indices_seqs_human_poliovirus3 = [];

for kk = 1:length(header)
    
    if ~isempty(strfind(header{kk},'poliovirus')) || ...
            ~isempty(strfind(header{kk},'Poliovirus'))
        indices_seqs_poliovirus = ...
            [indices_seqs_poliovirus kk];
    end
    
    if ~isempty(strfind(header{kk},'Human poliovirus')) || ...
            ~isempty(strfind(header{kk},'Human Poliovirus')) || ...
            ~isempty(strfind(header{kk},'human poliovirus')) || ...
            ~isempty(strfind(header{kk},'human Poliovirus'))
        indices_seqs_human_poliovirus = ...
            [indices_seqs_human_poliovirus kk];
    end
    
    if ~isempty(strfind(header{kk},'Human poliovirus 1')) || ...
            ~isempty(strfind(header{kk},'Human Poliovirus 1')) || ...
            ~isempty(strfind(header{kk},'human poliovirus 1')) || ...
            ~isempty(strfind(header{kk},'human Poliovirus 1'))
        indices_seqs_human_poliovirus1 = ...
            [indices_seqs_human_poliovirus1 kk];
    elseif ~isempty(strfind(header{kk},'Human poliovirus 2')) || ...
            ~isempty(strfind(header{kk},'Human Poliovirus 2')) || ...
            ~isempty(strfind(header{kk},'human poliovirus 2')) || ...
            ~isempty(strfind(header{kk},'human Poliovirus 2'))
        indices_seqs_human_poliovirus2 = ...
            [indices_seqs_human_poliovirus2 kk];
    elseif ~isempty(strfind(header{kk},'Human poliovirus 3')) || ...
            ~isempty(strfind(header{kk},'Human Poliovirus 3')) || ...
            ~isempty(strfind(header{kk},'human poliovirus 3')) || ...
            ~isempty(strfind(header{kk},'human Poliovirus 3'))
        indices_seqs_human_poliovirus3 = ...
            [indices_seqs_human_poliovirus3 kk];
    end
    
end

Nseq_total = length(seqs);
Nseq_poliovirus = length(indices_seqs_poliovirus);
Nseq_human_poliovirus = length(indices_seqs_human_poliovirus);
Nseq_human_poliovirus1 = length(indices_seqs_human_poliovirus1);
Nseq_human_poliovirus2 = length(indices_seqs_human_poliovirus2);
Nseq_human_poliovirus3 = length(indices_seqs_human_poliovirus3);

figure(1);
bar(1:3,[Nseq_human_poliovirus1, Nseq_human_poliovirus2, Nseq_human_poliovirus3], 0.3, ...
    'FaceColor',darkgray,'EdgeColor',darkgray)
xlabel('Serotype')
ylabel('Number of sequences')
grid off
axis([0 4 0 3500])
post_proc_fig

%% Separating sequences of human poliovirus 1,2,and 3, and making fasta
%% files

for kk = 1:Nseq_human_poliovirus1
    Human_poliovirus1(kk).Header = header{indices_seqs_human_poliovirus1(kk)};
    Human_poliovirus1(kk).Sequence = seqs{indices_seqs_human_poliovirus1(kk)};
end

% delete Human_poliovirus1.fasta
fastawrite('Human_poliovirus1.fasta',Human_poliovirus1)
   
for kk = 1:Nseq_human_poliovirus2
    Human_poliovirus2(kk).Header = header{indices_seqs_human_poliovirus2(kk)};
    Human_poliovirus2(kk).Sequence = seqs{indices_seqs_human_poliovirus2(kk)};
end

% delete Human_poliovirus2.fasta
fastawrite('Human_poliovirus2.fasta',Human_poliovirus2)

for kk = 1:Nseq_human_poliovirus3
    Human_poliovirus3(kk).Header = header{indices_seqs_human_poliovirus3(kk)};
    Human_poliovirus3(kk).Sequence = seqs{indices_seqs_human_poliovirus3(kk)};
end

% delete Human_poliovirus3.fasta
fastawrite('Human_poliovirus3.fasta',Human_poliovirus3)
