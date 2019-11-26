function Pmatrix = peaks_geographical_analysis(input_data_file,data_peaks_rank,data_country_vp1,data_year_vp1)

load(input_data_file)
load(data_peaks_rank)
load(data_country_vp1)
load(data_year_vp1)

run startup.m

%%
for kk = 1:no_of_peak_seqs_vp1
    year_seqs_peak_seq{kk} = year_vp1(indices_seqs_that_converged_to_peak_seq_vp1{kk});
    country_seqs_peak_seq{kk} = country(indices_seqs_that_converged_to_peak_seq_vp1{kk});
end

%Adding years of valley 10 manually based on the associated report
year_seqs_peak_seq{10}(1:17) = 1992;
year_seqs_peak_seq{10}(18:35) = 1993;

clear no_unique_countries
for kk = 1:no_of_peak_seqs_vp1
    unique_countries{kk} = unique(country_seqs_peak_seq{kk});
    for mm = 1:length(unique_countries{kk})
        temp = 0;        
        for nn = 1:length(country_seqs_peak_seq{kk})            
            if strcmp(unique_countries{kk}{mm},country_seqs_peak_seq{kk}{nn})==1
                temp = temp + 1;
                year_country{kk}{mm}(temp) = year_seqs_peak_seq{kk}(nn);                
            end
        end
        no_unique_countries{kk}(mm) = temp;
    end
end

%% Distribution of the number of sequences in the ten most-populous peaks
for kk = 1:10 %top 10 peaks
    Pmatrix{kk} = cell(length(no_unique_countries{kk}),2);
    for mm = 1:length(no_unique_countries{kk})
        Pmatrix{kk}{mm,1} = unique_countries{kk}{mm};
        Pmatrix{kk}{mm,2} = no_unique_countries{kk}(mm);
    end
end