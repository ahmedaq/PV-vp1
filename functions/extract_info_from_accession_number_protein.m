function [title,journal,country,year,gi] = ...
    extract_info_from_accession_number_protein(accession_number)

%Input: A cell array consisting of accession numbers

%Output: title of paper, journal, country, and collection year of the
%strain associated with the accession number

%% Getting complete info from genbank using accession number

title = cell(length(accession_number),1);
journal = cell(length(accession_number),1);
country = cell(length(accession_number),1);
% gi = zeros(length(accession_number),1);
% year = zeros(length(accession_number),4);
for kk = 1:length(accession_number)
%     clear data
    data = getgenpept(accession_number{kk}); %If the accession number are of protein strains
%     data = getgenbank(accession_number{kk}); %If the accession number are of nucleotide strains
%     gi(kk) = data.GI;
    title{kk} = data.Reference{1}.Title;
    journal{kk} = data.Reference{1}.Journal;
    datafeatures = data.Features;
    row_country = zeros(1,size(datafeatures,1));
    row_date = zeros(1,size(datafeatures,1));
    for mm = 1:size(datafeatures,1)
        row_country(mm) = ~isempty(strfind(datafeatures(mm,:),'country'));
        row_date(mm) = ~isempty(strfind(datafeatures(mm,:),'collection_date'));
    end
    row_country_value = find(row_country==1);
    row_date_value = find(row_date==1);
    country_index_quotes = find(datafeatures(row_country_value,:)=='"');
    if isempty(country_index_quotes)
        country{kk} = [];
    else
        country{kk} = datafeatures(row_country_value,country_index_quotes(1)+1:country_index_quotes(2)-1);
    end
    collectiondate_index_quotes = find(datafeatures(row_date_value,:)=='"');
    if isempty(collectiondate_index_quotes)
        year(kk,:) = '0000';
    else
        if length(collectiondate_index_quotes(1):collectiondate_index_quotes(2)) == 6
            year(kk,:) = datafeatures(row_date_value,collectiondate_index_quotes(1)+1:collectiondate_index_quotes(2)-1);
        else
            year(kk,:) = datafeatures(row_date_value,collectiondate_index_quotes(2)-4:collectiondate_index_quotes(2)-1);
        end
    end
    kk
end

