function [parametersModel ] = sequenceAnalyzer(input_file)
%% This function takes the gene sequence and returns the ratio between the specific tRNA_copyNumber (i)/ mean_tRNACopyNumber for each codon

try
    whole_gene = fastaread(input_file);
catch
    whole_gene.Sequence =input_file;
end

typeOfTag = [];
tagPositions = [];
geneSequence = [];

%% Defining the distinct tag sequences
T_SunTag = 'EELLSKNYHLENEVARLKK';
T_Flag = 'DYKDDDDK';
T_Hemagglutinin = 'YPYDVPDYA';

try
    %% Detecting all the possible ORFs in the gene
    orfPositions = seqshoworfs(whole_gene.Sequence,'MinimumLength',100, 'nodisplay', 'true','Frames',1);
    
    %% Read each possible ORF and trasnlate those sequences from DNA to protein
    counter = 1;
    for i =1: length (orfPositions)
        for j =1: length (orfPositions(i).Start)
            try
                temporal_ORF{counter} = whole_gene.Sequence (orfPositions(i).Start(j):orfPositions(i).Stop(j)-1);
                temporal_Proteins {counter} =   nt2aa(temporal_ORF{counter});
                counter = counter+1;
            catch
            end
        end
    end
    
    %% Find if one of those proteins contains the tags.
    
    position_SunTag = strfind(temporal_Proteins,T_SunTag);
    position_Flag = strfind(temporal_Proteins,T_Flag);
    position_Hemagglutinin = strfind(temporal_Proteins,T_Hemagglutinin);
    
    %# find empty cells in the cell array
    emptyCells_SunTag = cellfun(@isempty,position_SunTag);
    emptyCells_Flag = cellfun(@isempty,position_Flag);
    emptyCells_Hemagglutinin = cellfun(@isempty,position_Hemagglutinin);
    
    % selecting the Number of Protein that contains the Tag region
    selected_counter_SunTag = find(~emptyCells_SunTag,1);
    selected_counter_Flag = find(~emptyCells_Flag,1);
    selected_counter_Hemagglutinin = find(~emptyCells_Hemagglutinin,1);
    
    %# remove empty cells
    position_SunTag(emptyCells_SunTag) = [];
    position_Flag(emptyCells_Flag) = [];
    position_Hemagglutinin(emptyCells_Hemagglutinin) = [];
    
    
    %% Finding the type of Tags
    if isempty(position_SunTag) ==0
        typeOfTag = 'SunTag';
        tagPositions = position_SunTag{1};
        geneSequence = upper(temporal_ORF{selected_counter_SunTag});
    end
    
    if isempty(position_Flag) == 0
        typeOfTag = 'Flag';
        tagPositions = position_Flag{1};
        geneSequence = upper(temporal_ORF{selected_counter_Flag});
    end
    
    if isempty(position_Hemagglutinin) ==0
        typeOfTag = 'HA';
        tagPositions = position_Hemagglutinin{1};
        geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
    end
    
    if isempty(position_Hemagglutinin) ==0 && isempty(position_Flag) == 0
        typeOfTag = 'HA_Flag';
        tagPositions = [position_Hemagglutinin{1}];
        geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
    end
    
    if isempty(position_Hemagglutinin) ==0 && isempty(position_SunTag) ==0
        typeOfTag = 'HA_SunTag';
        tagPositions = [position_SunTag{1}];
        geneSequence = upper(temporal_ORF{selected_counter_Hemagglutinin});
    end
    
    generated_geneSequence.Sequence= geneSequence;
catch
    
    
end
codons = geneSequence;
geneLength = length(codons)/3;
generated_geneSequence.Sequence= geneSequence;

%% Elongation constant.
k_elongation=  zeros (1,geneLength-1 ); % this removes the stop codon.
genomicCopyNumber;
tRNA_copyNumber = zeros (1,geneLength );

%%  Separating sequence in codons.
counter = 1;
for i =1 : geneLength
    separated_codons(i,:) = codons(counter : counter + 2);
    counter = counter + 3;
end

for i = 1 : geneLength
    field = separated_codons(i,:);
    tRNA_copyNumber (i) = strGenCopy.(field) ;
end
% calculating the mean values in the structure
mean_tRNACopyNumber = mean (mean(reshape(struct2array(strGenCopy),numel(fieldnames(strGenCopy)),[]),2));
for i = 1 : geneLength-1
    k_elongation(i)= (tRNA_copyNumber (i)/ mean_tRNACopyNumber);
end
parametersModel = k_elongation;
end

