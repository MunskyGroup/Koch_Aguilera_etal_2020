

%% Frequency of codon usage. http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
strGenCopy = struct( ...
'TTT', 17.6, 'TCT', 15.2, 'TAT', 12.2,  'TGT', 10.6, 'TTC', 20.3,...  
'TCC', 17.7, 'TAC', 15.3, 'TGC', 12.6, 'TTA', 7.7,   'TCA', 12.2,...  
'TAA', 1.0, 'TGA', 1.6, 'TTG', 12.9,  'TCG',  4.4,  'TAG', 0.8,...  
'TGG', 13.2,'CTT', 13.2, 'CCT', 17.5, 'CAT', 10.9,  'CGT', 4.5,...
'CTC', 19.6, 'CCC', 19.8, 'CAC', 15.1,  'CGC', 10.4, 'CTA',  7.2,... 
'CCA', 16.9, 'CAA', 12.3, 'CGA',  6.2,  'CTG', 39.6, 'CCG',  6.9,...
'CAG', 34.2, 'CGG', 11.4, 'ATT', 16.0,  'ACT', 13.1, 'AAT', 17.0,...  
'AGT', 12.1, 'ATC', 20.8, 'ACC', 18.9, 'AAC', 19.1,  'AGC', 19.5,...
'ATA',  7.5,  'ACA', 15.1,  'AAA', 24.4,  'AGA', 12.2, 'ATG', 22.0, ... 
'ACG', 6.1, 'AAG', 31.9,  'AGG', 12.0, 'GTT', 11.0,  'GCT', 18.4,...  
'GAT', 21.8, 'GGT', 10.8,'GTC', 14.5,  'GCC', 27.7,  'GAC', 25.1,...  
'GGC', 22.2, 'GTA',  7.1, 'GCA', 15.8,  'GAA', 29.0,  'GGA', 16.5,...
'GTG', 28.1,  'GCG', 7.4, 'GAG', 39.6, 'GGG', 16.5);


% T = struct2table(strGenCopy);


strGenCopy_Fast = struct( ...
'GCT', 27.7,  'GCC', 27.7,'GCA', 27.7,  'GCG', 27.7, ... % A
'CGT', 12.2, 'CGC', 12.2, 'CGA', 12.2, 'CGG', 12.2,'AGA', 12.2, 'AGG', 12.2,  ... % R
'AAT', 19.1,'AAC', 19.1,  ... % N 
'GAT', 25.1,'GAC', 25.1, ...  % D
'TGT', 12.6,  'TGC', 12.6, ... % C
'CAA', 34.2, 'CAG', 34.2, ... % Q
'GAA', 39.6, 'GAG', 39.6, ... % E
'GGT', 22.2,  'GGC', 22.2, 'GGA', 22.2, 'GGG', 22.2, ... % G
'CAT', 15.1, 'CAC', 15.1, ... % H
'ATT', 20.8, 'ATC', 20.8, 'ATA', 20.8, ... % I
'TTA', 39.6, 'TTG', 39.6,'CTT', 39.6,'CTC', 39.6, 'CTA', 39.6, 'CTG', 39.6,... % L
'AAA', 31.9, 'AAG', 31.9, ... % K
'ATG', 22.0,...   % M
'TTT', 20.3, 'TTC', 20.3, ...   % F
'CCT', 19.8,  'CCC', 19.8, 'CCA', 19.8, 'CCG', 19.8, ... % P
'TCT', 19.5, 'TCC', 19.5,'TCA', 19.5, 'TCG',  19.5, 'AGT', 19.5,'AGC', 19.5, ... % S
'ACT', 18.9, 'ACC', 18.9,  'ACA', 18.9, 'ACG', 18.9,... % T
'TGG', 13.2,  ... % W
'TAT', 15.3, 'TAC', 15.3, ... % Y
'GTT', 28.1,'GTC', 28.1, 'GTA', 28.1, 'GTG', 28.1, ... % V
'TAA', 1.6, 'TAG', 1.6, 'TGA', 1.6); % STOP

strGenCopy_Slow = struct( ...
'GCT', 7.4,  'GCC', 7.4,'GCA', 7.4,  'GCG', 7.4, ... % A
'CGT', 4.5, 'CGC', 4.5, 'CGA', 4.5, 'CGG', 4.5,'AGA', 4.5, 'AGG', 4.5,  ... %R
'AAT', 17.0,'AAC', 17.0,  ... %N
'GAT', 21.8,'GAC', 21.8, ... %D
'TGT', 10.6,  'TGC', 10.6, ... %C
'CAA', 12.3, 'CAG', 12.3, ... %Q
'GAA', 29.0, 'GAG', 29.0, ... % E
'GGT', 10.8,  'GGC', 10.8, 'GGA', 10.8, 'GGG', 10.8, ... %G
'CAT', 10.9, 'CAC', 10.9, ... %H
'ATT', 7.5, 'ATC', 7.5, 'ATA', 7.5, ... %I
'TTA', 7.2, 'TTG', 7.2,'CTT', 7.2,'CTC', 7.2, 'CTA', 7.2, 'CTG', 7.2,... %L
'AAA', 24.4, 'AAG', 24.4, ... %K
'ATG', 22.0,... %M
'TTT', 17.6, 'TTC', 17.6, ... % F
'CCT', 6.9,  'CCC', 6.9, 'CCA', 6.9, 'CCG',  6.9, ... %P
'TCT', 4.4, 'TCC', 4.4,'TCA', 4.4, 'TCG', 4.4, 'AGT', 4.4,'AGC', 4.4, ... %S
'ACT', 6.1, 'ACC', 6.1,  'ACA', 6.1, 'ACG', 6.1,...  %T
'TGG', 13.2,  ... % W
'TAT', 12.2, 'TAC', 12.2, ... % Y
'GTT', 7.1,'GTC', 7.1, 'GTA', 7.1, 'GTG', 7.1, ... % V
'TAA', 0.8, 'TAG', 0.8,  'TGA', 0.8);  % STOP CODON


ratio_Sentitivity_Fast_Slow = ([strGenCopy_Fast.GCT / strGenCopy_Slow.GCT, ... %A
strGenCopy_Fast.CGT / strGenCopy_Slow.CGT, ... %R
strGenCopy_Fast.AAT / strGenCopy_Slow.AAT, ... %N
strGenCopy_Fast.GAT / strGenCopy_Slow.GAT, ... %D
strGenCopy_Fast.TGT / strGenCopy_Slow.TGT, ... %C
strGenCopy_Fast.CAA / strGenCopy_Slow.CAA, ... %Q
strGenCopy_Fast.GAA / strGenCopy_Slow.GAA, ... %E
strGenCopy_Fast.GGT / strGenCopy_Slow.GGT, ... %G
strGenCopy_Fast.CAT / strGenCopy_Slow.CAT, ... %H
strGenCopy_Fast.ATT / strGenCopy_Slow.ATT, ... %I
strGenCopy_Fast.TTA / strGenCopy_Slow.TTA, ... %L
strGenCopy_Fast.AAA / strGenCopy_Slow.AAA, ... %K
strGenCopy_Fast.ATG / strGenCopy_Slow.ATG, ... %M
strGenCopy_Fast.TTT / strGenCopy_Slow.TTT, ... %F
strGenCopy_Fast.CCT / strGenCopy_Slow.CCT, ... %P
strGenCopy_Fast.TCT / strGenCopy_Slow.TCT, ... %S
strGenCopy_Fast.ACT / strGenCopy_Slow.ACT, ... %T
strGenCopy_Fast.TGG / strGenCopy_Slow.TGG, ... %W
strGenCopy_Fast.TAT / strGenCopy_Slow.TAT, ... %Y
strGenCopy_Fast.GTT / strGenCopy_Slow.GTT]); %V

