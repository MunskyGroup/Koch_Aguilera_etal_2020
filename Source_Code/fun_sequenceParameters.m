function sequenceParameters = fun_sequenceParameters
%% This function stores all the information about the DNA sequence.
% It transforms the DNA sequence into a vector of ratios between tRNA_copyNumber/mean_tRNACopyNumber
% The function stores the tag positions for CAP and IRES genes.
% It calculates codon_time. Defined as sum(1./ratio_ke).

sequenceParameters.ke_ires = sequenceAnalyzer('Sequence_IRES_Kif18b.txt'); % this is a vector containing the ratio between the specific tRNA_copyNumber (i)/ mean_tRNACopyNumber for each codon
sequenceParameters.ke_cap = sequenceAnalyzer('Sequence_CAP_KDM5B.txt'); % this is a vector containing the ratio between the specific tRNA_copyNumber (i)/ mean_tRNACopyNumber for each codon
sequenceParameters.probePosition_CAP =[2,11,20,196,206,218,228,300,309,318];
sequenceParameters.probePosition_IRES = [9,33,57,81,105,129,153,177,201,225,249,273,297,321,345,369,393,417,441,465,489,513,537,561];
sequenceParameters.codon_time_CAP = sum(1./sequenceParameters.ke_cap);
sequenceParameters.codon_time_IRES = sum(1./sequenceParameters.ke_ires);
end
