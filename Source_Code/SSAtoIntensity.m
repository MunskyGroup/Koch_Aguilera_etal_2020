% function [ IntensityVectors_1G,IntensityVectors_2G,pre_numOfRibosomes_1G,pre_numOfRibosomes_2G] = SSAtoIntensity(RibsomePositions,probePosition_1G,probePosition_2G,k_elongation_1G )
function [IntensityVectors_1G,IntensityVectors_2G] = SSAtoIntensity(RibsomePositions,sequenceParameters)

probePosition_1G =sequenceParameters.probePosition_CAP;
probePosition_2G = sequenceParameters.probePosition_IRES;
k_elongation_1G = sequenceParameters.ke_cap;

geneLength_1G = length(k_elongation_1G);
numberOftimePoints_withoutBurningTime = size(RibsomePositions,2);

%% Saving ribosome position in time
X_output_1G = RibsomePositions;
X_output_1G (X_output_1G>geneLength_1G) = 0;
X_output_2G = RibsomePositions;
X_output_2G (X_output_2G<geneLength_1G) = 0;
X_output_2G = (X_output_2G)- geneLength_1G;
X_output_2G (X_output_2G<0) = 0;

%% Saving the intensity vectors in a cell array.
IntensityVectors_1G = zeros (1,numberOftimePoints_withoutBurningTime );
IntensityVectors_2G  = zeros (1,numberOftimePoints_withoutBurningTime );

for tp =1: numberOftimePoints_withoutBurningTime
    IntensityVectors_1G(tp) = sum(sum(X_output_1G(:,tp)>=probePosition_1G));
    IntensityVectors_2G(tp) = sum(sum(X_output_2G(:,tp)>=probePosition_2G));
end

%% Maxing Intensity vectors a row vector of zeros in case the vector is zero
if numberOftimePoints_withoutBurningTime>1
    if length(IntensityVectors_1G)==1
        IntensityVectors_1G = zeros(1,numberOftimePoints_withoutBurningTime);
    end
    if  length(IntensityVectors_2G)==1
        IntensityVectors_2G = zeros(1,numberOftimePoints_withoutBurningTime);
    end
    IntensityVectors_1G = [IntensityVectors_1G./length(probePosition_1G)]';
    IntensityVectors_2G = [IntensityVectors_2G./length(probePosition_2G)]';
end

end

