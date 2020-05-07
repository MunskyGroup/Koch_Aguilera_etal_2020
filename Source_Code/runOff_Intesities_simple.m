function [sim_CAP_ONLY, sim_IRES_ONLY, error_CAP_ONLY, error_IRES_ONLY] = runOff_Intesities_simple(intensityVector_1G, intensityVector_2G)
nSpots=size(intensityVector_1G,2); % or 50 exp repetitions
nRepetitions_ForError = size(intensityVector_1G,1);
norm1= mean (intensityVector_1G(:,1));
norm2 = mean (intensityVector_2G(:,1));
% Normalize vectores if they are non zero.
if norm1>0 && norm2>0
    IntensityVectors_1G_HT = zeros(nRepetitions_ForError,nSpots);
    IntensityVectors_2G_HT = zeros(nRepetitions_ForError,nSpots);
    for i =1: size(intensityVector_1G,1)
        IntensityVectors_1G_HT(i,:) =intensityVector_1G(i,:) ./ norm1;
        IntensityVectors_2G_HT(i,:) = intensityVector_2G(i,:) ./ norm2;
    end
    %% Calculating mean values
    int1 =  mean(IntensityVectors_1G_HT);
    int2 = mean(IntensityVectors_2G_HT);
    
    %% second normalization with respect to first 6 points
    sim_CAP_ONLY = int1./ mean(int1(1:6)) ;
    sim_IRES_ONLY = int2./ mean(int2(1:6)) ;
    
    %% Calculating SEM
    er_int1 = std(IntensityVectors_1G_HT)./ sqrt(nRepetitions_ForError);
    er_int2 = std(IntensityVectors_2G_HT)./ sqrt(nRepetitions_ForError);
    error_CAP_ONLY =  er_int1;
    error_IRES_ONLY = er_int2;
else
    % Prealocate vectores to zeros if the means are zero.
    sim_CAP_ONLY = zeros(1,nSpots);
    sim_IRES_ONLY = zeros(1,nSpots);
    error_CAP_ONLY = zeros(1,nSpots);
    error_IRES_ONLY = zeros(1,nSpots);
end
end