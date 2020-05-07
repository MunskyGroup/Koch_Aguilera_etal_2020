function [sim_CAP_ONLY, sim_IRES_ONLY, error_CAP_ONLY, error_IRES_ONLY] = intesities_selection_HT(intensityVector_1G, intensityVector_2G,pointOfNormalization)
nTimePoints=size(intensityVector_1G,2); % or 50 exp repetitions
nRepetitions = size(intensityVector_1G,1);
int_1Gm = zeros (nRepetitions,nTimePoints);
int_2Gm = zeros (nRepetitions,nTimePoints);
for i =1: nRepetitions
    for j =1: nTimePoints
        if intensityVector_1G(i,j) >0 %&& intensityVector_2G(i,j) == 0
            int_1Gm (i,j) = intensityVector_1G(i,j);
        end
        if  intensityVector_2G(i,j) >0
            int_2Gm (i,j) = intensityVector_2G(i,j);
        end
    end
end
int_CAP_ONLY = sum(int_1Gm);
int_IRES_ONLY = sum(int_2Gm);
n1 =  mean(mean(int_CAP_ONLY(:,1:pointOfNormalization)));
n2 =  mean(mean(int_IRES_ONLY(:,1:pointOfNormalization)));
sim_CAP_ONLY = int_CAP_ONLY./n1;
sim_IRES_ONLY  = int_IRES_ONLY./n2;

%% calculating error by bootstraping
nbootstr =1e3;
nSelection =1000;
int_1Gm_C{nbootstr} =[];
int_2Gm_C{nbootstr} =[];
parfor kk =1: nbootstr
    index_BS =    randi([1 nRepetitions],[1, nSelection]);
    intensityVector_1G_bs = intensityVector_1G(index_BS,:);
    intensityVector_2G_bs = intensityVector_2G(index_BS,:);
    int_1Gm_C{kk} = zeros (nSelection,nTimePoints);
    int_2Gm_C{kk} = zeros (nSelection,nTimePoints);
    for i =1: nSelection
        for j =1: nTimePoints
            if intensityVector_1G_bs(i,j) >0 %&& intensityVector_2G_bs(i,j) == 0
                int_1Gm_C{kk}(i,j) = intensityVector_1G_bs(i,j);
            end
            if  intensityVector_2G_bs(i,j) >0
                int_2Gm_C{kk} (i,j) = intensityVector_2G_bs(i,j);
            end
        end
    end
    int_CAP_ONLY = sum(int_1Gm_C{kk});
    int_IRES_ONLY = sum(int_2Gm_C{kk});
    n1 =  mean(mean(int_CAP_ONLY(:,1:pointOfNormalization)));
    n2 =  mean(mean(int_IRES_ONLY(:,1:pointOfNormalization)));
    sim_CAP_ONLY_bs(kk,:) = int_CAP_ONLY./n1;
    sim_IRES_ONLY_bs(kk,:)  = int_IRES_ONLY./n2;
end
error_CAP_ONLY = std(sim_CAP_ONLY_bs);
error_IRES_ONLY = std(sim_IRES_ONLY_bs);