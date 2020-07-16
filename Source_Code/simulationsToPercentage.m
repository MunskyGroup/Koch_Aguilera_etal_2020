function [sim_1G,sim_2G,sim_BG,sim_NT, err_sim_1G,err_sim_2G,err_sim_BG,err_sim_NT]  = simulationsToPercentage (numberOfSpots,StrData,observationTime,plottingCondition)
%% bootstraping
n_bootstrapping = 50;
Sampling_numberOfSpots =300; % representative number of spots per cell

vec_sim_0F = zeros(1,n_bootstrapping);
vec_sim_1F = zeros(1,n_bootstrapping);
vec_sim_BF = zeros(1,n_bootstrapping);
vec_sim_NT = zeros(1,n_bootstrapping);
for j =1: n_bootstrapping
    numberOfTrajectoriesIn1G = 0;
    numberOfTrajectoriesIn2G = 0;
    numberOfTrajectoriesInBG = 0;
    numberOfTrajectoriesNonTrans = 0;
    
    for i =1:Sampling_numberOfSpots
        ii= randi(numberOfSpots);
        if max(StrData(1,ii).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,ii).trajectory_2(end-observationTime:end)) == 0
            numberOfTrajectoriesIn1G = numberOfTrajectoriesIn1G + 1;
        end
        if max(StrData(1,ii).trajectory_2(end-observationTime:end)) > 0  && max(StrData(1,ii).trajectory_1(end-observationTime:end)) == 0
            numberOfTrajectoriesIn2G = numberOfTrajectoriesIn2G + 1;
        end
        if max(StrData(1,ii).trajectory_1(end-observationTime:end)) > 0 && max(StrData(1,ii).trajectory_2(end-observationTime:end)) > 0
            numberOfTrajectoriesInBG = numberOfTrajectoriesInBG + 1;
        end
        if max(StrData(1,ii).trajectory_1(end-observationTime:end)) == 0 && max(StrData(1,ii).trajectory_2(end-observationTime:end)) == 0
            numberOfTrajectoriesNonTrans = numberOfTrajectoriesNonTrans + 1;
        end
    end
    vec_sim_0F(j) = (numberOfTrajectoriesIn1G./Sampling_numberOfSpots)*100;
    vec_sim_1F(j) = (numberOfTrajectoriesIn2G./Sampling_numberOfSpots)*100;
    vec_sim_BF(j) = (numberOfTrajectoriesInBG./Sampling_numberOfSpots)*100;
    vec_sim_NT(j) = (numberOfTrajectoriesNonTrans./Sampling_numberOfSpots)*100;
end

if plottingCondition ==1
    
end

err_sim_1G = std(vec_sim_0F);
err_sim_2G = std(vec_sim_1F);
err_sim_BG = std(vec_sim_BF);
err_sim_NT = std(vec_sim_NT);

sim_1G = mean(vec_sim_0F);
sim_2G = mean(vec_sim_1F);
sim_BG = mean(vec_sim_BF);
sim_NT = mean(vec_sim_NT);

end