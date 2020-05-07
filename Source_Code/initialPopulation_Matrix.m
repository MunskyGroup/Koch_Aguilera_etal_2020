function [initPop] = initialPopulation_Matrix (Model_type)
mainFolder = 'Best_Fits_For_Each_Model';
mainName = 'Fit_OP_MT_';
All_init_Pops = [];
for Model_type = 1:14
    fn =  [mainFolder,'/' ,mainName,num2str(Model_type),'/','op_val.mat'];
    if ispc
        fn =  [mainFolder,'\' ,mainName,num2str(Model_type),'\','op_val.mat'];
    end
    load(fn);
    initPop = x;
    All_init_Pops = [All_init_Pops;initPop];
end
initPop = All_init_Pops;
end
