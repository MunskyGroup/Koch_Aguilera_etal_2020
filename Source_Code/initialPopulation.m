function [initPop] = initialPopulation (Model_type)
mainFolder = 'Best_Fits_For_Each_Model';
mainName = 'Fit_OP_MT_';
fn =  [mainFolder,'/' ,mainName,num2str(Model_type),'/','op_val.mat'];
if ispc
    fn =  [mainFolder,'\' ,mainName,num2str(Model_type),'\','op_val.mat'];
end
load(fn);
initPop = x;

end
