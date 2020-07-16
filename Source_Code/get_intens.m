% function that converts ribosome 
function Intens = get_intens(dt,dtmax,t_vs_position)
dt = dt(dt>0&dt<dtmax);
A = repmat(t_vs_position,length(dt),1);  % times at which probes would bind.
B = repmat(dt,1,length(t_vs_position));  % current times since initiation
Intens = sum(B(:)>A(:));  % sum of probes that have been passed but where the protein has not completed.
end
