% Thresholds the correlation matrix to retain only positive correlations, setting all other values to zero.
% inputs : data = correlation matrix with weights
% output : matrix with only positive correlations,
function [pos_wmatrix]=posweights(data)
pos_wmatrix=zeros(size (data,1), size(data,2));
for j=1:(size (data, 1)*size(data,2));
    if data(j)>0;
        pos_wmatrix(j)=data(j);
    else data(j)<=0;
        pos_wmatrix(j)=0;
    end
end