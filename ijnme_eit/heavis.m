function y=heavis(x)
%-------------------------------------------------------------------%
% compute the heavisde function: y=1 if x\geq0, otherwise x=0;      %
% JIN Bangti (kimbtsing@yahoo.com.cn), Feb. 2, 2009                 %
%-------------------------------------------------------------------%

y=ones(size(x));
y(find(x<0))=0;

return