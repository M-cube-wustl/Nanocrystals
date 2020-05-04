function [  ] = write_NC_2_xlsx( fn, B_WT,B_ST, krange, WTebs, STebs, WTmeds, STmeds )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = numel(krange);
B_WT_vec = [B_WT, zeros(1,N-2)];
B_ST_vec = [B_ST, zeros(1,N-2)];
A = [krange', WTebs', STebs', WTmeds', STmeds', B_WT_vec', B_ST_vec'];
Acell = num2cell(A);
hcell = {'k','WTeb','STeb','WTmed','STmed','B_WT','B_ST'};
A = [hcell;Acell]
range = ['A1:G',num2str(N+1)];
xlswrite(fn,A,1,range)

end

