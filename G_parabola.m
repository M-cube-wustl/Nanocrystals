function [ G] = G_parabola( G0, Gf, N, is )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
G = 1*(Gf-G0)/N^2*is.^2+G0;


end
