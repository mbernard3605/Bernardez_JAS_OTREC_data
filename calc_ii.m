function [ii] = calc_ii(satent,z)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



ii=mean(satent(z<=3000 & z>=1000,:,:))-mean(satent(z<=7000 & z>=5000,:,:));


end

