function [sq,sfb,sft] = calc_sq(qv,qv_sat,pres)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%pressure should be in units of Pa
height_index=find(pres>=60000);


g=9.81;
sft=trapz(pres(height_index(end)+1:end),-qv(height_index(end)+1:end,:)/g)./...
    trapz(pres,-qv_sat/g);
sfb=trapz(pres(1:height_index(end)),-qv(1:height_index(end),:)/g)./...
    trapz(pres,-qv_sat/g);

sq=sfb-sft;

end

