function [sf] = calc_sf(qv,qv_sat,pres)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

g=9.81;
sf=trapz(pres,-qv/g)./trapz(pres,-qv_sat/g);

end

