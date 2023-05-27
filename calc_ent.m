function [ent] = calc_ent(T,qv,pres)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Tf=273.15;
pr=1e5;
Rd=287;
cp=1004;
L=2.5e6;

theta=T.*(pr./pres).^(Rd/cp);
sd=cp*log(theta/Tf);

ent=sd+L*qv/Tf;

end

