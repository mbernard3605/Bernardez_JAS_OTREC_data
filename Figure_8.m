%Script that plots the base of figure 8

close all;
clear variables;
clc;

%Set parameters of toy model

RH=.70;
dtdz=-1.4/200;
z=0:200:16000;
z=z(1:78);
T=zeros(length(z),1);T1=T;T2=T; 
T(1)=300;
load('pres.mat');

%Calculate temperature profile
for i = 2:length(z)
T(i)=T(i-1)+dtdz*(z(i)-z(i-1));    
end
T=T(1:78);


%%


qvstar=calc_qvstar(T,presmean*100);

qv=RH*qvstar;
qv(1:5)=(RH+.1)*qvstar(1:5);

[s,h,hsat,hp,tp]= zero_plume_buoyancy_multi_plume(T,qv,z',presmean*100,4000);
tdiff=tp-T;
ztop=hp>=hsat;
tpert=-exp(-(z-5000).^2./2000.^2);
T4=T+tpert';
qvstar4=calc_qvstar(T4,presmean*100);

qv4=RH*qvstar4;
qv4(1:5)=(RH+.1)*qvstar4(1:5);

[s4,h4,hsat4,hp4,tp4]= zero_plume_buoyancy_multi_plume(T4,qv,z',presmean*100,4000);
tdiff4=tp4-T4;
ztop4=hp4>=hsat4;

%Set the levels and parameters for the anomalies
     c3=2000;
     c2=1000;
     c1=1000;
     %a1=500/(c1*sqrt(2*pi));
     %a2=1000/(c2*sqrt(2*pi));
     b1=4000;
     b2=6000;
     b3=5000;

    
     %Calculate moisture anomalies that are added, both anomalies are
     %gaussians and the MDC is constructed to remove and add equal moisture
     qpert1 = exp(-(z-b1).^2./c1.^2);
     a1=(.98-RH)/max(qpert1);
     qpert1=qpert1.*a1;
     qvpert1=qpert1'.*qvstar;
     qanom1 = trapz(presmean,-qvpert1/9.81);
     qpert2 = exp(-(z-b2).^2./c2.^2);
     a2=qanom1./trapz(presmean,qpert2'.*qvstar/9.81);
     qpert2 = a2*qpert2;
     qvpert2=qpert2'.*qvstar;
     qanom2 = trapz(presmean,-qvpert2/9.81);
     
      qpert3 = exp(-(z-b3).^2./c3.^2);
     a3=(.9-RH)/max(qpert3);
     qpert3=qpert3.*a3;
     qvpert3=qpert3'.*qvstar;


    

     MDC_pos_anom=qvpert1+qvpert2;
     SF_pos_anom=qvpert3;

     
     
%Add the moisture anomalies to the base moisture profile
qv2=qv+MDC_pos_anom;
qv3=qv+SF_pos_anom;

%Calculate the plume buoyancies
[s2,h2,hsat2,hp2,tp2]= zero_plume_buoyancy_multi_plume(T,qv2,z',presmean*100,4000);
ztop2=hp2>=hsat2;
tdiff2=tp2-T;
[s3,h3,hsat3,hp3,tp3]= zero_plume_buoyancy_multi_plume(T,qv3,z',presmean*100,4000);
ztop3=hp3>=hsat3;
tdiff3=tp3-T;

%Find the 800 hPa and 400 hPa levels to display
plot_level1=find(presmean<=800,1);
plot_level2=find(presmean<=400,1);

%%
% plot figure 8
figure('units','normalized','outerposition',[0.0375 0 0.75 .95]);


subplot(2,3,1),
plot(h3,presmean,'b',hsat3,presmean,'r',hp3(ztop3),presmean(ztop3),'k-.','linewidth',3);
hold on;

p=plot(h,presmean,'b',hsat,presmean,'r',hp(ztop),presmean(ztop),'k','linewidth',1);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','fontsize',16,'xlim',[3.16e5 3.42e5],'yticklabels',{},'xticklabels',{});
title('High SF');
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
yline(presmean(find(presmean<=600,1)),'--','linewidth',1.5);
ylabel('Pressure');

subplot(2,3,2)
p2=plot(tdiff,presmean,'k',tdiff3,presmean,'k-.','linewidth',3);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','xlim',[0 10],'fontsize',16,'yticklabels',{},'xticklabels',{});
title('Vertical Motion proxy');
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
yline(presmean(find(presmean<=600,1)),'--','linewidth',1.5);
legend(p,'MSE','Saturated MSE','Plume MSE','location','northwest');
legend(p2,'Default','Anom. added')


subplot(2,3,3)
p2=plot(diff(tdiff)./tdiff(1:end-1)-diff(tdiff3)./tdiff3(1:end-1),presmean(1:end-1),'k-.','linewidth',3);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','fontsize',16,'yticklabels',{},'xticklabels',{});
hold on;
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
plot([-15 15],[presmean(find(presmean<=600,1)) presmean(find(presmean<=600,1))],'--k','linewidth',1.5);
plot(zeros(length(presmean),1),presmean);
inBetween1=[diff(tdiff)./tdiff(1:end-1)-diff(tdiff3)./tdiff3(1:end-1) zeros(length(tdiff)-1,1)];
inBetween1(isnan(inBetween1))=0;
p1=[presmean(1:end-1) flipud(presmean(1:end-1))];
fill(real(inBetween1),p1,'b');
p3(1)=plot([0 0],[presmean(1) presmean(1)],'b','linewidth',10);
set(gca,'xlim',[-.075 .075]);
title('detrainment difference')
xlabel('$\Delta \frac{1}{w} \frac{\partial w}{\partial p} = -\Delta \delta$', 'interpreter', 'latex','fontsize',24);



subplot(2,3,4),
plot(h2,presmean,'b',hsat2,presmean,'r',hp2(ztop2),presmean(ztop2),'k--','linewidth',3);
hold on;

p=plot(h,presmean,'b',hsat,presmean,'r',hp(ztop),presmean(ztop),'k','linewidth',1);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','fontsize',16,'xlim',[3.16e5 3.42e5],'yticklabels',{},'xticklabels',{});
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
yline(presmean(find(presmean<=600,1)),'--','linewidth',1.5);
xlabel('Energy J/kg');
title('High MDC');
ylabel('Pressure');

subplot(2,3,5)
p2=plot(tdiff,presmean,'k',tdiff2,presmean,'k--','linewidth',3);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','xlim',[0 10],'fontsize',16,'yticklabels',{},'xticklabels',{});
xlabel('Plume buoyancy','fontweight','bold');
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
yline(presmean(find(presmean<=600,1)),'--','linewidth',1.5);
legend(p2,'Default','Anom. added')

subplot(2,3,6)
p2=plot(diff(tdiff)./tdiff(1:end-1)-diff(tdiff2)./tdiff2(1:end-1),presmean(1:end-1),'k--','linewidth',3);
set(gca,'ylim',[presmean(plot_level2) presmean(plot_level1)],'ydir','reverse','fontsize',16,'yticklabels',{},'xticklabels',{});
hold on;
yline(presmean(plot_level1),'linewidth',1.5);
yline(presmean(plot_level2),'linewidth',1.5);
plot([-15 15],[presmean(find(presmean<=600,1)) presmean(find(presmean<=600,1))],'--k','linewidth',1.5);
plot(zeros(length(presmean),1),presmean);
inBetween2=[diff(tdiff)./tdiff(1:end-1)-diff(tdiff2)./tdiff2(1:end-1) zeros(length(tdiff)-1,1)];
inBetween2(isnan(inBetween2))=0;
p1=[presmean(1:end-1) flipud(presmean(1:end-1))];
fill(real(inBetween2),p1,'b');
p3(1)=plot([0 0],[presmean(1) presmean(1)],'b','linewidth',10);
set(gca,'xlim',[-.075 .075]);
xlabel('       <-Drying Moistening->');
print(gcf,'-dpng','Schematic_one_II_temp.png');


