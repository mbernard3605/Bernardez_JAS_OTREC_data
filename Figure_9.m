clear variables;
%close all;
clc;

%load the EOFS from ERA5
load('ERA5_EOFS.mat');
weights_old=weights;
lds_ERA=lds(:,1:2)./repmat(weights_old.^0.5,1,2);


load('anglecolormap2.mat');
testmap=[anglemap(33:end,:); anglemap(1:32,:)];
colors=linspace(-180,180,64);
g=9.81;
Rd=287;
cp=1005;
Lv=2.5e6;
%load Predict data





load('OTREC_B1.mat');
ang_ERA_qv_binned_OTREC_B1 = ang_ERA;
o1_ERA_qv_binned_OTREC_B1 = o1_ERA;
o2_ERA_qv_binned_OTREC_B1 = o2_ERA;
DCIN_qv_binned_OTREC_B1 = DCIN;
II_qv_binned_OTREC_B1 = iiShape;
SF_qv_binned_OTREC_B1 = sfShape;
SFt_qv_binned_OTREC_B1 = sftShape;
SFb_qv_binned_OTREC_B1 = sfbShape;
Tb_qv_binned_OTREC_B1 = Tb;
sst_qv_binned_OTREC_B1 = sstShape;
omega_binn_OTREC_B1 = omegaShapeFull;
pres_B1=presShape;


load('OTREC_B2.mat');
ang_ERA_qv_binned_OTREC_B2 = ang_ERA;
o1_ERA_qv_binned_OTREC_B2 = o1_ERA;
o2_ERA_qv_binned_OTREC_B2 = o2_ERA;
DCIN_qv_binned_OTREC_B2 = DCIN;
II_qv_binned_OTREC_B2 = iiShape;
SF_qv_binned_OTREC_B2 = sfShape;
SFt_qv_binned_OTREC_B2 = sftShape;
SFb_qv_binned_OTREC_B2 = sfbShape;
Tb_qv_binned_OTREC_B2 = Tb;
sst_qv_binned_OTREC_B2 = sstShape;
omega_binn_OTREC_B2 = omegaShapeFull;
pres_B2=presShape;

load('OTREC_B1_zbp_start_ent.mat');
Temp_OTREC_B1=tempShape;
MR_OTREC_B1=mrShape;
z_OTREC_B1=z;
p_OTREC_B1=presShape;

load('OTREC_B2_zbp_start_ent.mat');
Temp_OTREC_B2=tempShape;
MR_OTREC_B2=mrShape;
z_OTREC_B2=z;
p_OTREC_B2=presShape;

% 
% ang_ERA=[ang_ERA_PREDICT ang_ERA_qv_binned_OTREC_B1 ang_ERA_qv_binned_OTREC_B2];
% O1=[o1_ERA_PREDICT o1_ERA_qv_binned_OTREC_B1 o1_ERA_qv_binned_OTREC_B2];
% O2=[o2_ERA_PREDICT o2_ERA_qv_binned_OTREC_B1 o2_ERA_qv_binned_OTREC_B2];
% DCIN=[DCIN_PREDICT DCIN_qv_binned_OTREC_B1 DCIN_qv_binned_OTREC_B2];
% II=[II_PREDICT II_qv_binned_OTREC_B1 II_qv_binned_OTREC_B2];
% SF=[SF_PREDICT SF_qv_binned_OTREC_B1 SF_qv_binned_OTREC_B2];
% SQ=[SFb_PREDICT./SF_PREDICT SFb_qv_binned_OTREC_B1./SF_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2./SF_qv_binned_OTREC_B2]*100;
% SQ2=[SFb_PREDICT-SFt_PREDICT SFb_qv_binned_OTREC_B1-SFt_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2-SFt_qv_binned_OTREC_B2];
% SQ3=[(SFb_PREDICT-SFt_PREDICT)./SF_PREDICT (SFb_qv_binned_OTREC_B1-SFt_qv_binned_OTREC_B1)./SF_qv_binned_OTREC_B1 (SFb_qv_binned_OTREC_B2-SFt_qv_binned_OTREC_B2)./SF_qv_binned_OTREC_B2];
% SFt=[SFt_PREDICT SFt_qv_binned_OTREC_B1 SFt_qv_binned_OTREC_B2];
% SFb=[SFb_PREDICT SFb_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2];
% Tb=[Tb_PREDICT Tb_qv_binned_OTREC_B1 Tb_qv_binned_OTREC_B2];
% sst=[sst_PREDICT sst_qv_binned_OTREC_B1 sst_qv_binned_OTREC_B2];
% %sst=sst-nanmean(sst);
% sst_anom2=[sst_PREDICT-nanmean(sst_PREDICT) sst_qv_binned_OTREC_B1-nanmean(sst_qv_binned_OTREC_B1) sst_qv_binned_OTREC_B2-nanmean(sst_qv_binned_OTREC_B2)];



ang_ERA=[ang_ERA_qv_binned_OTREC_B1 ang_ERA_qv_binned_OTREC_B2];
O1_ERA=[o1_ERA_qv_binned_OTREC_B1 o1_ERA_qv_binned_OTREC_B2];
O2_ERA=[o2_ERA_qv_binned_OTREC_B1 o2_ERA_qv_binned_OTREC_B2];
DCIN=[DCIN_qv_binned_OTREC_B1 DCIN_qv_binned_OTREC_B2];
II=[II_qv_binned_OTREC_B1 II_qv_binned_OTREC_B2];
SF=[SF_qv_binned_OTREC_B1 SF_qv_binned_OTREC_B2];
SQ=[SFb_qv_binned_OTREC_B1./SF_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2./SF_qv_binned_OTREC_B2]*100;
MDC=[SFb_qv_binned_OTREC_B1-SFt_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2-SFt_qv_binned_OTREC_B2];
SQ3=[(SFb_qv_binned_OTREC_B1-SFt_qv_binned_OTREC_B1)./SF_qv_binned_OTREC_B1 (SFb_qv_binned_OTREC_B2-SFt_qv_binned_OTREC_B2)./SF_qv_binned_OTREC_B2];
SFt=[SFt_qv_binned_OTREC_B1 SFt_qv_binned_OTREC_B2];
SFb=[SFb_qv_binned_OTREC_B1 SFb_qv_binned_OTREC_B2];
Tb=[Tb_qv_binned_OTREC_B1 Tb_qv_binned_OTREC_B2];
SST=[sst_qv_binned_OTREC_B1 sst_qv_binned_OTREC_B2];
sst_anom2=[sst_qv_binned_OTREC_B1-nanmean(sst_qv_binned_OTREC_B1) sst_qv_binned_OTREC_B2-nanmean(sst_qv_binned_OTREC_B2)];
omega=[omega_binn_OTREC_B1 omega_binn_OTREC_B2];
pres=[pres_B1 pres_B2]/100;
Temp=[Temp_OTREC_B1 Temp_OTREC_B2];
MR=[MR_OTREC_B1 MR_OTREC_B2];
z=[z_OTREC_B1 z_OTREC_B2];
z=nanmean(z(1:78,:),2)';
angles=linspace(-180,180,size(testmap,1));
colors_plot=interp1(angles,testmap,ang_ERA,[],'extrap');
% 
% tp=zeros(size(Temp));hp=tp;zlcl=zeros(1,size(Temp,2));
% for i = 1:size(tempShape,2)
%     [tp(:,i),hp(:,i),zlcl(i)]=zero_plume_buoyancy_split(Temp(:,i),MR(:,i),z(1:78)'*1000,pres(:,i)*100);
% end
%%
%plot the plume buoyancies based upon a mean temperature and moisture
%profile
presmean=nanmean(pres,2)*100;
cwv=trapz(presmean,-MR/g);

satqv=calc_qvstar(Temp,presmean);






presmean=nanmean(pres,2)*100;
cwv=trapz(presmean,-MR/g);
satMR=calc_qvstar(Temp,presmean);


SFbin1=quantile(SF,2/3);
%cwvbin=quantile(cwv,0.75);
MDC_bin1=quantile(MDC,2/3);
MDC_bin2=quantile(MDC,1/3);
T1=nanmean(Temp(:,(MDC>MDC_bin1 & SF>SFbin1)),2);
T2=nanmean(Temp(:,((MDC>MDC_bin2 & MDC<=MDC_bin1) & SF>SFbin1)),2);
qvstar2=calc_qvstar(T2,presmean);

T3=nanmean(Temp(:,(MDC<=MDC_bin2 & SF>SFbin1)),2);



qv1=nanmean(MR(:,(MDC>MDC_bin1 & SF>SFbin1)),2);
qv2=nanmean(MR(:,((MDC>MDC_bin2 & MDC<=MDC_bin1) & SF>SFbin1)),2);
qv3=nanmean(MR(:,(MDC<=MDC_bin2 & SF>SFbin1)),2);
qvstar1=calc_qvstar(T1,presmean);
qvstar3=calc_qvstar(T3,presmean);


omega1=nanmean(omega(:,(MDC>MDC_bin1 & SF>SFbin1)),2);
omega2=nanmean(omega(:,((MDC>MDC_bin2 & MDC<=MDC_bin1) & SF>SFbin1)),2);
omega3=nanmean(omega(:,(MDC<=MDC_bin2 & SF>SFbin1)),2);
omega_highq=nanmean(omega(:,SF>SFbin1),2);
highq_T=nanmean(Temp(:,SF>SFbin1),2);

highq_q=nanmean(MR(:,SF>SFbin1),2);
highq_satq=calc_qvstar(highq_T,presmean);
highq_rh=highq_q./highq_satq;
meant=nanmean(Temp,2);
meanq=nanmean(MR,2);
meansatq=calc_qvstar(meant,presmean);
meanrh=meanq./meansatq;


[smean,hmean,hsatmean,hpmean,tpmean]= zero_plume_buoyancy_multi_plume(highq_T,highq_q,z(1:78)'*1000,presmean,2000);
bmean=tpmean-highq_T;
qtenmean=diff(bmean);

[s1,h1,hsat1,hp1_mean,tp1_mean]= zero_plume_buoyancy_multi_plume(T1,qv1,z(1:78)'*1000,presmean,2000);
b1=tp1_mean-T1;
qten1=diff(b1);
[s2,h2,hsat2,hp2_mean,tp2_mean]= zero_plume_buoyancy_multi_plume(T2,qv2,z(1:78)'*1000,presmean,2000);
b2=tp2_mean-T2;
qten2=diff(b2);
[s3,h3,hsat3,hp3_mean,tp3_mean]= zero_plume_buoyancy_multi_plume(T3,qv3,z(1:78)'*1000,presmean,2000);
b3=tp3_mean-T3;
qten3=diff(b3);
[s4,h4,hsat4,hp4_mean,tp4_mean]= zero_plume_buoyancy_multi_plume(T1,qvstar1.*highq_rh,z(1:78)'*1000,presmean,2000);
b4=tp4_mean-T1;
qten4=diff(b4);
[s5,h5,hsat5,hp5_mean,tp5_mean]= zero_plume_buoyancy_multi_plume(T2,qvstar2.*highq_rh,z(1:78)'*1000,presmean,2000);
b5=tp5_mean-T2;
qten5=diff(b5);
[s6,h6,hsat6,hp6_mean,tp6_mean]= zero_plume_buoyancy_multi_plume(T3,qvstar3.*highq_rh,z(1:78)'*1000,presmean,2000);
b6=tp6_mean-T3;
qten6=diff(b6);
[s7,h7,hsat7,hp7_mean,tp7_mean]= zero_plume_buoyancy_multi_plume(highq_T,qv1,z(1:78)'*1000,presmean,2000);
b7=tp7_mean-highq_T;
qten7=diff(b7);
[s8,h8,hsat8,hp8_mean,tp8_mean]= zero_plume_buoyancy_multi_plume(highq_T,qv2,z(1:78)'*1000,presmean,2000);
b8=tp8_mean-highq_T;
qten8=diff(b8);
[s9,h9,hsat9,hp9_mean,tp9_mean]= zero_plume_buoyancy_multi_plume(highq_T,qv3,z(1:78)'*1000,presmean,2000);
b9=tp9_mean-highq_T;
qten9=diff(b9);

%%


figure('units','normalized','outerposition',[0.0375 0 0.78125 .675]);
subplot(1,4,1),
plot(omega1,presmean,'r--',omega2,presmean,'k-',omega3,presmean,'B-.','linewidth',3);
set(gca,'ylim',[20000 max(nanmean(presShape,2))],'ydir','reverse','xdir','reverse','fontsize',16,'xlim',[-0.5 0.01],'yticklabel',[200:100:1000]);
ylabel('pressure [hPa]');xlabel('Omega [Pa/s]');
title('Vertical Motion');

subplot(1,4,2),
plot(b1-bmean,presmean,'r--',b2-bmean,presmean,'k-',b3-bmean,presmean,'B-.','linewidth',3);
set(gca,'ylim',[20000 max(nanmean(presShape,2))],'ydir','reverse','xlim',[-0.75 0.75],'fontsize',16,'yticklabel',[]);
xlabel('Buoyancy [K]');
title({'Full Variability'});
legend('High MDC','Med MDC','Low MDC','position',[0.300 0.800 0.010 0.10]);

subplot(1,4,3),
plot(b4-bmean,presmean,'r--',b5-bmean,presmean,'k-',b6-bmean,presmean,'B-.','linewidth',3);
set(gca,'ylim',[20000 max(nanmean(presShape,2))],'ydir','reverse','xlim',[-0.75 0.75],'fontsize',16,'yticklabel',[]);
xlabel('Buoyancy [K]');
title({'T Variability'});

subplot(1,4,4),
plot(b7-bmean,presmean,'r--',b8-bmean,presmean,'k-',b9-bmean,presmean,'B-.','linewidth',3);
set(gca,'ylim',[20000 max(nanmean(presShape,2))],'ydir','reverse','xlim',[-0.75 0.75],'fontsize',16,'yticklabel',[]);
xlabel('Buoyancy [K]');
title({'q Variability'});

print(gcf,'-djpeg','Figure_9.jpg');