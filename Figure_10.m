%Script that plots figure 10

clear variables;
close all;

%load the colormap
load('anglecolormap5.mat');

rawPath = 'Raw/';
folders = dir(rawPath);

plotDir = '../Plots/ecmwf/';

load('ERA5_EOFS.mat');
eof=lds./weights.^0.5;

load('OTREC_B1.mat');
ang_ERA_OTREC_B1 = ang_ERA;
o1_ERA_OTREC_B1 = o1_ERA;
o2_ERA_OTREC_B1 = o2_ERA;
DCIN_OTREC_B1 = DCIN;
II_OTREC_B1 = iiShape;
SF_OTREC_B1 = sfShape;
SFt_OTREC_B1 = sftShape;
SFb_OTREC_B1 = sfbShape;
Tb_OTREC_B1 = Tb;
sst_OTREC_B1 = sstShape;
omega_binn_OTREC_B1 = omegaShapeFull;
pres_B1=presShape;


load('OTREC_B2.mat');
ang_ERA_OTREC_B2 = ang_ERA;
o1_ERA_OTREC_B2 = o1_ERA;
o2_ERA_OTREC_B2 = o2_ERA;
DCIN_OTREC_B2 = DCIN;
II_OTREC_B2 = iiShape;
SF_OTREC_B2 = sfShape;
SFt_OTREC_B2 = sftShape;
SFb_OTREC_B2 = sfbShape;
Tb_OTREC_B2 = Tb;
sst_OTREC_B2 = sstShape;
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

ang_ERA_OTREC=[ang_ERA_OTREC_B1 ang_ERA_OTREC_B2];
O1_ERA_OTREC=[o1_ERA_OTREC_B1 o1_ERA_OTREC_B2];
O2_ERA_OTREC=[o2_ERA_OTREC_B1 o2_ERA_OTREC_B2];
DCIN_OTREC=[DCIN_OTREC_B1 DCIN_OTREC_B2];
II_OTREC=[II_OTREC_B1 II_OTREC_B2];
SF_OTREC=[SF_OTREC_B1 SF_OTREC_B2];
MDC_OTREC=[SFb_OTREC_B1-SFt_OTREC_B1 SFb_OTREC_B2-SFt_OTREC_B2];
SST_OTREC=[sst_OTREC_B1 sst_OTREC_B2];
sst_anom2=[sst_OTREC_B1-nanmean(sst_OTREC_B1) sst_OTREC_B2-nanmean(sst_OTREC_B2)];
omega_OTREC=[omega_binn_OTREC_B1 omega_binn_OTREC_B2];
pres=[pres_B1 pres_B2]/100;
Temp_OTREC=[Temp_OTREC_B1 Temp_OTREC_B2];
MR_OTREC=[MR_OTREC_B1 MR_OTREC_B2];
z_OTREC=[z_OTREC_B1 z_OTREC_B2];

presmean=nanmean(pres,2)*100;
presmid=(presmean(2:end)+presmean(1:end-1))/2;
dp=diff(presmean);
W_EOF=interp1(level,eof(:,1:2),presmean,[],'extrap');

Temp_mean_OTREC=nanmean(Temp_OTREC,2);
Temp_detrended_OTREC=Temp_OTREC-Temp_mean_OTREC;
MR_OTREC(MR_OTREC<=1e-4)=0;
MR_mean=nanmean(MR_OTREC,2);
MR_detrended=MR_OTREC-MR_mean;
%MR_detrended(MR_detrended<=1e-4)=0;



numberData=119+175;


mdc_avg=NaN(numberData*numel(4:length(folders)-3),5);sf_avg=NaN(numberData*numel(4:length(folders)-3),5);ii_avg=NaN(numberData*numel(4:length(folders)-3),5);
qv_avg=NaN(15,numberData*numel(4:length(folders)-3),5);Temp_avg=NaN(15,numberData*numel(4:length(folders)-3),5);omega_avg=NaN(15,numberData*numel(4:length(folders)-3),5);
sst_avg=NaN(numberData*numel(4:length(folders)-3),5);ang_avg=NaN(numberData*numel(4:length(folders)-3),5);
o1_avg=NaN(numberData*numel(4:length(folders)-3),5);o2_avg=NaN(numberData*numel(4:length(folders)-3),5);
lats_avg=NaN(numberData*numel(4:length(folders)-3),5);lons_avg=NaN(numberData*numel(4:length(folders)-3),5);
counter1=0;
counter2=0;
for folderIndex = 4:length(folders)-3
savePath = [rawPath folders(folderIndex).name '/'];


load([savePath 'thermodynamics_new.mat']);

if size(mdc_all,2)==2
    
counter1=counter1+1;
indexes=[(counter1-1)*numberData+1:counter1*numberData];
for forecastIndex = 1:2

lats_avg(indexes,forecastIndex*2)=lats_all(:,forecastIndex);
lons_avg(indexes,forecastIndex*2)=lons_all(:,forecastIndex);
o1_avg(indexes,forecastIndex*2)=o1_all(:,forecastIndex);
o2_avg(indexes,forecastIndex*2)=o2_all(:,forecastIndex);
mdc_avg(indexes,forecastIndex*2)=mdc_all(:,forecastIndex);
sf_avg(indexes,forecastIndex*2)=sf_all(:,forecastIndex);
ii_avg(indexes,forecastIndex*2)=ii_all(:,forecastIndex);
qv_avg(:,indexes,forecastIndex*2)=qv_all(:,:,forecastIndex);
Temp_avg(:,indexes,forecastIndex*2)=T_all(:,:,forecastIndex);
omega_avg(:,indexes,forecastIndex*2)=omega_all(:,:,forecastIndex);
sst_avg(indexes,forecastIndex*2)=sst_all(:,forecastIndex);
ang_avg(indexes,forecastIndex*2)=ang_all(:,forecastIndex);

    
end
elseif size(mdc_all,2)==3
    
counter2=counter2+1;
indexes=[(counter2-1)*numberData+1:counter2*numberData];
   for forecastIndex = 1:3


lats_avg(indexes,forecastIndex*2-1)=lats_all(:,forecastIndex);
lons_avg(indexes,forecastIndex*2-1)=lons_all(:,forecastIndex);
o1_avg(indexes,forecastIndex*2-1)=o1_all(:,forecastIndex);
o2_avg(indexes,forecastIndex*2-1)=o2_all(:,forecastIndex);
mdc_avg(indexes,forecastIndex*2-1)=mdc_all(:,forecastIndex);
sf_avg(indexes,forecastIndex*2-1)=sf_all(:,forecastIndex);
ii_avg(indexes,forecastIndex*2-1)=ii_all(:,forecastIndex);
qv_avg(:,indexes,forecastIndex*2-1)=qv_all(:,:,forecastIndex);
Temp_avg(:,indexes,forecastIndex*2-1)=T_all(:,:,forecastIndex);
omega_avg(:,indexes,forecastIndex*2-1)=omega_all(:,:,forecastIndex);
sst_avg(indexes,forecastIndex*2-1)=sst_all(:,forecastIndex);
ang_avg(indexes,forecastIndex*2-1)=ang_all(:,forecastIndex);
       
       
     
    
    
   end
end




end


%%
%plot 1: Correlaton between QVstar and QV
presmean=nanmean([pres_B1 pres_B2],2);


SFlims_OTREC=quantile(SF_OTREC,[1/3 2/3]);
SFindex3_OTREC=SF_OTREC>SFlims_OTREC(2);
SFindex1_OTREC=SF_OTREC<SFlims_OTREC(1);
SFindex2_OTREC=SF_OTREC<=SFlims_OTREC(2) & SF_OTREC>=SFlims_OTREC(1);
MDClims_OTREC=quantile(MDC_OTREC,[1/3 2/3]);
MDCindex3_OTREC=MDC_OTREC>MDClims_OTREC(2);
MDCindex1_OTREC=MDC_OTREC<MDClims_OTREC(1);
MDCindex2_OTREC=MDC_OTREC<=MDClims_OTREC(2) & MDC_OTREC>=MDClims_OTREC(1);

for forecastIndex = 5
close all;
   ang_for = ang_avg(:,forecastIndex);
   mdc_for = mdc_avg(:,forecastIndex);
   sf_for = sf_avg(:,forecastIndex);
   [ang_sort,ang_index]=sort(ang_for);
qv_for=qv_avg(:,:,forecastIndex);
%qv_for(:,sum(qv_for,1)==0)=NaN;
T_for=Temp_avg(:,:,forecastIndex);
%T_for(:,sum(T_for,1)==0)=NaN;
omega_for=omega_avg(:,:,forecastIndex);
%omega_for(:,sum(omega_for,1)==0)=NaN;
qvstar_for=calc_qvstar(T_for,Pressure);
qvstar_OTREC=calc_qvstar(Temp_OTREC,presmean);
rh_for=qv_for./qvstar_for;



%%
t_anom=T_for-nanmean(T_for,2);
rh_anom=rh_for-nanmean(rh_for,2);

start=1;angles_plot=[];omega_plot=[];qv_plot=[];t_plot=[];rh_plot=[];
for  angleIndex = -180:1:179
    
   angles_plot(start)=nanmean(ang_for(ang_for>angleIndex & ang_for<(angleIndex+1))); 
   omega_plot(start,:)=nanmean(omega_for(:,ang_for>angleIndex & ang_for<(angleIndex+1)),2);
   rh_plot(start,:)=nanmean(rh_anom(:,ang_for>angleIndex & ang_for<(angleIndex+1)),2);
   t_plot(start,:)=nanmean(t_anom(:,ang_for>angleIndex & ang_for<(angleIndex+1)),2);
start=start+1;
end




figure('units','normalized','outerposition',[.1 .1 .75 .75]);
[c,h] = contourf(angles_plot,Pressure,rh_plot',[-1 -.15:.01:.15]);
set(h,'linecolor','none');set(gca,'ydir','reverse','ylim',[20000 100000],'xlim',[-180 180]);
colormap(redblue);colorbar;
set(gca,'fontsize',24);
title('RH anomaly composited by top-heaviness');
xticks([-180:30:180]);
caxis([-.15 .15]);
print(gcf,'-dpng',['../Plots/RH_anom_by_top-heaviness_' num2str(forecastIndex) '.png']);

figure('units','normalized','outerposition',[.1 .1 .75 .75]);
[c,h] = contourf(angles_plot,Pressure,t_plot',[-100 -1:.1:1]);
set(h,'linecolor','none');set(gca,'ydir','reverse','ylim',[20000 100000],'xlim',[-180 180]);
colormap(redblue);colorbar;
set(gca,'fontsize',24);
title('T anomaly composited by top-heaviness');
xticks([-180:30:180]);
caxis([-1 1]);
print(gcf,'-dpng',['../Plots/T_anom_by_top-heaviness_' num2str(forecastIndex) '.png']);







end
