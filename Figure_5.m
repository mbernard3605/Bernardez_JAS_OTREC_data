%Plots figure 5
clear variables;


%load the colormap
load('anglecolormap5.mat');

rawPath = '../Raw/';
folders = dir(rawPath);

plotDir = '../Raw/';

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

%combine the two boxes
ang_ERA_OTREC=[ang_ERA_OTREC_B1 ang_ERA_OTREC_B2];
DCIN_OTREC=[DCIN_OTREC_B1 DCIN_OTREC_B2];
II_OTREC=[II_OTREC_B1 II_OTREC_B2];
MDC_OTREC=[SFb_OTREC_B1-SFt_OTREC_B1 SFb_OTREC_B2-SFt_OTREC_B2];

%plot the instability index and MDC
figure, set(gcf,'units','normalized','outerposition',[0.05 0.05 0.6 0.75]);
scatter(II_OTREC,MDC_OTREC,5,'filled');
set(gca,'fontsize',16,'xlim',[-5 55],'ylim',[.35 .7]);
xlabel('Instbility Index [J kg^{-1} K^{-1}]','fontsize',16);
ylabel('Moisture Dipole Coeffecient','fontsize',16);
axis square;
print(gcf,'-djpeg','MDC_vs_II_OTREC_Scatter.jpg');
