%script that generates the fourth figure for the OTREC data analysis for the
%OTREC presentation in July

%Computes the O2/O1 Ratio from 3dvar data from OTREC
clear variables;
close all;

%load the era5 eofs and the colormap
load('ERA5_EOFS.mat');
lds_ERA=lds;
load('anglecolormap2.mat');
anglesindex=linspace(-180,180,64);
anglemap=[anglemap(33:end,:); anglemap(1:32,:)];

%constants
Cp=1005;
Rd=287;
g=9.81;

%rename to save headaches later
weights_old=weights;

%file directory
fileDir = 'Box2/';

%get the list of  files
files = dir([fileDir '*.nc']);

%structures to put the shaped data in
omegaShapeFull = [];
xShapeFull = [];
yShapeFull = [];
tempShape = [];
entShape = [];
satentShape = [];
mrShape = [];
satmrShape = [];
iiShape=[];
sstShape=[];
sstgradShape=[];
sfShape=[];
sfbShape=[];
sftShape=[];
uShape=[];
vShape=[];
divShape=[];
vortShape=[];
presShape = [];

%load the data
for fileIndex=1:length(files)
filename = files(fileIndex).name;
fullname = [fileDir filename];

%get the area of each flight
 lon = nc_varget(fullname,'lon');
 lat = nc_varget(fullname,'lat');
 [x{fileIndex},y{fileIndex}]=meshgrid(lat,lon);

z = nc_varget(fullname,'z');
rho = nc_varget(fullname,'rho');
pres = nc_varget(fullname,'pres')*100;
%rho = repmat(nc_varget(fullname,'rho'),1,size(pres,2),size(pres,3));
mflux = nc_varget(fullname,'mflux');

%calculate the pressure velocity
omega=-mflux.*9.81;

%Reshape all the matrices to be nPressureXntime 
pshape=reshape(pres,size(omega,1),[]);
temp = reshape(nc_varget(fullname,'temp'),size(omega,1),[])+273.15;
mr = reshape(nc_varget(fullname,'mr'),size(omega,1),[])./1000;
satmr = reshape(nc_varget(fullname,'satmr'),size(omega,1),[])./1000;
sst=nc_varget(fullname,'sst');
sst_grad=calc_grad(lon,lat,sst);
sst=reshape(sst,[],1);
sst_grad=reshape(sst_grad,[],1);
div = reshape(nc_varget(fullname,'div'),size(omega,1),[]);
vort = reshape(nc_varget(fullname,'absvort'),size(omega,1),[]);
u = reshape(nc_varget(fullname,'u'),size(omega,1),[]);
v = reshape(nc_varget(fullname,'v'),size(omega,1),[]);


[~,sfb,sft] = calc_sq(mr,satmr,nanmean(pshape,2));
sf=sfb+sft;
ent = calc_ent(temp,mr,pshape);
satent = calc_ent(temp,satmr,pshape);
ii = calc_ii(satent,z*1000);


omegaShape=reshape(omega,size(omega,1),[]);
xShape = reshape(x{fileIndex},[],1);
yShape = reshape(y{fileIndex},[],1);
nanindex=find(~isnan(mean(omegaShape,1)));
omegaShapeReal=omegaShape(:,nanindex);
xShapeReal = xShape(nanindex);
yShapeReal = yShape(nanindex);
nanIndices{fileIndex}=nanindex;
counts(fileIndex)=size(omegaShapeReal,1);

%remove the top few layers becuase they have places where pressure is
%monotonic
omegaShapeFull=cat(2,omegaShapeFull,omegaShape(1:78,nanindex));
xShapeFull=cat(2,xShapeFull,xShapeReal');
yShapeFull=cat(2,yShapeFull,yShapeReal');
tempShape = cat(2,tempShape,temp(1:78,nanindex));
entShape = cat(2,entShape,ent(1:78,nanindex));
satentShape = cat(2,satentShape,satent(1:78,nanindex));
mrShape = cat(2,mrShape,mr(1:78,nanindex));
satmrShape = cat(2,satmrShape,satmr(1:78,nanindex));
iiShape = cat(2,iiShape,ii(nanindex));
sstShape = cat(2,sstShape,sst(nanindex)');
sstgradShape = cat(2,sstgradShape,sst_grad(nanindex)');
sfShape = cat(2,sfShape,sf(nanindex));
sfbShape = cat(2,sfbShape,sfb(nanindex));
sftShape = cat(2,sftShape,sft(nanindex));
presShape = cat(2,presShape,pshape(1:78,nanindex));
uShape = cat(2,uShape,u(1:78,nanindex));
vShape = cat(2,vShape,v(1:78,nanindex));
divShape = cat(2,divShape,div(1:78,nanindex));
vortShape = cat(2,vortShape,vort(1:78,nanindex));







end

%calculate the weights and weight the vertical motion
presmean=mean(presShape,2);
dp=diff(presShape,1,2);
weightsShapeFull=[dp dp(:,end)]./sum([dp dp(:,end)]);
weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);
omega2=omegaShapeFull.*(weightsShapeFull.^0.5);

%interpolate and project the vertical motion onto the EOFs to get the PCs
height_index=find(level>=presmean(end));
lds_ERA=lds_ERA(height_index,:);
lds_ERA(:,1)=-lds_ERA(:,1);
omega2_interp=interp1(presmean,omega2,level(height_index),[],'extrap');
pcs_ERA=lds_ERA'*omega2_interp;
lds_ERA_interp=interp1(level(height_index),lds_ERA,presmean,[],'extrap');
lds_ERA_weight=lds_ERA_interp./(mean(weightsShapeFull,2).^0.5);
a=sqrt(trapz(presmean,lds_ERA_weight.^2)./(presmean(end)-presmean(1)));
lds_ERA_norm=lds_ERA_weight./a;
pcs_ERA_norm=pcs_ERA.*a';
o1_ERA=pcs_ERA_norm(1,:);
o2_ERA=pcs_ERA_norm(2,:);





%Calculate omega PCs
[lam,lds,pcs,per]=eof_routine(omega2',2);
for i = 1:size(lds,2)
dat1(:,i)=(lds(:,i))./(mean(weightsShapeFull,2).^0.5);
dat(:,i)=dat1(:,i).*sqrt((presmean(end)-presmean(1))./trapz(presmean,dat1(:,i).^2));
pc(:,i)=pcs(:,i)./sqrt((presmean(end)-presmean(1))./trapz(presmean,dat1(:,i).^2));
end
dat(:,1)=dat(:,1);
pc(:,1)=pc(:,1);




%clear omega2 weightsShapeFull










%calculate the natural and ERA5 top-heaviness angles 
ang = atan2d(pc(:,2),pc(:,1));



ang_ERA = atan2d(pcs_ERA_norm(2,:),pcs_ERA_norm(1,:));
r_ERA=sqrt(pcs_ERA_norm(2,:).^2+pcs_ERA_norm(1,:).^2);
%ang_test = atan2d(test(2,:),test(1,:));



  
evengrid=find(mod(xShapeFull,1)==0 & mod(yShapeFull,1)==0);



%sort the angles and get the index for that sort
[ang_ERA_sorted,era_index]=sort(ang_ERA);
[ang_sorted,sorted_index]=sort(ang);
%get the color values for the top-heaviness angle
ang_ERA_colors = interp1(anglesindex,anglemap,ang_ERA_sorted,[],'extrap');
ang_colors = interp1(anglesindex,anglemap,ang_sorted,[],'extrap');

%Create the smoothing window that we use
windowSize=ceil(length(ang_ERA)/90);
windowSize=45;
a=find(ang_ERA_sorted>=-45 & ang_ERA_sorted<=45);
b=gausswin(windowSize)/sum(gausswin(windowSize));
c=ones(windowSize,1);

sstSorted=sstShape(sorted_index);
sstgradSorted=sstgradShape(sorted_index);
iiSorted=iiShape(sorted_index);
sfSorted=sfShape(sorted_index);
sqSorted=sfbShape(sorted_index);



%Equation 7 from Ahmed et al 2020 
wb=.52;
wl=1-wb;
Pii=(presmean./100000).^(Rd/Cp);
b=find(presmean/100>=850);b=b(end);
c=find(presmean/100>=500);c=c(end);
entb=mean(entShape(1:b,:),1);
entl=mean(entShape(b+1:c,:),1);
entt=mean(entShape(c+1:end,:),1);
satentb=mean(satentShape(1:b,:),1);
satentl=mean(satentShape(b+1:c,:),1);
satentt=mean(satentShape(c+1:end,:),1);
Piib=mean(Pii(1:b));
Piil=mean(Pii(b+1:c));
Piit=mean(Pii(c+1:end));
plusentl=satentl-entl;

z_use=z(1:78);

Tb1=g*wb*(entb-satentl)./satentl;
Tb2=g*wl*plusentl./satentl;
Tb=Tb1-Tb2;
Tc=-0.015;
alpha=20.6;
p=alpha*(Tb-Tc).*heaviside(Tb-Tc);

DCIN=mean(satentShape(find(z_use<=2 & z_use>=1.5),:))-mean(entShape(find(z_use<=1),:));


%%

save('OTREC_B2.mat','xShape','yShape','divShape','presShape','omegaShapeFull','sstShape','o1_ERA','o2_ERA','ang_ERA','DCIN','iiShape','Tb','sfShape','sfbShape','sftShape','sstgradShape');
return;




