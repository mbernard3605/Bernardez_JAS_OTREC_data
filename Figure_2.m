%Script that plots part of figure 2


%Computes the O2/O1 Ratio from 3dvar data from OTREC
clear variables;
close all;

%load the era5 eofs and the colormap
load('ERA5_EOFS.mat');
lds_ERA=lds;
load('anglecolormap2.mat');
anglesindex=linspace(-180,180,64);
anglemap=[anglemap(33:end,:); anglemap(1:32,:)];

%rename to save headaches later
weights_old=weights;

%file directory
fileDir = '../merged_files/';

%get the list of  files
files = dir([fileDir '*.nc']);

%structures to put the shaped data in
omegaShapeFull = [];
xShapeFull = [];
yShapeFull = [];
presShape = [];


%load the data
for fileIndex=1:length(files)
filename = files(fileIndex).name;
fullname = [fileDir filename];

%get the area of each flight
 lon = nc_varget(fullname,'lon');
 lat = nc_varget(fullname,'lat');
 [x{fileIndex},ang_ERA_sorted{fileIndex}]=meshgrid(lat,lon);

z = nc_varget(fullname,'z');
pres = nc_varget(fullname,'pres')*100;
rho = nc_varget(fullname,'rhotherm');
w = nc_varget(fullname,'w');

%calculate the pressure velocity
omega=-w.*rho*9.81;

%Reshape all the matrices to be nPressureXntime 
pshape=reshape(pres,[],size(w,3));



w_col = squeeze(nanmean(nanmean(w,1),2));

omegaShape=reshape(omega,[],size(omega,3));
xShape = reshape(x{fileIndex},[],1);
yShape = reshape(ang_ERA_sorted{fileIndex},[],1);
nanindex=find(~isnan(mean(omegaShape,2)));
omegaShapeReal=omegaShape(nanindex,:);
xShapeReal = xShape(nanindex);
yShapeReal = xShape(nanindex);
nanIndices{fileIndex}=nanindex;
counts(fileIndex)=size(omegaShapeReal,1);

%remove the top few layers becuase they have places where pressure is
%monotonic
omegaShapeFull=cat(1,omegaShapeFull,omegaShape(nanindex,1:78));
xShapeFull=cat(1,xShapeFull,xShapeReal);
yShapeFull=cat(1,yShapeFull,yShapeReal);
presShape = cat(1,presShape,pshape(nanindex,1:78));



end

%calculate the weights and weight the vertical motion
presmean=mean(presShape);
dp=diff(presShape,1,2);
weightsShapeFull=[dp dp(:,end)]./sum([dp dp(:,end)]);
weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);
omega2=omegaShapeFull.*(weightsShapeFull.^0.5);

%interpolate and project the vertical motion onto the EOFs to get the PCs
height_index=find(level>=presmean(end));
lds_ERA=lds_ERA(height_index,:);
lds_ERA(:,1)=-lds_ERA(:,1);
omega2_interp=interp1(presmean,omega2',level(height_index),[],'extrap');
pcs_ERA=lds_ERA'*omega2_interp;
lds_ERA_interp=interp1(level(height_index),lds_ERA,presmean,[],'extrap');
lds_ERA_weight=lds_ERA_interp./(mean(weightsShapeFull,1)'.^0.5);
a=sqrt(trapz(presmean,lds_ERA_weight.^2)./(presmean(end)-presmean(1)));
lds_ERA_norm=lds_ERA_weight./a;
pcs_ERA_norm=pcs_ERA.*a';





%calculate top-heaviness angles

ang_ERA = atan2d(pcs_ERA_norm(2,:),pcs_ERA_norm(1,:));

%ang_test = atan2d(test(2,:),test(1,:));



  
evengrid=find(mod(xShapeFull,1)==0 & mod(yShapeFull,1)==0);



%sort the angles and get the index for that sort
[ang_ERA_sorted,era_index]=sort(ang_ERA);
%get the color values for the top-heaviness angle
ang_ERA_colors = interp1(anglesindex,anglemap,ang_ERA_sorted,[],'extrap');

%Create the smoothing window that we use
windowSize=45;
a=find(ang_ERA_sorted>=-45 & ang_ERA_sorted<=45);
b=gausswin(windowSize)/sum(gausswin(windowSize));

%%
%Composite the vertical motion by the top-heaivness angle
start=1;angles_plot=[];omega_plot=[];
for  angleIndex = -180:1:179
    
   angles_plot(start)=nanmean(ang_ERA(ang_ERA>angleIndex & ang_ERA<(angleIndex+1))); 
   omega_plot(start,:)=nanmean(omegaShapeFull(ang_ERA>angleIndex & ang_ERA<(angleIndex+1),:));

start=start+1;
end

close all;

figure('units','normalized','outerposition',[.1 .1 .75 .75]);
[c,h] = contourf(angles_plot,presmean(1:72),omega_plot(:,1:72)',-.6:.1:.6);
set(h,'linecolor','none');set(gca,'ydir','reverse');
colormap(redblue);colorbar;
set(gca,'fontsize',24);
title('Vertical Motion composited by top-heaviness');
xticks([-180:30:180]);
caxis([-.6 .6])
print(gcf,'-dpng','../omega_by_top-heaviness.png');

