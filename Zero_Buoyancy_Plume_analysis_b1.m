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
fileDir = 'Box1/';

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
mrShape = cat(2,mrShape,mr(1:78,nanindex));
presShape = cat(2,presShape,pshape(1:78,nanindex));


end

% tp=zeros(size(tempShape));hp=tp;
% for i = 1:size(tempShape,2)
%     [tp(:,i),hp(:,i),zlcl(i)]=zero_plume_buoyancy_split(tempShape(:,i),mrShape(:,i),z(1:78)*1000,presShape(:,i));
% end



save('OTREC_B1_zbp_start_ent.mat','xShape','yShape','tempShape','mrShape','presShape','z');




