function [sst_grad] = calc_grad(x,y,sst)
%calculate the gradient of SST's

nx=numel(x);
ny=numel(y);


a=6.371e6;%radius of earth

dlon=deg2rad(x(3:end)-x(1:end-2))*a/2;
dlon=repmat([dlon(1); dlon; dlon(end)]',ny,1);
dlat=a*deg2rad(y(3:end)-y(1:end-2)).*cosd(y(2:end-1))/2;
dlat=repmat(interp1(y(2:end-1),dlat,y,'spline','extrap'),1,nx);

dsst_dx=zeros(size(sst));
dsst_dy=zeros(size(sst));


dsst_dy(2:end-1,:) = (sst(3:end,:)-sst(1:end-2,:))./(2*dlat(2:end-1,:)); 
dsst_dy(1,:)= (sst(2,:)-sst(1,:))./(dlat(1,:));
dsst_dy(end,:)= (sst(end,:)-sst(end-1,:))./dlat(end,:);
dsst_dx(:,2:end-1) = (sst(:,3:end)-sst(:,1:end-2))./(2*dlon(:,2:end-1));
dsst_dx(:,1)= (sst(:,2)-sst(:,1))./(dlon(:,1));
dsst_dx(:,end)= (sst(:,end)-sst(:,end-1))./dlon(:,end);
sst_grad=dsst_dy+dsst_dy;


end

