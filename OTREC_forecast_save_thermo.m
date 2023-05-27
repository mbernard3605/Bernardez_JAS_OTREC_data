clear variables;
close all;


addpath('~/Documents/mexnc/');
addpath('~/Documents/snctools/');
%close all;

%load colormap and eofs
load('anglecolormap2.mat');
load('ERA5_EOFS.mat');
lds_ERA=lds;
%anglemap=flip(anglemap);
anglemap=[anglemap(33:end,:); anglemap(1:32,:)];
colormap(anglemap);
weights_old=weights;


rawPath = 'Raw/';
plotDir = 'Raw/';
plotDir1 = 'Raw/';
folders = dir(rawPath);
load coastlines
		     %folderIndex=3;
             load('anglecolormap5.mat');
             
             
box2lon=[-89:0.5:-86 ones(size(3:0.5:11))*-86 -86:-0.5:-89 ones(size(3:0.5:11))*-89];
box2lat=[ones(size(-94:0.5:-91))*3 3:0.5:11 ones(size(-94:0.5:-91))*11 11:-0.5:3];
box1alon=[-81:0.5:-78 ones(size(4:0.5:7))*-78 -78:-0.5:-81 ones(size(4:0.5:7))*-81];
box1alat=[ones(size(-81:0.5:-78))*4 4:0.5:7 ones(size(-78:-0.5:-81))*7 4:0.5:7];
box1blon=[-79.5:0.5:-81 -81:0.5:-83.5 -83.5:0.5:-82 -82:0.5:-79.5];
box1blat=[11:0.5:13 13:0.5:11.5 11.5:0.5:9.5 9.5:0.5:11];
             
             
 started=0;           
 runningCount=0;
for folderIndex = 3:length(folders)-3
    runningCount=runningCount+1
            folderName = folders(folderIndex).name;
		    files = dir([rawPath folderName '/*.nc']);
    if isempty(files)
        continue;
    else
        
    
        
        
        modelRun = files(1).name(end-11:end-4);
        savePath = [plotDir folders(folderIndex).name '/'];
if ~exist(savePath, 'dir')
mkdir(savePath)
end
        
        %if ~exist([savePath 'thermodynamics.mat'])
        modelhour= str2num(modelRun(5:6));
        fileHour = str2num(folderName(5:6));
        
        (datenum([2019 str2num(modelRun(1:2)) str2num(modelRun(3:4)) str2num(modelRun(5:6)) 0 0])-...
            datenum([2019 str2num(folderName(1:2)) str2num(folderName(3:4)) str2num(folderName(5:6)) 0 0]))*24;
        
        
        if modelhour >= 10
            Prefix = ['model.ECMWF_diurnal_average.2019' folderName '00.00' num2str(modelhour) '_'];
        else
            Prefix = ['model.ECMWF_diurnal_average.2019' folderName '00.0' num2str(modelhour) '_'];
        end
        
        if fileHour > 0
            filePath1 = [rawPath folders(folderIndex).name '/' files(1).name];
            filePath2 = [rawPath folders(folderIndex).name '/' files(5).name];
            filePath3 = [rawPath folders(folderIndex).name '/' files(9).name];
            
            
            
            
            Pressure = nc_varget(filePath1,'lv_ISBL2')*100;
            
            
            % if started==0
            % Pressure = nc_varget(filePath,'lv_ISBL2')*100;
            % PressureInd=find(Pressure>=10000);
            % Pressure=[0; Pressure(PressureInd)];
            %
            % started=1;
            % end
            
            
            W{1} = nc_varget(filePath1,'W_GDS0_ISBL');
            W{2} = nc_varget(filePath2,'W_GDS0_ISBL');
            W{3} = nc_varget(filePath3,'W_GDS0_ISBL');
            T{1} = nc_varget(filePath1,'T_GDS0_ISBL');
            T{2} = nc_varget(filePath2,'T_GDS0_ISBL');
            T{3} = nc_varget(filePath3,'T_GDS0_ISBL');
            qv{1} = nc_varget(filePath1,'Q_GDS0_ISBL');
            qv{2} = nc_varget(filePath2,'Q_GDS0_ISBL');
            qv{3} = nc_varget(filePath3,'Q_GDS0_ISBL');
            z{1} = nc_varget(filePath1,'GH_GDS0_ISBL');
            z{2} = nc_varget(filePath2,'GH_GDS0_ISBL');
            z{3} = nc_varget(filePath3,'GH_GDS0_ISBL');
            
            
            sst{1}=nc_varget(filePath1,'SSTK_GDS0_SFC');
            landMask=(sst{1}>=280);
            sst{2} = nc_varget(filePath2,'SSTK_GDS0_SFC');
            sst{3} = nc_varget(filePath3,'SSTK_GDS0_SFC');
            
            
            
            %end
            lon=nc_varget(filePath1,'g0_lon_1');
            lat=nc_varget(filePath1,'g0_lat_0');
            [x,y]=meshgrid(lon,lat);
                            box1alons=ismember(lon,box1alon);
                box1blons=ismember(lon,box1blon);
                box2lats=ismember(lat,box2lat);
                box1alats=ismember(lat,box1alat);
                box1blats=ismember(lat,box1blat);
                box2lons=ismember(lon,box2lon);
                
                a=box1alons | box1blons;
                b=box1alats | box1blats;
                [xa,xb]=meshgrid(a,b);
                [ya,yb]=meshgrid(box2lons,box2lats);
                box1=(xa & xb & landmask);%box1(box1==0)=NaN;
                box2=(ya & yb & landmask);%box2(box2==0)=NaN;
                count=sum(sum(box1+box2));
            
            
            lons_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            lats_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            mdc_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            sf_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            qv_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            T_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            z_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            omega_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            ii_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            sst_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            ang_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            o1_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            o2_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            
            for forecastIndex = 1:3
                
                
                Wtemp=W{forecastIndex}.*permute(repmat(landMask,1,1,15),[3 1 2]);
                Ttemp=T{forecastIndex}.*permute(repmat(landMask,1,1,15),[3 1 2]);
                qvtemp=qv{forecastIndex}.*permute(repmat(landMask,1,1,15),[3 1 2]);
                ztemp=z{forecastIndex}.*permute(repmat(landMask,1,1,15),[3 1 2]);
                ssttemp=sst{forecastIndex}.*landMask;
                omega = reshape(Wtemp,size(Wtemp,2)*size(Wtemp,3)*size(Wtemp,4),size(Wtemp,1));
                
                dp = diff(Pressure,1,1);
                
                %w = dp./sum(dp);
                %omega=detrend(DataShaped,'constant');
                
                weightsShapeFull=[dp; dp(end)]'./sum([dp; dp(end)]');
                weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);
                omega2=omega.*(weightsShapeFull.^0.5);
                
                %interpolate and project the vertical motion onto the EOFs to get the PCs
                height_index=find(level>=Pressure(end));
                lds_ERA(:,1)=-lds_ERA(:,1);
                omega2_interp=interp1(Pressure,omega2',level,[],'extrap');
                pcs_ERA=lds_ERA'*omega2_interp;
                lds_ERA_interp=interp1(level,lds_ERA,Pressure,[],'extrap');
                lds_ERA_weight=lds_ERA_interp./(nanmean(weightsShapeFull,1)'.^0.5);
                a=sqrt(trapz(Pressure,lds_ERA_weight.^2)./(Pressure(end)-Pressure(1)));
                lds_ERA_norm=lds_ERA_weight./a;
                pcs_ERA_norm=pcs_ERA.*a';
                
                o1=pcs_ERA_norm(1,:);
                o2=pcs_ERA_norm(2,:);
                
                o1 = reshape(o1,size(W{1},2),size(W{1},3));                
                o2 = reshape(o2,size(W{1},2),size(W{1},3));
               
                
                ang = atan2d(o2,o1);
                ang = reshape(ang,size(W{1},2),size(W{1},3));
                
                presmid=(Pressure(2:end)+Pressure(1:end-1))/2;
                dp=diff(Pressure);
                
                
                
                o1filt=imgaussfilt(o1);
                o2filt=imgaussfilt(o2);
                
                angfilt = atan2d(o2filt,o1filt);
                tp=zeros(size(Ttemp));hp=tp;senv=tp;henv=tp;henv_sat=tp;
                qvtemp=qv{forecastIndex};
                satqv=zeros(size(qvtemp));
                ent = zeros(size(qvtemp));
                satent=ent; sf=zeros(size(qvtemp,2),size(qvtemp,3));
                sfb=sf;sft=sf;ii=sf;
                for iIndex = 1:size(Ttemp,2)
                    for jIndex = 1:size(Ttemp,3)
                        satqv(:,iIndex,jIndex)=calc_qvstar(Ttemp(:,iIndex,jIndex),Pressure);
                        [~,sfb(iIndex,jIndex),sft(iIndex,jIndex)] = calc_sq(flip(qvtemp(:,iIndex,jIndex),1),flip(satqv(:,iIndex,jIndex),1),flip(Pressure,1));
                        ent(:,iIndex,jIndex) = calc_ent(flip(Ttemp(:,iIndex,jIndex),1),flip(qvtemp(:,iIndex,jIndex),1),flip(Pressure,1));
                        satent(:,iIndex,jIndex) = calc_ent(flip(Ttemp(:,iIndex,jIndex),1),flip(satqv(:,iIndex,jIndex),1),flip(Pressure,1));
                        ii(iIndex,jIndex) = calc_ii(satent(:,iIndex,jIndex),ztemp(:,iIndex,jIndex));
                        
                        
                        
                        
                        
                    end
                end
                
                
                
                sf=sfb+sft;
                mdc=sfb-sft;
                
                
                
                lons_all(:,forecastIndex)=cat(1,reshape(x(box1),[],1),reshape(x(box2),[],1));
                lats_all(:,forecastIndex)=cat(1,reshape(y(box1),[],1),reshape(y(box2),[],1));
                mdc_all(:,forecastIndex)=cat(1,reshape(mdc(box1),[],1),reshape(mdc(box2),[],1));
                ii_all(:,forecastIndex)=cat(1,reshape(ii(box1),[],1),reshape(ii(box2),[],1));
                sf_all(:,forecastIndex)=cat(1,reshape(sf(box1),[],1),reshape(sf(box2),[],1));
                qv_all(:,:,forecastIndex)=cat(2,reshape(qvtemp(:,box1),size(qvtemp,1),[],1),reshape(qvtemp(:,box2),size(qvtemp,1),[],1));
                T_all(:,:,forecastIndex)=cat(2,reshape(Ttemp(:,box1),size(Ttemp,1),[],1),reshape(Ttemp(:,box2),size(qvtemp,1),[],1));
                z_all(:,:,forecastIndex)=cat(2,reshape(ztemp(:,box1),size(ztemp,1),[],1),reshape(ztemp(:,box2),size(qvtemp,1),[],1));
                omega_all(:,:,forecastIndex)=cat(2,reshape(Wtemp(:,box1),size(Wtemp,1),[],1),reshape(Wtemp(:,box2),size(Wtemp,1),[],1));
                sst_all(:,forecastIndex)=cat(1,reshape(ssttemp(box1),[],1),reshape(ssttemp(box2),[],1));
                ang_all(:,forecastIndex)=cat(1,reshape(ang(box1),[],1),reshape(ang(box2),[],1));
                o1_all(:,forecastIndex)=cat(1,reshape(o1(box1),[],1),reshape(o1(box2),[],1));
                o2_all(:,forecastIndex)=cat(1,reshape(o2(box1),[],1),reshape(o2(box2),[],1));
                
                
                                
                

                
            end
           
            save([savePath 'thermodynamics_new.mat'],'Pressure','lons_all','lats_all','mdc_all','sf_all','ii_all','qv_all','T_all','omega_all','ang_all','sst_all','o1_all','o2_all');
            
        else
            filePath1 = [rawPath folders(folderIndex).name '/' files(3).name];
            filePath2 = [rawPath folders(folderIndex).name '/' files(7).name];
            
            
            Pressure = nc_varget(filePath1,'lv_ISBL2')*100;
            
            W{1} = nc_varget(filePath1,'W_GDS0_ISBL');
            W{2} = nc_varget(filePath2,'W_GDS0_ISBL');
            T{1} = nc_varget(filePath1,'T_GDS0_ISBL');
            T{2} = nc_varget(filePath2,'T_GDS0_ISBL');
            qv{1} = nc_varget(filePath1,'Q_GDS0_ISBL');
            qv{2} = nc_varget(filePath2,'Q_GDS0_ISBL');
            z{1} = nc_varget(filePath1,'GH_GDS0_ISBL');
            z{2} = nc_varget(filePath2,'GH_GDS0_ISBL');
            
            
            sst{1}=nc_varget(filePath1,'SSTK_GDS0_SFC');
            landMask=(sst{1}>=280);
            %landMask(~landMask)=NaN;
            sst{2} = nc_varget(filePath2,'SSTK_GDS0_SFC');
            
            
            
            %end
            lon=nc_varget(filePath1,'g0_lon_1');
            lat=nc_varget(filePath1,'g0_lat_0');
            [x,y]=meshgrid(lon,lat);
                            box1alons=ismember(lon,box1alon);
                box1blons=ismember(lon,box1blon);
                box2lats=ismember(lat,box2lat);
                box1alats=ismember(lat,box1alat);
                box1blats=ismember(lat,box1blat);
                box2lons=ismember(lon,box2lon);
                
                a=box1alons | box1blons;
                b=box1alats | box1blats;
                [xa,xb]=meshgrid(a,b);
                [ya,yb]=meshgrid(box2lons,box2lats);
                box1=(xa & xb & landMask);%box1(box1==0)=NaN;
                box2=(ya & yb & landMask);%box2(box2==0)=NaN;
                count=sum(sum(box1+box2));
            
            
            lons_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            lats_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            mdc_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            sf_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            qv_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            T_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            z_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            omega_all=zeros(size(W{1},1),(sum(sum(box1))+sum(sum(box2))),2);
            ii_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            sst_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            ang_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            o1_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            o2_all=zeros((sum(sum(box1))+sum(sum(box2))),2);
            
            for forecastIndex = 1:2
                Wtemp=W{forecastIndex};
                Ttemp=T{forecastIndex};
                qvtemp=qv{forecastIndex};
                ztemp=z{forecastIndex};
                ssttemp=sst{forecastIndex};
                
                lons_all(:,forecastIndex)=cat(1,reshape(x(box1),[],1),reshape(x(box2),[],1));
                lats_all(:,forecastIndex)=cat(1,reshape(y(box1),[],1),reshape(y(box2),[],1));
               
                qv_all(:,:,forecastIndex)=cat(2,reshape(qvtemp(:,box1),size(qvtemp,1),[],1),reshape(qvtemp(:,box2),size(qvtemp,1),[],1));
                T_all(:,:,forecastIndex)=cat(2,reshape(Ttemp(:,box1),size(Ttemp,1),[],1),reshape(Ttemp(:,box2),size(Ttemp,1),[],1));
                z_all(:,:,forecastIndex)=cat(2,reshape(ztemp(:,box1),size(ztemp,1),[],1),reshape(ztemp(:,box2),size(ztemp,1),[],1));
                omega_all(:,:,forecastIndex)=cat(2,reshape(Wtemp(:,box1),size(Wtemp,1),[],1),reshape(Wtemp(:,box2),size(Wtemp,1),[],1));
                sst_all(:,forecastIndex)=cat(1,reshape(ssttemp(box1),[],1),reshape(ssttemp(box2),[],1));

                
                omega = reshape(Wtemp,size(Wtemp,2)*size(Wtemp,3),size(Wtemp,1));
                
                dp = diff(Pressure,1,1);
                
                %w = dp./sum(dp);
                %omega=detrend(DataShaped,'constant');
                
                weightsShapeFull=[dp; dp(end)]'./sum([dp; dp(end)]');
                weightsShapeFull(weightsShapeFull<=0)=-weightsShapeFull(weightsShapeFull<=0);
                omega2=omega.*(weightsShapeFull.^0.5);
                
                %interpolate and project the vertical motion onto the EOFs to get the PCs
                height_index=find(level>=Pressure(end));
                lds_ERA(:,1)=-lds_ERA(:,1);
                omega2_interp=interp1(Pressure,omega2',level,[],'extrap');
                pcs_ERA=lds_ERA'*omega2_interp;
                lds_ERA_interp=interp1(level,lds_ERA,Pressure,[],'extrap');
                lds_ERA_weight=lds_ERA_interp./(nanmean(weightsShapeFull,1)'.^0.5);
                a=sqrt(trapz(Pressure,lds_ERA_weight.^2)./(Pressure(end)-Pressure(1)));
                lds_ERA_norm=lds_ERA_weight./a;
                pcs_ERA_norm=pcs_ERA.*a';
                
                o1=pcs_ERA_norm(1,:);
                o2=pcs_ERA_norm(2,:);
                

                
                o1 = reshape(o1,size(W{1},2),size(W{1},3));                
                o2 = reshape(o2,size(W{1},2),size(W{1},3));
                
                ang = atan2d(o2,o1);
            
                ang = reshape(ang,size(W{1},2),size(W{1},3));
                
                
                o1filt=imgaussfilt(o1);
                o2filt=imgaussfilt(o2);
                
                angfilt = atan2d(o2filt,o1filt);
                
                tp=zeros(size(Ttemp));hp=tp;senv=tp;henv=tp;henv_sat=tp;
                satqv=zeros(size(qvtemp));
                ent = zeros(size(qvtemp));
                satent=ent; sf=zeros(size(qvtemp,2),size(qvtemp,3));
                sfb=sf;sft=sf;ii=sf;
                for iIndex = 1:size(Ttemp,2)
                    for jIndex = 1:size(Ttemp,3)
                        satqv(:,iIndex,jIndex)=calc_qvstar(Ttemp(:,iIndex,jIndex),Pressure);
                        [~,sfb(iIndex,jIndex),sft(iIndex,jIndex)] = calc_sq(flip(qvtemp(:,iIndex,jIndex),1),flip(satqv(:,iIndex,jIndex),1),flip(Pressure,1));
                        ent(:,iIndex,jIndex) = calc_ent(flip(Ttemp(:,iIndex,jIndex),1),flip(qvtemp(:,iIndex,jIndex),1),flip(Pressure,1));
                        satent(:,iIndex,jIndex) = calc_ent(flip(Ttemp(:,iIndex,jIndex),1),flip(satqv(:,iIndex,jIndex),1),flip(Pressure,1));
                        ii(iIndex,jIndex) = calc_ii(satent(:,iIndex,jIndex),ztemp(:,iIndex,jIndex));
                        
                        
                        
                        
                        
                    end
                end
                
                
                
                sf=sfb+sft;
                mdc=sfb-sft;
                
                
                
                
                                
                

                
            end
           
            save([savePath 'thermodynamics_new.mat'],'Pressure','lons_all','lats_all','mdc_all','sf_all','ii_all','qv_all','T_all','omega_all','ang_all','sst_all','o1_all','o2_all');

            
        end
        

    end
   
end






