%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Script to load NOAA OISST Hi Res Sea Surface Temperature Data,       %
%    calculate anomalies (obs - mean in reference period) and normalised  %
%    anomalies (obs - mean, scaled by standard deviation in reference     %
%    period), and create plots of these within region of interest         %
%                                                                         %
%    Key dependencies: cbrewer function for colormap                      %
%                                                                         %
%  Author - Bee Berx                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%% USER DEFINED VARIABLES 
% set region of interest, data folder location and climatology ref period
lon_extent = [-80 50];
lat_extent = [0 90];

OISST_datafolder = [getenv('isilon'),'\Data_External\NOAA_oisst.v2.highres\'];

clim_ref_period = [1991 2020];

%% set folder structure for output
if ~exist([pwd,'\OISST_AnomSST'],'dir')
    mkdir(['.\OISST_AnomSST\'])
end
if ~exist([pwd,'\OISST_NormAnomSST'],'dir')
    mkdir(['.\OISST_NormAnomSST\'])
end

%% Load OISST data
[OISST,OISST_time,OISST_lon,OISST_lat] = ...
    fun_get_OISST_timeseries(OISST_datafolder,[datenum(1981,1,1),floor(now)],lon_extent,lat_extent);
OISST_tvec = datevec(OISST_time);
OISST(OISST<-100)=NaN;

%% Calculate Climatology
[OIClimMean,OIClimSdev] = fun_get_OISST_climatology(OISST_tvec,OISST,clim_ref_period);

%% Calculate anomalies and normalised anomalies
OI_Anom = NaN.*zeros(size(OISST,1),size(OISST,2),size(OISST,3));
OI_NormAnom = NaN.*zeros(size(OISST,1),size(OISST,2),size(OISST,3));
for mm=1:12;
    idxOI =  find(OISST_tvec(:,2)==mm);
    OI_Anom(:,:,idxOI) = OISST(:,:,idxOI) - repmat(OIClimMean(:,:,mm),1,1,length(idxOI));
    OI_NormAnom(:,:,idxOI) = (OISST(:,:,idxOI) - repmat(OIClimMean(:,:,mm),1,1,length(idxOI)))./ ...
        repmat(OIClimSdev(:,:,mm),1,1,length(idxOI));
end

ButtonName = questdlg('Process all years or specific year?', ...
    'Years Question', ...
    'All', 'Range', 'Specific', 'Specific');
tim_chck=datevec(floor(now));
switch ButtonName
    case 'All'
        selectyear=1981:tim_chck(1);
    case 'Range'
        prompt={'Begin Year:','End Year:'};
        name='Input Year Range';
        numlines=1;
        defaultanswer={num2str(1981),num2str(tim_chck(1))};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        selectyear = str2num(answer{1}):str2num(answer{2});
        clear name prompt answer numlines defaultanswer
    case 'Specific'
        prompt={'Selected Year:'};
        name='Input Year';
        numlines=1;
        defaultanswer={num2str(tim_chck(1))};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        selectyear = str2num(answer{1});
        clear name prompt answer numlines defaultanswer
end % switch
clear tim_chck ButtonName

for yy=1:length(selectyear)
    OI_Anom_select = OI_Anom(:,:,OISST_tvec(:,1)==selectyear(yy));
    OI_NormAnom_select = OI_NormAnom(:,:,OISST_tvec(:,1)==selectyear(yy));


    %% mapping parameters
    %cmap = flipud(cbrewer('div','RdYlBu',22));
    cmap = cat(1,flipud(cbrewer('seq','PuBu',10)),cbrewer('seq','YlOrRd',10));
    cmap(cmap>1)=1;cmap(cmap<0)=0;

    pos12 = [0.02 0.68 0.225 0.3;0.255 0.68 0.225 0.3;0.49 0.68 0.225 0.3;0.725 0.68 0.225 0.3;...
        0.02 0.36 0.225 0.3;0.255 0.36 0.225 0.3;0.49 0.36 0.225 0.3;0.725 0.36 0.225 0.3;...
        0.02 0.04 0.225 0.3;0.255 0.04 0.225 0.3;0.49 0.04 0.225 0.3;0.725 0.04 0.225 0.3];


    %% ANOM SST OI SST relative OISST
    close all;
    figure(1)
    for mm=1:12
        [tim_chck]=intersect(OISST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(OISST_lon(:))-0.025;OISST_lon+0.025];
        lat2plot = [min(OISST_lat(:))-0.025;OISST_lat+0.025];
        heat2plot = OI_Anom_select(:,:,mm);

        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        dlon = ceil(double(max(lon2plot)-min(lon2plot))/5);
        dlat = ceil(double(max(lat2plot)-min(lat2plot))/5);

        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        ax = gca;

        caxis([-2.5 2.5])
        colormap(cmap(:,:))
        bplot_coastGLOB
        set(ax,'position',pos12(mm,:));
        if ismember(mm,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-2.5:0.5:2.5],'yticklabel',sprintf('% -3.1f\n',[-2.5:0.5:2.5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,'^o C','fontsize',10,'fontname','arial')
            set(hc,'position',[0.955 0.06 0.015 0.90])
        end
        if ~ismember(mm,[9,10,11,12])
            set(ax,'xticklabel',[])
        end
        if ~ismember(mm,[1,5,9])
            set(ax,'yticklabel',[])
        end
        text(20,20,datestr(datenum(selectyear(yy),mm,1),'mmm-yy'),'VerticalAlignment','middle','HorizontalAlignment','center')
        % title(datestr(datenum(2023,mm,1),'mmm-yy'))
        %
    end
    fun_savepngL(gcf,['.\OISST_AnomSST\SST_Maps_OISST_AnomSST_',num2str(selectyear(yy)),'_NAtl.png'])

    %% Standardised ANOM SST OI SST relative OISST
    close all;
    figure(1)
    for mm=1:12
        [tim_chck]=intersect(OISST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(OISST_lon(:))-0.025;OISST_lon+0.025];
        lat2plot = [min(OISST_lat(:))-0.025;OISST_lat+0.025];
        heat2plot = OI_NormAnom_select(:,:,mm);

        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        dlon = ceil(double(max(lon2plot)-min(lon2plot))/5);
        dlat = ceil(double(max(lat2plot)-min(lat2plot))/5);

        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        ax = gca;

        caxis([-5 5])
        colormap(cmap(:,:))
        bplot_coastGLOB
        set(ax,'position',pos12(mm,:));
        if ismember(mm,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-5:1:5],'yticklabel',sprintf('% -3.1f\n',[-5:1:5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,{'St. Dev.','Units'},'fontsize',10,'fontname','arial')
            set(hc,'position',[0.965 0.06 0.01 0.875])
        end
        if ~ismember(mm,[9,10,11,12])
            set(ax,'xticklabel',[])
        end
        if ~ismember(mm,[1,5,9])
            set(ax,'yticklabel',[])
        end
        text(20,20,datestr(datenum(selectyear(yy),mm,1),'mmm-yy'),'VerticalAlignment','middle','HorizontalAlignment','center')
        % title(datestr(datenum(2023,mm,1),'mmm-yy'))
        %
    end
    fun_savepngL(gcf,['.\OISST_NormAnomSST\SST_Maps_OISST_NormAnomSST_',num2str(selectyear(yy)),'_NAtl.png'])


    %% ANOM SST OI SST relative OISST  -SPNA only
    close all;
    figure(1)

    for mm=1:12
        [tim_chck]=intersect(OISST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(OISST_lon(:))-0.025;OISST_lon+0.025];
        lat2plot = [min(OISST_lat(:))-0.025;OISST_lat+0.025];
        heat2plot = OI_Anom_select(:,:,mm);

        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        dlon = ceil(double(max(lon2plot)-min(lon2plot))/5);
        dlat = ceil(double(max(lat2plot)-min(lat2plot))/5);

        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        ax = gca;
        caxis([-2.5 2.5])
        colormap(cmap(:,:))
        bplot_coastGLOB
        set(ax,'position',pos12(mm,:));
        if ismember(mm,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-2.5:0.5:2.5],'yticklabel',sprintf('% -3.1f\n',[-2.5:0.5:2.5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,'^o C','fontsize',10,'fontname','arial')
            set(hc,'position',[0.965 0.06 0.01 0.875])
        end
        if ~ismember(mm,[9,10,11,12])
            set(ax,'xticklabel',[])
        end
        if ~ismember(mm,[1,5,9])
            set(ax,'yticklabel',[])
        end
        set(ax,'xlim',[-50 15],'ylim',[30 65])
        text(5,32.5,datestr(datenum(selectyear(yy),mm,1),'mmm-yy'),'VerticalAlignment','middle','HorizontalAlignment','center')
    end
    fun_savepngL(gcf,['.\OISST_AnomSST\SST_Maps_OISST_AnomSST_',num2str(selectyear(yy)),'_SPNA.png'])

    %% ANOM SST OI SST relative OISST  -SPNA only
    close all;
    figure(1)

    for mm=1:12
        [tim_chck]=intersect(OISST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(OISST_lon(:))-0.025;OISST_lon+0.025];
        lat2plot = [min(OISST_lat(:))-0.025;OISST_lat+0.025];
        heat2plot = OI_NormAnom_select(:,:,mm);

        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        dlon = ceil(double(max(lon2plot)-min(lon2plot))/5);
        dlat = ceil(double(max(lat2plot)-min(lat2plot))/5);

        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        ax = gca;

        caxis([-5 5])
        colormap(cmap(:,:))
        bplot_coastGLOB
        set(ax,'position',pos12(mm,:));
        if ismember(mm,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-5:1:5],'yticklabel',sprintf('% -3.1f\n',[-5:1:5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,{'St. Dev.','Units'},'fontsize',10,'fontname','arial')
            set(hc,'position',[0.965 0.06 0.01 0.875])
        end

        if ~ismember(mm,[9,10,11,12])
            set(ax,'xticklabel',[])
        end
        if ~ismember(mm,[1,5,9])
            set(ax,'yticklabel',[])
        end
        set(ax,'xlim',[-50 15],'ylim',[30 65])
        text(5,32.5,datestr(datenum(selectyear(yy),mm,1),'mmm-yy'),'VerticalAlignment','middle','HorizontalAlignment','center')
    end
    fun_savepngL(gcf,['.\OISST_NormAnomSST\SST_Maps_OISST_NormAnomSST_',num2str(selectyear(yy)),'_SPNA.png'])
end
return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS SUPPORTING THE ABOVE CODE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SST,SST_time,sellon,sellat] = fun_get_OISST_timeseries(data_folder,time_choice,lon_choice,lat_choice)
% function fun_get_OISST_timeseries
%
% Extract time series of OI SST data from NOAA OI SST Hi Res product.
% Requires most up to date netcdf files to be saved in data_folder
% 
%
%
% USE:
%       [SST,time,lon,lat] = fun_get_OISST_timeseries(data_folder,time_choice,lon_choice,lat_choice)
%
% INPUT:
%       time_choice   = start and end time for field (matlab datenum format)
%       lon_choice    = western & eastern extent of longitude for field [-180 180]
%       lat_choice    = southern & northern extent of longitude for field [-90 90]
%
%
% OUTPUT:
%       SST (time x 1)        = monthly mean SST (deg C) from OISST
%       time                  = time stamp from datafile
%       lon                   = longitude of grid extracted (nearest long)
%       lat                   = latitude of grid extracted (nearest lat)
%
% DEPENDENCIES:
%   The function needs access to the following
%       DATA: data in NOAA format in data_folder
%
%
% EXAMPLE:
%       [SST,time,lon,lat] = fun_get_OISST_timeseries(data_folder, ...
%               [datenum(2022,1,1),floor(now)],[-25 15],[45 65]);
%

NOAA_OISSTDir = data_folder;

tim_chck = datevec(floor(now));
Years=1981:tim_chck(1);clear tim_chck
nyrs=length(Years);

%SST_time = [datenum(1981,1,1):1:datenum(2024,12,31)]';
SST_time = [sort(repmat(Years',12,1)),repmat([1:12]',nyrs,1),repmat(15,nyrs*12,1)];

idxValid = intersect(find(datenum(SST_time)>=time_choice(1)),find(datenum(SST_time)<=time_choice(2)));

SST_time = SST_time(idxValid,:);

% Loop through data year by year
for iyear=1:nyrs
    if ~ismember(Years(iyear),SST_time(:,1));continue;end
    fname = {[NOAA_OISSTDir 'sst.day.mean.' num2str(Years(iyear)) '.nc']};
    if iyear==1
        lat = ncread(fname{1},'lat');
        lon = ncread(fname{1},'lon');
        lon(lon>180)=lon(lon>180)-360;
        ilat1 = intersect(find(lat>=lat_choice(1)),find(lat<=lat_choice(2)));
        ilon2 = intersect(find(lon>=0),find(lon<=lon_choice(2)));
        ilon1 = intersect(find(lon>=lon_choice(1)),find(lon<=0));

        sellon1  = double(ncread(fname{1},'lon',[ilon1(1)],[length(ilon1)]));
        sellon2  = double(ncread(fname{1},'lon',[ilon2(1)],[length(ilon2)]));
        sellon = cat(1,sellon1,sellon2);
        sellon(sellon>180)=sellon(sellon>180)-360;

        sellat  = double(ncread(fname{1},'lat',[ilat1(1)],[length(ilat1)]));
        lsmask = cat(1,ncread([NOAA_OISSTDir,'lsmask.oisst.nc'],'lsmask',[ilon1(1),ilat1(1),1],[length(ilon1),length(ilat1),1]),...
            ncread([NOAA_OISSTDir,'lsmask.oisst.nc'],'lsmask',[ilon2(1),ilat1(1),1],[length(ilon2),length(ilat1),1]));
        lsmask(lsmask==0) = NaN;lsmask = double(lsmask);

        SST = NaN.*zeros(size(lsmask,1),size(lsmask,2),size(SST_time,1));

    end
    time_orig = datevec(datenum(1800,1,1)+ncread(fname{1},'time'));

    sst_in = cat(1,ncread(fname{1},'sst',[ilon1(1),ilat1(1),1],[length(ilon1),length(ilat1),Inf]),...
        ncread(fname{1},'sst',[ilon2(1),ilat1(1),1],[length(ilon2),length(ilat1),Inf]));
    sst_in = double(sst_in);
    sst_in(sst_in==-9.969209968386869e+36)=NaN;

    time_mon = unique(time_orig(:,[1:2]),"rows");

    for nn=1:size(time_mon,1)
        idx = intersect(find(time_orig(:,1)==time_mon(nn,1)),...
            find(time_orig(:,2)==time_mon(nn,2)));
        if length(idx)~=eomday(time_mon(nn,1),time_mon(nn,2));continue;end
        [~,~,idx2] = intersect(time_mon(nn,1:2),SST_time(:,1:2),'rows');
        SST(:,:,idx2) = mean(sst_in(:,:,idx).*repmat(lsmask,1,1,length(idx)),3);
    end
end

SST_time=datenum(SST_time);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [OIClimMean,OIClimSdev] = fun_get_OISST_climatology(OISST_tvec,OISST,clim_ref_period)
% sst climatology - mean and st dev in ref period

OIClimMean = NaN.*zeros(size(OISST,1),size(OISST,2),12);
OIClimSdev = NaN.*zeros(size(OISST,1),size(OISST,2),12);
for mm=1:12
    ClimIdx  = intersect(intersect(find(OISST_tvec(:,1)>=clim_ref_period(1)),...
        find(OISST_tvec(:,1)<=clim_ref_period(2))),...
        find(OISST_tvec(:,2)==mm));
    OIClimMean(:,:,mm) = mean(OISST(:,:,ClimIdx),3);
    OIClimSdev(:,:,mm) = std(OISST(:,:,ClimIdx),[],3);
    clear ClimIdx
end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_savepngL(figh,filename)
% function fun_savepngL
% save figure in landscape


set(figh,'paperorientation','landscape','papertype','a4','paperpositionmode','auto',...
    'paperunits','centimeters','paperposition',[0.6 0.6 28.4 19.7])
print(figh, '-dpng', '-r300', filename)
end