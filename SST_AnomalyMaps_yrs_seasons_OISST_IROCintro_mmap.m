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
addpath('.\cbrewer')

%% USER DEFINED VARIABLES
% set region of interest, data folder location and climatology ref period
%lon_extent = [-20,15];
%lat_extent = [45,65];
lon_extent = [-80 30];
lat_extent = [20 85];

OISST_datafolder = ['I:\Data_External\NOAA_oisst.v2.highres\'];

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

%% Select year or year range to plot
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

%% Work through selected years and make plots
for yy=1:length(selectyear)
    yidx = find(OISST_tvec(:,1)==selectyear(yy))-1;

    if length(yidx)~=12;disp(['insufficient months in year ' sprintf('%4d',selectyear(yy))]);continue;end

    OI_Anom_select = OI_Anom(:,:,yidx);
    OI_NormAnom_select = OI_NormAnom(:,:,yidx);

    OI_ann_Anom_select = mean(OI_Anom_select,3);
    OI_ann_NormAnom_select = mean(OI_NormAnom_select,3);


    OI_sea_Anom_select = NaN.*OI_Anom_select(:,:,1:4);
    OI_sea_NormAnom_select = NaN.*OI_NormAnom_select(:,:,1:4);
    sea_idx = [1,2,3;4,5,6;7,8,9;10,11,12];
    for ss=1:4;
        OI_sea_Anom_select(:,:,ss) = mean(OI_Anom_select(:,:,sea_idx(ss,:)),3);
        OI_sea_NormAnom_select(:,:,ss) = mean(OI_NormAnom_select(:,:,sea_idx(ss,:)),3);
    end

    %% mapping parameters
    %cmap = flipud(cbrewer('div','RdYlBu',22));
    cmap = cat(1,flipud(cbrewer('seq','PuBu',10)),cbrewer('seq','YlOrRd',10));
    cmap(cmap>1)=1;cmap(cmap<0)=0;

    tmp= cbrewer('div','PiYG',3);
    amap = cat(1,flipud(cbrewer('seq','Blues',6)),tmp(2,:),tmp(2,:),cbrewer('seq','Reds',6));
    amap(amap>1)=1;amap(amap<0)=0;
    clear tmp

    pl = 0.45;
    ph = 0.3;
    pos5_1_4 = [0.02 0.04+2*(ph+0.005) 2*pl+0.005 ph;
        0.02 0.04+1*(ph+0.005) pl ph;
        0.02+1*(pl+0.005) 0.04+1*(ph+0.005) pl ph;
        0.02 0.04 pl ph;
        0.02+1*(pl+0.005) 0.04 pl ph];


    %% ANOM SST OI SST relative OISST
    OI_toPlot_Anom = cat(3,OI_ann_Anom_select,OI_sea_Anom_select);
    OI_toPlot_NormAnom = cat(3,OI_ann_NormAnom_select,OI_sea_NormAnom_select);
    OI_toPlot_Titles = {'Annual','DJF','MAM','JJA','SON'};

    lon2plot = [min(OISST_lon(:))-0.025;OISST_lon+0.025];
    lat2plot = [min(OISST_lat(:))-0.025;OISST_lat+0.025];
    dlon = ceil(double(max(lon2plot)-min(lon2plot))/5);
    dlat = ceil(double(max(lat2plot)-min(lat2plot))/5);

    %% ANOM SST OI SST relative OISST
    close all;
    for ss=1:5;
        ax = axes('position',pos5_1_4(ss,:));hold on
        heat2plot = OI_toPlot_Anom(:,:,ss);
        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        caxis([-2.5 2.5])
        colormap(cmap(:,:))
        bplot_coastGLOB
        set(ax,'position',pos5_1_4(ss,:),'xlim',lon_extent,'ylim',lat_extent,'Box','on');

        if ismember(ss,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-2.5:0.5:2.5],'yticklabel',sprintf('% -3.1f\n',[-2.5:0.5:2.5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,'^o C','fontsize',10,'fontname','arial')
            set(hc,'position',[0.955 0.04 0.015 3*ph-0.02])
        end
        if ismember(ss,[2,3])
            set(ax,'xticklabel',[])
        end
        if ismember(ss,[3,5])
            set(ax,'yticklabel',[])
        end
        if ismember(ss,[1])
            set(ax,'xaxislocation','top')
        end
        if ss==1
            text(lon_extent(1)+abs(diff(lon_extent))*0.925,lat_extent(1)+abs(diff(lat_extent))*0.01,{num2str(selectyear(yy));OI_toPlot_Titles{ss}},'VerticalAlignment','bottom','HorizontalAlignment','center')
        else
            text(lon_extent(1)+abs(diff(lon_extent))*0.925,lat_extent(1)+abs(diff(lat_extent))*0.01,OI_toPlot_Titles{ss},'VerticalAlignment','bottom','HorizontalAlignment','center')
        end
    end
    fun_savepngL(gcf,['.\OISST_AnomSST\SST_Maps_OISST_AnomSST_',num2str(selectyear(yy)),'_YearSeasons_IROC.png'])
    %print(gcf, '-dpng', '-r300',['.\OISST_AnomSST\SST_Maps_OISST_AnomSST_',num2str(selectyear(yy)),'_YearSeasons_IROC.png'])

    %% Standardised ANOM SST OI SST relative OISST
    close all;
    for ss=1:5;
        ax = axes('position',pos5_1_4(ss,:));hold on
        heat2plot = OI_toPlot_NormAnom(:,:,ss);
        heat2plot = cat(1,heat2plot,NaN.*heat2plot(1,:));
        heat2plot = cat(2,heat2plot,heat2plot(:,1).*NaN);
        pcolor(lon2plot,lat2plot,heat2plot');shading flat;
        caxis([-3.5 3.5])
        colormap(amap(:,:))

        bplot_coastGLOB
        set(ax,'position',pos5_1_4(ss,:),'xlim',lon_extent,'ylim',lat_extent,'Box','on');

        if ismember(ss,[1])
            [hc]=colorbar(ax,'eastoutside');%([.05 .9],.05,CS,CH,'endpiece','no','axfrac',.025,'levels','set','fontsize',10,'fontname','arial');
            set(hc,'ytick',[-3.5:0.5:3.5],'yticklabel',sprintf('% -3.1f\n',[-3.5:0.5:3.5]'))%,'yticklabelrotation',0,'ticklength',[0.01     0.05])
            title(hc,{'St. Dev.','Units'},'fontsize',10,'fontname','arial')
            set(hc,'position',[0.955 0.04 0.015 3*ph-0.02])
        end
        if ismember(ss,[2,3])
            set(ax,'xticklabel',[])
        end
        if ismember(ss,[3,5])
            set(ax,'yticklabel',[])
        end
        if ismember(ss,[1])
            set(ax,'xaxislocation','top')
        end
        if ss==1
            text(lon_extent(1)+abs(diff(lon_extent))*0.925,lat_extent(1)+abs(diff(lat_extent))*0.01,{num2str(selectyear(yy));OI_toPlot_Titles{ss}},'VerticalAlignment','bottom','HorizontalAlignment','center')
        else
            text(lon_extent(1)+abs(diff(lon_extent))*0.925,lat_extent(1)+abs(diff(lat_extent))*0.01,OI_toPlot_Titles{ss},'VerticalAlignment','bottom','HorizontalAlignment','center')
        end
    end
    fun_savepngL(gcf,['.\OISST_NormAnomSST\SST_Maps_OISST_NormAnomSST_',num2str(selectyear(yy)),'_YearSeasons_IROC.png'])
    %print(gcf, '-dpng', '-r300',['.\OISST_NormAnomSST\SST_Maps_OISST_NormAnomSST_',num2str(selectyear(yy)),'_YearSeasons_IROC.png'])

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


fname = {[NOAA_OISSTDir 'sst.mon.mean.nc']};
iname = {[NOAA_OISSTDir 'icec.mon.mean.nc']};

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

time_orig = datevec(datenum(1800,1,1)+ncread(fname{1},'time'));

sst_in = cat(1,ncread(fname{1},'sst',[ilon1(1),ilat1(1),1],[length(ilon1),length(ilat1),Inf]),...
    ncread(fname{1},'sst',[ilon2(1),ilat1(1),1],[length(ilon2),length(ilat1),Inf]));
sst_in = double(sst_in);
sst_in(sst_in==-9.969209968386869e+36)=NaN;

ice_in = cat(1,ncread(iname{1},'icec',[ilon1(1),ilat1(1),1],[length(ilon1),length(ilat1),Inf]),...
    ncread(iname{1},'icec',[ilon2(1),ilat1(1),1],[length(ilon2),length(ilat1),Inf]));
ice_in = double(ice_in);
ice_in(ice_in==-9.969209968386869e+36)=NaN;

icemask = abs(isnan(ice_in));
icemask(icemask==0)=NaN;

clear ice_in;

sst_in_masked = sst_in.* icemask .*repmat(lsmask,1,1,size(sst_in,3));

idxValid = intersect(find(datenum(time_orig)>=time_choice(1)),find(datenum(time_orig)<=time_choice(2)));

SST_time = datenum(time_orig(idxValid,:));
SST = sst_in_masked(:,:,idxValid);

clear ice_in icemask time_mon sst_in sst_in_masked Tout Tin


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
    mask = sum(~isnan(OISST(:,:,ClimIdx)),3);
    mask(mask<25)=NaN;mask(~isnan(mask))=1;%masks out where there are less than 25 of 30 years
    OIClimMean(:,:,mm) = mean(OISST(:,:,ClimIdx),3,'omitnan').*mask;
    OIClimSdev(:,:,mm) = std(OISST(:,:,ClimIdx),[],3,'omitnan').*mask;
    clear ClimIdx mask
end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_savepngL(figh,filename)
% function fun_savepngL
% save figure in landscape


set(figh,'paperorientation','landscape','papertype','a4','paperpositionmode','auto',...
    'paperunits','centimeters','paperposition',[0.25 0.25 (29.7-0.5) (21-0.5)])
print(figh, '-dpng', '-r300', filename)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bplot_coastGLOB(plotspec)
%
% BPLOT_COASTGLOB
%
% Plots WORLD coastline on current figure.  If no input
% arguments, the function will plot the coastline as a black line.
%
% Coastline file from m_map toolbox - https://www-old.eoas.ubc.ca/~rich/map.html

if nargin<1
    plotspec = '-k';
end

load(['.\m_coasts.mat'])

figure(gcf)
hold on
plot(ncst(:,1),ncst(:,2),plotspec)
end