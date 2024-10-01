%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Script to load Hadley Centre HADISST v1 Sea Surface Temperature Data,%
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
lon_extent = [-80 50];
lat_extent = [0 90];

HADSST_datafolder = ['I:\Data_External\HADISST1\Datafiles\RAW\'];

clim_ref_period = [1991 2020];

%% set folder structure for output
if ~exist([pwd,'\HADSST_AnomSST'],'dir')
    mkdir(['.\HADSST_AnomSST\'])
end
if ~exist([pwd,'\HADSST_NormAnomSST'],'dir')
    mkdir(['.\HADSST_NormAnomSST\'])
end

%% Load HADSST data
[HADSST,HADSST_time,HADSST_lon,HADSST_lat] = ...
    fun_get_Hadley_SST_timeseries(HADSST_datafolder,[datenum(1870,1,1),floor(now)],...
    lon_extent,lat_extent);
HADSST_tvec = datevec(HADSST_time);
HADSST(HADSST<-100)=NaN;

%% Calculate Climatology
[HADClimMean,HADClimSdev] = fun_get_Hadley_SST_climatology(HADSST_tvec,HADSST,clim_ref_period);

%% Calculate anomalies and normalised anomalies
HAD_Anom = NaN.*zeros(size(HADSST,1),size(HADSST,2),size(HADSST,3));
HAD_NormAnom = NaN.*zeros(size(HADSST,1),size(HADSST,2),size(HADSST,3));
for mm=1:12;
    idxHAD =  find(HADSST_tvec(:,2)==mm);
    HAD_Anom(:,:,idxHAD) = HADSST(:,:,idxHAD) - repmat(HADClimMean(:,:,mm),1,1,length(idxHAD));
    HAD_NormAnom(:,:,idxHAD) = (HADSST(:,:,idxHAD) - repmat(HADClimMean(:,:,mm),1,1,length(idxHAD)))./ ...
        repmat(HADClimSdev(:,:,mm),1,1,length(idxHAD));
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
    HAD_Anom_select = HAD_Anom(:,:,HADSST_tvec(:,1)==selectyear(yy));
    HAD_NormAnom_select = HAD_NormAnom(:,:,HADSST_tvec(:,1)==selectyear(yy));


    %% mapping parameters
    %cmap = flipud(cbrewer('div','RdYlBu',22));
    cmap = cat(1,flipud(cbrewer('seq','PuBu',10)),cbrewer('seq','YlOrRd',10));
    cmap(cmap>1)=1;cmap(cmap<0)=0;

    pos12 = [0.02 0.68 0.225 0.3;0.255 0.68 0.225 0.3;0.49 0.68 0.225 0.3;0.725 0.68 0.225 0.3;...
        0.02 0.36 0.225 0.3;0.255 0.36 0.225 0.3;0.49 0.36 0.225 0.3;0.725 0.36 0.225 0.3;...
        0.02 0.04 0.225 0.3;0.255 0.04 0.225 0.3;0.49 0.04 0.225 0.3;0.725 0.04 0.225 0.3];


    %% ANOM SST HAD SST relative HADSST
    close all;
    figure(1)
    for mm=1:12
        [tim_chck]=intersect(HADSST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(HADSST_lon(:))-0.025;HADSST_lon+0.025];
        lat2plot = [min(HADSST_lat(:))-0.025;HADSST_lat+0.025];
        heat2plot = HAD_Anom_select(:,:,mm);

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
    fun_savepngL(gcf,['.\HADSST_AnomSST\SST_Maps_HADSST_AnomSST_',num2str(selectyear(yy)),'_NAtl.png'])

    %% Standardised ANOM SST HAD SST relative HADSST
    close all;
    figure(1)
    for mm=1:12
        [tim_chck]=intersect(HADSST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(HADSST_lon(:))-0.025;HADSST_lon+0.025];
        lat2plot = [min(HADSST_lat(:))-0.025;HADSST_lat+0.025];
        heat2plot = HAD_NormAnom_select(:,:,mm);

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
    fun_savepngL(gcf,['.\HADSST_NormAnomSST\SST_Maps_HADSST_NormAnomSST_',num2str(selectyear(yy)),'_NAtl.png'])


    %% ANOM SST HAD SST relative HADSST  -SPNA only
    close all;
    figure(1)

    for mm=1:12
        [tim_chck]=intersect(HADSST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(HADSST_lon(:))-0.025;HADSST_lon+0.025];
        lat2plot = [min(HADSST_lat(:))-0.025;HADSST_lat+0.025];
        heat2plot = HAD_Anom_select(:,:,mm);

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
    fun_savepngL(gcf,['.\HADSST_AnomSST\SST_Maps_HADSST_AnomSST_',num2str(selectyear(yy)),'_SPNA.png'])

    %% ANOM SST HAD SST relative HADSST  -SPNA only
    close all;
    figure(1)

    for mm=1:12
        [tim_chck]=intersect(HADSST_tvec(:,1:2),[selectyear(yy),mm],'rows');
        if isempty(tim_chck);clear tim_chck;continue;else;clear tim_chck;end;

        subplot(3,4,mm)

        lon2plot = [min(HADSST_lon(:))-0.025;HADSST_lon+0.025];
        lat2plot = [min(HADSST_lat(:))-0.025;HADSST_lat+0.025];
        heat2plot = HAD_NormAnom_select(:,:,mm);

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
    fun_savepngL(gcf,['.\HADSST_NormAnomSST\SST_Maps_HADSST_NormAnomSST_',num2str(selectyear(yy)),'_SPNA.png'])
end
return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS SUPPORTING THE ABOVE CODE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SST,SST_TimeFromFile,lon,lat] = fun_get_Hadley_SST_timeseries(data_folder,time_choice,lon_choice,lat_choice)
% function fun_get_Hadley_SST_timeseries
%
% Extract time series of Hadley SST data from the UK Met Office HADSST1 product.
% Requires most up to date netcdf files to be saved in data_folder
%
%
% USE:
%       [SST,time,lon,lat] = fun_get_Hadley_SST_timeseries(time_choice,lon_choice,lat_choice)
%
% INPUT:
%       time_choice   = start and end time for field (matlab datenum format)
%       lon_choice    = western & eastern extent of longitude for field [-180 180]
%       lat_choice    = southern & northern extent of longitude for field [-90 90]
%
%
% OUTPUT:
%       SST (time x 1)        = monthly mean SST (deg C) from HADSST
%       time                  = time stamp from datafile
%       lon                   = longitude of grid extracted (nearest long)
%       lat                   = latitude of grid extracted (nearest lat)
%
% DEPENDENCIES:
%   The function needs access to the following
%       DATA: data saved in HADSST format (not compressed) in data_folder
%
%
% EXAMPLE:
%       [SST,time,lon,lat] = fun_get_Hadley_SST_timeseries(data_folder, ...
%               [datenum(2022,1,1),floor(now)],[-25 15],[45 65]);
%

HADSSTDir = data_folder;
fname = {[HADSSTDir,'HadISST_sst.nc'];[HADSSTDir,'HadISST1_SST_update.nc']};

if length(lon_choice)==1
    lon_choice = [lon_choice lon_choice];
end

if length(lat_choice)==1
    lat_choice = [lat_choice lat_choice];
end
lat_choice = [max(lat_choice) min(lat_choice)];

time_request = datevec([min(time_choice):1:max(time_choice)]);
time_request = unique(time_request(:,1:2),'rows','stable');

lat = ncread(fname{1},'latitude');
lon = ncread(fname{1},'longitude');
lon(lon>180)=lon(lon>180)-360;

[~,idxlon] = min(abs(lon-lon_choice));
[~,idxlat] = min(abs(lat-lat_choice));
clear lon lat
lon  = ncread(fname{1},'longitude',[idxlon(1)],[idxlon(2)-idxlon(1)+1]);
lat  = ncread(fname{1},'latitude',[idxlat(1)],[idxlat(2)-idxlat(1)+1]);

SST_TimeFromFile = NaN.*time_request(:,1);
SST = NaN.*zeros(size(lon,1),size(lat,1),size(time_request,1));

for ff=1:length(fname)
    time_orig = datevec(double(datenum(1870,1,1)+ncread(fname{ff},'time')));
    [~,itime,~] = intersect(time_orig(:,1:2),time_request(:,1:2),'rows');
    if isempty(itime);continue;end
    sst_in = ncread(fname{ff},'sst',[idxlon(1),idxlat(1),min(itime)],...
        [idxlon(2)-idxlon(1)+1,idxlat(2)-idxlat(1)+1,max(itime)-min(itime)+1]);
    [~,ia,ib]=intersect(time_orig(:,1:2),time_request(:,1:2),'rows');
    sst_in(sst_in==-1.000000015047466e+30)=NaN;
    SST(:,:,ib)=squeeze(sst_in);
    SST_TimeFromFile(ib,1) = datenum(time_orig(ia,:));
    clear time_orig itime sst_in ia ib
end
SST(:,:,isnan(SST_TimeFromFile))=[];
SST_TimeFromFile(isnan(SST_TimeFromFile))=[];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HADClimMean,HADClimSdev] = fun_get_Hadley_SST_climatology(HADSST_tvec,HADSST,clim_ref_period)
% sst climatology - mean and st dev in ref period

HADClimMean = NaN.*zeros(size(HADSST,1),size(HADSST,2),12);
HADClimSdev = NaN.*zeros(size(HADSST,1),size(HADSST,2),12);
for mm=1:12
    ClimIdx  = intersect(intersect(find(HADSST_tvec(:,1)>=clim_ref_period(1)),...
        find(HADSST_tvec(:,1)<=clim_ref_period(2))),...
        find(HADSST_tvec(:,2)==mm));
    HADClimMean(:,:,mm) = mean(HADSST(:,:,ClimIdx),3);
    HADClimSdev(:,:,mm) = std(HADSST(:,:,ClimIdx),[],3);
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