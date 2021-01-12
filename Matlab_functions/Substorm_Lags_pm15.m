load('Indicies_1997to2001.mat') %year,month,day,hour,minute,sec,SME,SML,SMU,SMEd,SMLd,SMUd,SMR,Bx,By,Bz,Vx,Vy,Vz 
load('substorm_info_116.mat') %Contains onset and peak timings for each substorm
for storm_num=[1:116]
    Substorm_Lags_1(storm_num,Indicies_4yrs,s_info,'OFF');
    storm_num
    datetime
end
for storm_num=[1,13,18,29,33]
    for surr=1:10
        Substorm_Lags_1(storm_num,Indicies_4yrs,s_info,'ON ',surr);
        storm_num
        datetime
        surr
    end
end

function[]=Substorm_Lags_1(storm_num,Indicies_4yrs,s_info,SURR,surrno)
    %Takes the information from the per day processed_w128 structures and
    %saves the information relevent to the substorm as a structure
    [Substorm]=relevent_information(storm_num,s_info,SURR);
    
    %Set up and calculate matrix for cross correlation matrix (no. stats x no.
    %stats x no. lags x length substorm), max lag cross correlation matrix and max lag values (no. stats x no.
    %stats x length substorm),
    CanonXCorr=nan(Substorm.nstats,Substorm.nstats,37,Substorm.length);
    Max_CXC=nan(Substorm.nstats,Substorm.nstats,Substorm.length);
    Max_Lags=nan(Substorm.nstats,Substorm.nstats,Substorm.length);
    for i=19:Substorm.length+18
        [CanonXCorr(:,:,:,i-18),Max_CXC(:,:,i-18),Max_Lags(:,:,i-18)]=canonical_matrix_lagged(Substorm,i);
    end
    
    % Add correlation information to substorm structure
    Substorm.CanonXCorr=CanonXCorr;
    Substorm.Max_CXC=Max_CXC;
    Substorm.Max_Lags=Max_Lags;

    %Include indicies in substorm structure
    TIMES=Substorm.datetimes(1)-minutes(128+18):minutes(1):Substorm.datetimes(end)+minutes(18);
    [x,~]=find(Indicies_4yrs(:,1)==year(TIMES) & Indicies_4yrs(:,2)==month(TIMES) & Indicies_4yrs(:,3)==day(TIMES) & Indicies_4yrs(:,4)== hour(TIMES) & Indicies_4yrs(:,5)==minute(TIMES));
    Substorm.indicies=Indicies_4yrs(x,[7:9,103]); %SME,SML,SMU,SMR
    %Calculate Separation matrix
    Substorm.Geodetic_dist=distanceseparation(Substorm);
    
    if SURR=='OFF'
        save(['Substorms_v3/',num2str(storm_num),'.substorm_v3.mat'],'Substorm','-v7.3')
    else
        save(['Substorms_v3_surrogates/',num2str(storm_num),'.substorm_v3',num2str(surrno),'.mat'],'Substorm','-v7.3')
    end
end

function[cxc,Max_CXC,Max_Lags]=canonical_matrix_lagged(struct,t1)
%Detrends and calculated the canonical cross correlation for lagged time
%series. Outputs zero lag cross-correlation matrix (CCM), max lag CCM and
%a matrix of the corresponding maximum lags
    time=128;
    t2=t1+127;
    lags=18;
    nstat=struct.nstats;
    detrended=zeros(time,nstat,3,2*lags+1);
    for l=-lags:1:lags
        detrended(:,:,1,lags+l+1)=detrend(struct.n_nez(:,t1+l:t2+l)');
        detrended(:,:,2,lags+l+1)=detrend(struct.e_nez(:,t1+l:t2+l)');
        detrended(:,:,3,lags+l+1)=detrend(struct.z_nez(:,t1+l:t2+l)');
    end
    
    act_stat=sum(sum(isnan(squeeze(detrended(:,:,:,19))),3),1)==0;
    no_act_stat=sum(act_stat);
    act_stat=find(act_stat==1);
    deadstations=isnan(squeeze(sum(sum(sum(detrended,1),3),4)));

    
    cxc=nan(nstat,nstat,2*lags+1);
    Max_CXC=nan(nstat,nstat);
    Max_Lags=nan(nstat,nstat);
    
    for i=1:no_act_stat
        d1=squeeze(detrended(:,act_stat(1,i),:,19));
        for j=1:no_act_stat
            if deadstations(act_stat(1,j))==0
                for Lg=1:lags*2+1
                     d2=squeeze(detrended(:,act_stat(1,j),:,Lg));                
                     [~,~,r]=canoncorr(d1,d2);
                     cxc(act_stat(1,i),act_stat(1,j),Lg)=r(1);
                end
            XC=squeeze(cxc(act_stat(1,i),act_stat(1,j),:));
            [Max_CXC(act_stat(1,i),act_stat(1,j)),Max_Lags(act_stat(1,i),act_stat(1,j))]=Peak_Finder(XC);
            end
        end
    end      
end

function[Max_CXC,Max_Lags]=Peak_Finder(canxcor)
%Finds the peak of lagged canonical correlation. Conditions are set such that it must be a true peak and not simply a continuous upward trend or flat. 
    [maxi,mind]=max(canxcor(4:16+18,1));
    mind=mind+3;
    if mind==19
        Max_CXC=maxi;
        Max_Lags=mind;
    else
            [~,pind]=findpeaks(canxcor(3:35,1));
            peaks=-1;
            pind=pind+2;
            for i=pind'
                if canxcor(i-2,1)<=canxcor(i-1,1) && canxcor(i+2,1)<=canxcor(i+1,1) && canxcor(i,1)>peaks || canxcor(i-3,1)<=canxcor(i-1,1) && canxcor(i+3,1)<=canxcor(i+1,1) && canxcor(i,1)>peaks ||canxcor(i-2,1)<=canxcor(i,1) && canxcor(i-3,1)<canxcor(i,1) && canxcor(i+2,1)<=canxcor(i,1) && canxcor(i+3,1)<canxcor(i,1) && canxcor(i,1)>peaks  
                    peaks=canxcor(i,1);
                    indpeaks=i;
                end
            end
         if peaks==maxi
            Max_CXC=peaks;
            Max_Lags=indpeaks;
         elseif peaks>canxcor(19,1)
            Max_CXC=peaks;
            Max_Lags=indpeaks;
         else
            Max_CXC=canxcor(19,1);
            Max_Lags=19;
         end   

    end     
        Max_Lags=Max_Lags-19;
end


function[substorm]=relevent_information(storm_num,s_info,SURR)
    %Loads processed_w128 files per day and outputs the information relevent
    %to the substorm as a structure called 'substorm'. It saves timings from 5xlength of expansion phase.
    substorm_info.onset=s_info.onset(storm_num,1);
    substorm_info.peak=s_info.peak(storm_num,1);
    sub_start=substorm_info.onset-2*(substorm_info.peak-substorm_info.onset);
    sub_fin=substorm_info.peak+2*(substorm_info.peak-substorm_info.onset);
    substorm_times=sub_start:minutes(1):sub_fin;
    time_window=[hour(sub_start)*60+minute(sub_start)+1,hour(sub_fin)*60+minute(sub_fin)+1];
    
    %Length of substorm based on if they span one day or two
    length_substorm=length(substorm_times)+127+2*18;
    start=sub_start-minutes(127+18);
    start_min=hour(start)*60+minute(start)+1;
    start_year=year(start);
    finish=sub_fin+minutes(18);
    finish_min=hour(finish)*60+minute(finish)+1;
    finish_year=year(finish); 
    substorm_info.year=start_year;
    
    if time_window(1)<=(127+18) 
            load(['supermag/processed_w128/',num2str(substorm_info.year),'/',num2str(substorm_info.year),num2str(month(substorm_info.onset),'%02.f'),num2str(day(substorm_info.onset),'%02.f'),'.processed.mat'])
            main_day=struct2;
            clear('struct2');
            load(['supermag/processed_w128/',num2str(start_year),'/',num2str(start_year),num2str(month(start),'%02.f'),num2str(day(start),'%02.f'),'.processed.mat'])
            prev_day=struct2;
            clear('struct2');
            
            num_prev=-1*(time_window(1)-(128+18));

            if length([start_min:1440,1:finish_min])~=length([1:num_prev,num_prev+1:length_substorm])
                error('Lengths different')
            end
                
                substorm.time=nan(length_substorm,5);
                substorm.active=nan(main_day.nstats,length_substorm);
                substorm.mlat=nan(main_day.nstats,length_substorm);
                substorm.mlt=nan(main_day.nstats,length_substorm);
                substorm.glat=nan(main_day.nstats,length_substorm);
                substorm.glon=nan(main_day.nstats,length_substorm);
                substorm.names=main_day.names;
                substorm.n_nez=nan(main_day.nstats,length_substorm);
                substorm.e_nez=nan(main_day.nstats,length_substorm);
                substorm.z_nez=nan(main_day.nstats,length_substorm);
                substorm.nstats=main_day.nstats;
                substorm.daylit=nan(main_day.nstats,length_substorm);
                
                substorm.time(num_prev+1:end,:)=main_day.time(1:finish_min,:);
                substorm.active(:,num_prev+1:end)=main_day.active(:,1:finish_min);
                substorm.mlat(:,num_prev+1:end)=main_day.mlat(:,1:finish_min);
                substorm.mlt(:,num_prev+1:end)=main_day.mlt(:,1:finish_min);
                substorm.glat(:,num_prev+1:end)=main_day.glat(:,1:finish_min);
                substorm.glon(:,num_prev+1:end)=main_day.glon(:,1:finish_min);             
                substorm.n_nez(:,num_prev+1:end)=main_day.n_nez(:,1:finish_min);
                substorm.e_nez(:,num_prev+1:end)=main_day.e_nez(:,1:finish_min);
                substorm.z_nez(:,num_prev+1:end)=main_day.z_nez(:,1:finish_min);
                substorm.daylit(:,num_prev+1:end)=main_day.daylit(:,1:finish_min);
                
                index=nan(prev_day.nstats,main_day.nstats);
                    for i=1:main_day.nstats
                        index(:,i)=strcmp(main_day.names(i,:),prev_day.names);
                    end
                [r,c]=find(index==1);   
                
                substorm.time(1:num_prev,:)=prev_day.time(start_min:end,:);
                substorm.active(c,1:num_prev)=prev_day.active(r,start_min:end);
                substorm.mlat(c,1:num_prev)=prev_day.mlat(r,start_min:end);
                substorm.mlt(c,1:num_prev)=prev_day.mlt(r,start_min:end);
                substorm.glat(c,1:num_prev)=prev_day.glat(r,start_min:end);
                substorm.glon(c,1:num_prev)=prev_day.glon(r,start_min:end);     
                substorm.n_nez(c,1:num_prev)=prev_day.n_nez(r,start_min:end);
                substorm.e_nez(c,1:num_prev)=prev_day.e_nez(r,start_min:end);
                substorm.z_nez(c,1:num_prev)=prev_day.z_nez(r,start_min:end);
                substorm.daylit(c,1:num_prev)=prev_day.daylit(r,start_min:end);
                
                substorm.length=length_substorm-127-2*18;
                substorm.datetimes=substorm_times;
                
    elseif day(start)~=day(finish)
            load(['supermag/processed_w128/',num2str(substorm_info.year),'/',num2str(substorm_info.year),num2str(month(sub_start),'%02.f'),num2str(day(sub_start),'%02.f'),'.processed.mat'])
            main_day=struct2;
            clear('struct2');
    
            load(['supermag/processed_w128/',num2str(finish_year),'/',num2str(finish_year),num2str(month(finish),'%02.f'),num2str(day(finish),'%02.f'),'.processed.mat'])
            second_day=struct2;
            clear('struct2');

            length_mainday=length(start_min:1440);
            if length([start_min:1440,1:finish_min])~=length_substorm
                error('Lengths are not equal')
            end
                substorm.time=nan(length_substorm,5);
                substorm.active=nan(main_day.nstats,length_substorm);
                substorm.mlat=nan(main_day.nstats,length_substorm);
                substorm.mlt=nan(main_day.nstats,length_substorm);
                substorm.glat=nan(main_day.nstats,length_substorm);
                substorm.glon=nan(main_day.nstats,length_substorm);
                substorm.names=main_day.names;
                substorm.n_nez=nan(main_day.nstats,length_substorm);
                substorm.e_nez=nan(main_day.nstats,length_substorm);
                substorm.z_nez=nan(main_day.nstats,length_substorm);
                substorm.nstats=main_day.nstats;
                substorm.daylit=nan(main_day.nstats,length_substorm);
                
                substorm.time(1:length_mainday,:)=main_day.time(start_min:1440,:);
                substorm.active(:,1:length_mainday)=main_day.active(:,start_min:1440);
                substorm.mlat(:,1:length_mainday)=main_day.mlat(:,start_min:1440);
                substorm.mlt(:,1:length_mainday)=main_day.mlt(:,start_min:1440);
                substorm.glat(:,1:length_mainday)=main_day.glat(:,start_min:1440);
                substorm.glon(:,1:length_mainday)=main_day.glon(:,start_min:1440);             
                substorm.n_nez(:,1:length_mainday)=main_day.n_nez(:,start_min:1440);
                substorm.e_nez(:,1:length_mainday)=main_day.e_nez(:,start_min:1440);
                substorm.z_nez(:,1:length_mainday)=main_day.z_nez(:,start_min:1440);
                substorm.daylit(:,1:length_mainday)=main_day.daylit(:,start_min:1440);
                
                index=nan(second_day.nstats,main_day.nstats);
                    for i=1:main_day.nstats
                        index(:,i)=strcmp(main_day.names(i,:),second_day.names);
                    end
                [r,c]=find(index==1);   
                
                substorm.time(length_mainday+1:end,:)=second_day.time(1:finish_min,:);
                substorm.active(c,length_mainday+1:end)=second_day.active(r,1:finish_min);
                substorm.mlat(c,length_mainday+1:end)=second_day.mlat(r,1:finish_min);
                substorm.mlt(c,length_mainday+1:end)=second_day.mlt(r,1:finish_min);
                substorm.glat(c,length_mainday+1:end)=second_day.glat(r,1:finish_min);
                substorm.glon(c,length_mainday+1:end)=second_day.glon(r,1:finish_min);     
                substorm.n_nez(c,length_mainday+1:end)=second_day.n_nez(r,1:finish_min);
                substorm.e_nez(c,length_mainday+1:end)=second_day.e_nez(r,1:finish_min);
                substorm.z_nez(c,length_mainday+1:end)=second_day.z_nez(r,1:finish_min);
                substorm.daylit(c,length_mainday+1:end)=second_day.daylit(r,1:finish_min);
                
                substorm.length=length_substorm-127-2*18;
                substorm.datetimes=substorm_times;
                
    elseif day(start)==day(finish)
                
                load(['supermag/processed_w128/',num2str(substorm_info.year),'/',num2str(substorm_info.year),num2str(month(start),'%02.f'),num2str(day(start),'%02.f'),'.processed.mat'])
                main_day=struct2;
                clear('struct2');
                
                substorm.time=main_day.time(start_min:finish_min,:);
                substorm.active=main_day.active(:,start_min:finish_min);
                substorm.mlat=main_day.mlat(:,start_min:finish_min);
                substorm.mlt=main_day.mlt(:,start_min:finish_min);
                substorm.glat=main_day.glat(:,start_min:finish_min);
                substorm.glon=main_day.glon(:,start_min:finish_min);     
                substorm.n_nez=main_day.n_nez(:,start_min:finish_min);
                substorm.e_nez=main_day.e_nez(:,start_min:finish_min);
                substorm.z_nez=main_day.z_nez(:,start_min:finish_min);
                substorm.daylit=main_day.daylit(:,start_min:finish_min);
                substorm.nstats=main_day.nstats;
                substorm.names=main_day.names;
                substorm.length=length_substorm-127-2*18;
                substorm.datetimes=substorm_times;
    else
        error('Something odd is going on')
    end
    %Exclude stations which are only outside our range of interest (50-85
    %degrees lat, nightside, in darkness)
     High_stats=zeros(substorm.nstats,substorm.length);
     Night_stats=zeros(substorm.nstats,substorm.length);
    for i=19:substorm.length+18
            high_stats=substorm.mlat(:,i+127);
            high_stats(high_stats<60)=nan;
            high_stats(high_stats>75)=nan;
            high_stats(high_stats>=60 & high_stats<=75)=1;
            High_stats(:,i-18)=high_stats;  
            night_stats=substorm.mlt(:,i+127);
            night_stats(night_stats<18 &night_stats>6)=nan;
            night_stats(night_stats>=18)=1;
            night_stats(night_stats<=6)=1;
            Night_stats(:,i-18)=night_stats;  
    end

    Dark_stats=zeros(substorm.nstats,substorm.length);
    for i=19:substorm.length+18
    dark_stats=sum(substorm.daylit(:,i-18:i+127+18),2)/164;
    Dark_stats(:,i-18)=dark_stats;
    end
    Excluded_stats=(Dark_stats+High_stats+Night_stats)./3;
    Excluded_stats=nansum(Excluded_stats,2);
    Excluded_stats=find(Excluded_stats==0);
    
    substorm.active(Excluded_stats,:)=[];
    substorm.mlat(Excluded_stats,:)=[];
    substorm.mlt(Excluded_stats,:)=[];
    substorm.glat(Excluded_stats,:)=[];
    substorm.glon(Excluded_stats,:)=[];     
    substorm.n_nez(Excluded_stats,:)=[];
    substorm.e_nez(Excluded_stats,:)=[];
    substorm.z_nez(Excluded_stats,:)=[];
    substorm.daylit(Excluded_stats,:)=[];
    substorm.nstats=substorm.nstats-size(Excluded_stats,1);
    substorm.names(Excluded_stats,:)=[];
    
    substorm.timings=30.*(substorm.datetimes-substorm_info.onset)./(substorm_info.peak-substorm_info.onset);
    substorm.info=substorm_info;
    
    if SURR=='ON '
        for k=1:size(substorm.n_nez,1)
            substorm.n_nez(k,:)=IAAFTsur(substorm.n_nez(k,:));
            substorm.e_nez(k,:)=IAAFTsur(substorm.e_nez(k,:));
            substorm.z_nez(k,:)=IAAFTsur(substorm.z_nez(k,:));
        end
    end
end


function [ Distance ] = distanceseparation(Substorm)
    %Find the geodestic distance separations between all pairs of stations
    nstats=Substorm.nstats;
    Long=Substorm.glon(:,1);
    Lat=Substorm.glat(:,1);
    Distance=zeros(nstats,nstats);
    
    for i=1:nstats
        for j=i+1:nstats
            Distance(i,j)=vdist(Lat(i,1),Long(i,1),Lat(j,1),Long(j,1));
            Distance(j,i)=Distance(i,j);
        end

    end

end

function s = vdist(lat1,lon1,lat2,lon2)
% Michael Kleder (2017). Vectorized geodetic distance and azimuth on the WGS84 
% earth ellipsoid (https://www.mathworks.com/matlabcentral/fileexchange/8607-vectorized
% -geodetic-distance-and-azimuth-on-the-wgs84-earth-ellipsoid), MATLAB Central File Exchange. Retrieved Feburary 13, 2017.
% VDIST - compute distance between points on the WGS-84 ellipsoidal Earth
%         to within a few millimeters of accuracy using Vincenty's algorithm
%
% s = vdist(lat1,lon1,lat2,lon2)
%
% s = distance in meters
% lat1 = GEODETIC latitude of first point (degrees)
% lon1 = longitude of first point (degrees)
% lat2, lon2 = second point (degrees)
%
%  Original algorithm source:
%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%  April 1975, pp 88-93
%
% Notes: (1) Error correcting code, convergence failure traps, antipodal corrections,
%            polar error corrections, WGS84 ellipsoid parameters, testing, and comments
%            written by Michael Kleder, 2004.
%        (2) Vincenty describes his original algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (3) Essentially antipodal points are treated as exactly antipodal,
%            potentially reducing accuracy by a small amount.
%        (4) Failures for points exactly at the poles are eliminated by
%            moving the points by 0.6 millimeters.
%        (5) Vincenty's azimuth formulas are not implemented in this
%            version, but are available as comments in the code.
%        (6) The Vincenty procedure was transcribed verbatim by Peter Cederholm,
%            August 12, 2003. It was modified and translated to English by Michael Kleder.
%            Mr. Cederholm's website is http://www.plan.aau.dk/~pce/
%        (7) Code to test the disagreement between this algorithm and the
%            Mapping Toolbox spherical earth distance function is provided
%            as comments in the code. The maximum differences are:
%            Max absolute difference: 38 kilometers
%            Max fractional difference: 0.56 percent

% Input check:
if abs(lat1)>90 | abs(lat2)>90
    error('Input latitudes must be between -90 and 90 degrees, inclusive.')
end
% Supply WGS84 earth ellipsoid axis lengths in meters:
a = 6378137; % definitionally
b = 6356752.31424518; % computed from WGS84 earth flattening coefficient definition
% convert inputs in degrees to radians:
lat1 = lat1 * 0.0174532925199433;
lon1 = lon1 * 0.0174532925199433;
lat2 = lat2 * 0.0174532925199433;
lon2 = lon2 * 0.0174532925199433;
% correct for errors at exact poles by adjusting 0.6 millimeters:
if abs(pi/2-abs(lat1)) < 1e-10;
    lat1 = sign(lat1)*(pi/2-(1e-10));
end
if abs(pi/2-abs(lat2)) < 1e-10;
    lat2 = sign(lat2)*(pi/2-(1e-10));
end
f = (a-b)/a;
U1 = atan((1-f)*tan(lat1));
U2 = atan((1-f)*tan(lat2));
lon1 = mod(lon1,2*pi);
lon2 = mod(lon2,2*pi);
L = abs(lon2-lon1);
if L > pi
    L = 2*pi - L;
end
lambda = L;
lambdaold = 0;
itercount = 0;
while ~itercount | abs(lambda-lambdaold) > 1e-12  % force at least one execution
    itercount = itercount+1;
    if itercount > 50
        warning('Points are essentially antipodal. Precision may be reduced slightly.');
        lambda = pi;
        break
    end
    lambdaold = lambda;
    sinsigma = sqrt((cos(U2)*sin(lambda))^2+(cos(U1)*...
        sin(U2)-sin(U1)*cos(U2)*cos(lambda))^2);
    cossigma = sin(U1)*sin(U2)+cos(U1)*cos(U2)*cos(lambda);
    sigma = atan2(sinsigma,cossigma);
    alpha = asin(cos(U1)*cos(U2)*sin(lambda)/sin(sigma));
    cos2sigmam = cos(sigma)-2*sin(U1)*sin(U2)/cos(alpha)^2;
    C = f/16*cos(alpha)^2*(4+f*(4-3*cos(alpha)^2));
    lambda = L+(1-C)*f*sin(alpha)*(sigma+C*sin(sigma)*...
        (cos2sigmam+C*cos(sigma)*(-1+2*cos2sigmam^2)));
    % correct for convergence failure in the case of essentially antipodal points
    if lambda > pi
        warning('Points are essentially antipodal. Precision may be reduced slightly.');
        lambda = pi;
        break
    end
end
u2 = cos(alpha)^2*(a^2-b^2)/b^2;
A = 1+u2/16384*(4096+u2*(-768+u2*(320-175*u2)));
B = u2/1024*(256+u2*(-128+u2*(74-47*u2)));
deltasigma = B*sin(sigma)*(cos2sigmam+B/4*(cos(sigma)*(-1+2*cos2sigmam^2)...
    -B/6*cos2sigmam*(-3+4*sin(sigma)^2)*(-3+4*cos2sigmam^2)));
s = b*A*(sigma-deltasigma);

% % =====================================================================
% % Vicenty's azimuth calculation code is left unused:
% % (results in radians)
% % From point #1 to point #2
% a12 = atan2(cos(U2)*sin(lambda),cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda));
% if a12 < 0
%     a12 = a12+2*pi;
% end
% % from point #2 to point #1
% a21 = atan2(cos(U1)*sin(lambda),-sin(U1)*cos(U2)+cos(U1)*sin(U2)*cos(lambda));
% if a21 < 0
%     a21 = a21+pi;
% end
% if (L>0) & (L<pi)
%     a21 = a21 + pi;
% end

% % =====================================================================
% % Code to test the Mapping Toolbox spherical earth distance against
% % Vincenty's algorithm using random test points:
% format short g
% errmax=0;
% abserrmax=0;
% for i = 1:10000
%     llat = rand * 184-92;
%     tlat = rand * 184-92;
%     llon = rand * 364 - 2;
%     tlon = rand * 364 - 2;
%     llat = max(-90,min(90,llat)); % round to include occasional exact poles
%     tlat = max(-90,min(90,tlat));
%     llon = max(0,min(360,llon));
%     tlon = max(0,min(360,tlon));
%     % randomly test exact equator
%     if rand < .01
%         llat = 0;
%         llon = 0;
%     else
%         if rand < .01
%             llat = 0;
%         end
%         if rand < .01
%             tlat = 0;
%         end
%     end
%     dm = 1000*deg2km(distance(llat,llon,tlat,tlon));
%     dv = vdist(llat,llon,tlat,tlon);
%     abserr = abs(dm-dv);
%     if abserr < 1e-2 % disagreement less than a centimeter
%         err = 0;
%     else
%         err = abs(dv-dm)/dv;
%     end
%     errmax = max(err,errmax);
%     abserrmax = max(abserr,abserrmax);
%     %     if i==1 | rand > .99
%     disp([i dv dm err errmax abserrmax])
%     %     end
%     if err > .01
%         break
%     end
% end
end
function zM = IAAFTsur(xV,nsur)
% zM = IAAFTsur(xV,nsur)
% IAAFT: Iterated Amplitude Adjusted Fourier Transform surrogates
% This function generates 'nsur' IAAFT-surrogate time series 
% and stores them in the matrix 'zM' (columnwise). These surrogates
% are supposed to have the same amplitude distribution (marginal cdf) 
% and autocorrelation as the given time series 'xV'. 
% The IAAFT algorithm is proposed in 
% Schreiber, T. and Schmitz, A. (1996) "Improved Surrogate Data for 
% Nonlinearity Tests", Physical Review Letters, Vol 77, 635-638.
% The IAAFT is an improvement of the AAFT. Iteratively, it fits the 
% amplitudes and at each step improves the spectral phases and then 
% reorders the derived time series at each step until convergence of 
% both spectral density and amplitude distribution is reached. 
% The algorithm terminates if complete convergence (same reordering in 
% two consecutive steps) is succeeded or if the 'maxi' number of 
% iterations is reached. 
% INPUT
% - xV  : the given time series
% - nsur: the number of surrogate time series (default is 1)
% OUTPUT
% - zM  : the n x nsur matrix of 'nsur' IAAFT surrogate time series
%========================================================================
%     <IAAFTsur.m>, v 1.0 2010/02/11 22:09:14  Kugiumtzis & Tsimpiris
%     This is part of the MATS-Toolkit http://eeganalysis.web.auth.gr/

%========================================================================
% Copyright (C) 2010 by Dimitris Kugiumtzis and Alkiviadis Tsimpiris 
%                       <dkugiu@gen.auth.gr>

%========================================================================
% Version: 1.0

% LICENSE:
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program. If not, see http://www.gnu.org/licenses/>.

%=========================================================================
% Reference : D. Kugiumtzis and A. Tsimpiris, "Measures of Analysis of Time Series (MATS): 
% 	          A Matlab  Toolkit for Computation of Multiple Measures on Time Series Data Bases",
%             Journal of Statistical Software, in press, 2010

% Link      : http://eeganalysis.web.auth.gr/
%========================================================================= 
maxi = 1000;
if nargin == 1
    nsur = 1;
end
n = length(xV);
zM = NaN*ones(n,nsur);
if rem(n,2) == 0
    n2 = n/2;
else
    n2 = (n-1)/2;
end
% Fourier transform of the original and smoothing of spectrum
zV = fft(xV);
S = abs(zV); 

for isur=1:nsur
    % permutation of the original series 'xV'
    rV = xV(randperm(n));
    % First comparison of spectrum of 'rV' to that of 'xV' using the 
    % smoothed spectra.
    wV = fft(rV);
    RS = abs(wV);
    rphV = angle(wV);
    k=1;
    [xoV, iVold] = sort(xV); 
    converge = 0;
    while k<=maxi & converge == 0 
        wwV = S.*exp(rphV.*i); 
        tmpV = real(ifft(wwV));
        [tmpV, indV] = sort(tmpV);
        [tmpV,iVnew] = sort(indV);
        rV = xoV(iVnew);
        wV = fft(rV);
        rphV = angle(wV);
        if iVnew == iVold
            converge = 1;
        else
            iVold = iVnew;
            k=k+1;
        end
    end
    zM(:,isur) = rV;
end
end