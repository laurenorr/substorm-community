%Code to read ncdf files into a matlab structure. Data is saved as year month day in separate year folders between 1996 to
%2001. Each day contains a normalization matrix of what proportion of the network would be connected if the correlation
%threshold was t where t=0:0.001:1. This matrix will later be used for a
%station and event specific threshold.
clear
w=128; %window length
n=1;
date_time_big=[datetime(1996,12,23):datetime(1997,04,13),datetime(1997,08,22):datetime(1998,04,01),datetime(1998,08,31):datetime(1999,04,10),datetime(1999,07,28):datetime(2000,04,01),datetime(2000,10,17):datetime(2001,03,09),datetime(2001,07,25):datetime(2001,11,17)];
time_periods=[1:1067:1068];
date_time=date_time_big(1,time_periods(n):time_periods(n+1)-1);
prevday=date_time(1,1)-days(1);
YR=year(prevday);
[struct,SZA]=firstday(prevday,w);
for q=1:length(date_time)
    cdate=date_time(1,q);
    if cdate~=prevday+days(1)
        prevday=cdate-days(1);
        YR=year(prevday);
        [struct,SZA]=firstday(prevday,w); 
    end
        cyr=year(cdate);
        if cyr~=YR
            YR=cyr;
            load(['SZA_data/SZA',num2str(cyr),'.mat']);
            if cyr==1996
                SZA=SZA1996;
                clear('SZA1996')
            elseif cyr==1997
                SZA=SZA1997;
                clear('SZA1997')
            elseif cyr==1998
                SZA=SZA1998;
                clear('SZA1998')
            elseif cyr==1999
                SZA=SZA1999;
                clear('SZA1999')
            elseif cyr==2000
                SZA=SZA2000;
                clear('SZA2000')
            elseif cyr==2001
                SZA=SZA2001;
                clear('SZA2001')
            end
        end
        date_num=[num2str(cyr),num2str(month(cdate),'%02.f'),num2str(day(cdate),'%02.f')];
        [struct]=Normalization(struct,date_num,cdate,SZA,w);
%         q
        prevday=cdate;
end


function[struct,SZA]=firstday(date,w)
%Normalization function for the firstday, when there is no previous day to
%load
    cyr=year(date);
    load(['SZA_data/SZA',num2str(cyr),'.mat']); 
    if cyr==1996
        SZA=SZA1996;
        clear('SZA1996')
    elseif cyr==1997
        SZA=SZA1997;
        clear('SZA1997')
    elseif cyr==1998
        SZA=SZA1998;
        clear('SZA1998')
    elseif cyr==1999
        SZA=SZA1999;
        clear('SZA1999')
    elseif cyr==2000
        SZA=SZA2000;
        clear('SZA2000')
    elseif cyr==2001
        SZA=SZA2001;
        clear('SZA2001')
    end
    [struct]=ncdf_conver(['supermag/raw/',num2str(year(date)),'/',num2str(year(date)),num2str(month(date),'%02.f'),num2str(day(date),'%02.f'),'.raw-supermag.ncdf']);
    [struct]=SZA_definer(date,SZA,struct);
     Norm_Sum=zeros(struct.nstats,1001);
     Norm_dark_Sum=Norm_Sum;
     Active_Sum=zeros(struct.nstats,1);
     Active_Dark_Sum=Active_Sum;
     w1=w-1;
     for t=1312-w1:1440-w1
 
         cxc=canonical_matrix(struct,t,w);
         [Norm_Dist,active_stats]=Normalizer(cxc);
         Norm_Dist(isnan(Norm_Dist))=0;
         Norm_Sum=Norm_Sum+Norm_Dist;
         Active_Sum=Active_Sum+active_stats;
 
         Dark_stats=repmat(sum(struct.daylit(:,t:t+w1),2)/w,1,struct.nstats);
         cxc2=cxc.*Dark_stats.*Dark_stats';
         [Norm_Dist,active_stats]=Normalizer(cxc2);
         Norm_Dist(isnan(Norm_Dist))=0;
         Norm_dark_Sum=Norm_dark_Sum+Norm_Dist;
         Active_Dark_Sum=Active_Dark_Sum+active_stats;
     end
 
     struct.Norm=Norm_Sum./repmat(Active_Sum,1,1001);
     struct.Norm_Dark=Norm_dark_Sum./repmat(Active_Dark_Sum,1,1001);
%      save(['supermag/processed_w',num2str(w),'/',num2str(year(date)),'/',num2str(num),'.processed.mat'],'struct','-v7.3')

end

function [structure] = ncdf_conver(data)
    %Convert ncdf supermag files into structures

    %Time- [year, month, day, hour, minute]
    structure.time=[ncread(data,'time_yr'),ncread(data,'time_mo'),ncread(data,'time_dy'),ncread(data,'time_hr'),ncread(data,'time_mt')];

    %We only want mid-high latitude stations
    mlat_unordered=ncread(data,'mlat');
    mlt_unordered=ncread(data,'mlt');
    glat_unordered=ncread(data,'glat');
    glon_unordered=ncread(data,'glon');
    n_nez_unordered=ncread(data,'dbn_nez');
    e_nez_unordered=ncread(data,'dbe_nez');
    z_nez_unordered=ncread(data,'dbz_nez');
    Lowest_mlat=40;

    %names
    unordered_name=ncread(data,'id');
    unordered_name=squeeze(string(permute(unordered_name(1:3,:,:),[2 1 3]))); %delete empty characters
    names=unique(unordered_name);
    no_names=length(names);
    names_ordered=zeros(size(unordered_name,1),size(unordered_name,2));
    names_indexed=nan(no_names,1440);
    for i=1:no_names
        names_ordered(:,:)=strcmp(names(i,1),unordered_name);
        [r,c]=find(names_ordered(:,:)==1);
         if length(unique(c))==length(c)
                names_indexed(i,c)=r;
         end

    end

    mlat_ordered=nan(no_names,1440);
    mlt_ordered=nan(no_names,1440);
    glat_ordered=nan(no_names,1440);
    glon_ordered=nan(no_names,1440);
    n_nez_ordered=nan(no_names,1440);
    e_nez_ordered=nan(no_names,1440);
    z_nez_ordered=nan(no_names,1440);
    for q=1:no_names
        for w=1:1440
            if names_indexed(q,w)>0
                mlat_ordered(q,w)=mlat_unordered(names_indexed(q,w),w);
                mlt_ordered(q,w)=mlt_unordered(names_indexed(q,w),w);
                glat_ordered(q,w)=glat_unordered(names_indexed(q,w),w);
                glon_ordered(q,w)=glon_unordered(names_indexed(q,w),w);
                n_nez_ordered(q,w)=n_nez_unordered(names_indexed(q,w),w);
                e_nez_ordered(q,w)=e_nez_unordered(names_indexed(q,w),w);
                z_nez_ordered(q,w)=z_nez_unordered(names_indexed(q,w),w);
            end
        end
    end

    lat_choice=ones(size(mlat_ordered));
    lat_choice(mlat_ordered<Lowest_mlat)=nan;
    lat_choice(isnan(mlat_ordered))=nan;

    %No data stations have a massive value
    lat_choice(mlat_ordered>90)=nan;

    %Rows with no data in them are deleted to save space
    Empty_Rows=all(isnan(lat_choice),2);
    
    structure.active=lat_choice;
    structure.active(Empty_Rows,:)=[];

    structure.mlat=mlat_ordered.*lat_choice;
    structure.mlat(Empty_Rows,:)=[];

    structure.mlt=mlt_ordered.*lat_choice;
    structure.mlt(Empty_Rows,:)=[];

    structure.glat=glat_ordered.*lat_choice;
    structure.glat(Empty_Rows,:)=[];

    structure.glon=glon_ordered.*lat_choice;
    structure.glon(Empty_Rows,:)=[];

    structure.n_nez=n_nez_ordered.*lat_choice;
    structure.n_nez(Empty_Rows,:)=[];

    structure.e_nez=e_nez_ordered.*lat_choice;
    structure.e_nez(Empty_Rows,:)=[];

    structure.z_nez=z_nez_ordered.*lat_choice;
    structure.z_nez(Empty_Rows,:)=[];

    structure.names=names;
    structure.names(Empty_Rows,:)=[];
    
    structure.nstats=size(structure.active,1);
    

    % checks
    if sum(ncread(data,'time_extent'))~=60*1440
        msg='Time Extent are not all 1 minute';
        error(msg)
    elseif sum(ncread(data,'time_sc'))~=0
        msg='Seconds not all equal to zero';
        error(msg)
    end
end

function[str]=SZA_definer(date,sza,str)
%Get rid of stations in daylight i.e. where the cosine of the zenith number is positive.
     str.doy=day(date,'dayofyear');
     no_names=length(str.names);
     index=nan(length(sza.names),no_names);
        for i=1:no_names
            index(:,i)=strcmp(str.names(i,:),sza.names);
        end
     [r,c]=find(index==1);
      SZA_data=nan(no_names,1440);   
      SZA_data(c,:)=sza.data(r,(str.doy-1)*1440+1:str.doy*1440);
      str.daylit=heaviside(-cos(deg2rad(SZA_data))); 
      str.daylit(str.daylit==0)=nan;
end

function[cxc]=canonical_matrix(struct,t1,w)
%Calculate canonical correlation at zero lag
    time=w; %length(t1:t2);
    w1=w-1;
    t2=t1+w1;
    nstat=struct.nstats;
    detrended=zeros(time,nstat,3);
    detrended(:,:,1)=detrend(struct.n_nez(:,t1:t2)');
    detrended(:,:,2)=detrend(struct.e_nez(:,t1:t2)');
    detrended(:,:,3)=detrend(struct.z_nez(:,t1:t2)');

    act_stat=sum(sum(isnan(detrended),3),1)==0;
    no_act_stat=sum(act_stat);
    act_stat=find(act_stat==1);
    
    cxc=nan(nstat,nstat);
    for i=1:no_act_stat
        d1=squeeze(detrended(:,act_stat(1,i),:));
        for j=i+1:no_act_stat
             d2=squeeze(detrended(:,act_stat(1,j),:));                
             [~,~,r]=canoncorr(d1,d2);
             cxc(act_stat(1,i),act_stat(1,j))=r(1);
             cxc(act_stat(1,j),act_stat(1,i))=r(1);
        end
    end
end

function[Norm_Dist,active_stats]=Normalizer(cxc_mat)
%Calculate what proportion of network would be connected if the correlation
%threshold was t
    active=~isnan(cxc_mat);
    Possible_Connections=sum(active)';
    active_stats=Possible_Connections~=0;
    Norm_Dist=nan(length(cxc_mat),1001);
    for t=0:1/1000:1       
        Norm_Dist(:,round(t*1000+1))=sum(cxc_mat>=t,2)./Possible_Connections;
    end
end


function[struct_overlapping]=Overlapper(struct1,struct2,w)
%Minutes from between two days 
    w1=w-1;
    struct_overlapping.n_nez=nan(struct2.nstats,w1*2);
    struct_overlapping.e_nez=nan(struct2.nstats,w1*2);
    struct_overlapping.z_nez=nan(struct2.nstats,w1*2);
    struct_overlapping.daylit=nan(struct2.nstats,w1*2);
    struct_overlapping.mlat=nan(struct2.nstats,w1*2);
    struct_overlapping.n_nez(:,w:end)=struct2.n_nez(:,1:w1);
    struct_overlapping.e_nez(:,w:end)=struct2.e_nez(:,1:w1);
    struct_overlapping.z_nez(:,w:end)=struct2.z_nez(:,1:w1);
    struct_overlapping.daylit(:,w:end)=struct2.daylit(:,1:w1);
    struct_overlapping.mlat(:,w:end)=struct2.mlat(:,1:w1);

    index=nan(struct1.nstats,struct2.nstats);
        for i=1:struct2.nstats
            index(:,i)=strcmp(struct2.names(i,:),struct1.names);
        end
    [r,c]=find(index==1); 
    w2=w1-1;
    wind=1440-w2:1440;
    struct_overlapping.n_nez(c,1:w1)=struct1.n_nez(r,wind);
    struct_overlapping.e_nez(c,1:w1)=struct1.e_nez(r,wind);
    struct_overlapping.z_nez(c,1:w1)=struct1.z_nez(r,wind);
    struct_overlapping.daylit(c,1:w1)=struct1.daylit(r,wind);
    struct_overlapping.mlat(c,1:w1)=struct1.mlat(r,wind);
    struct_overlapping.nstats=struct2.nstats;
end

function[struct2]=Normalization(struct,num,date,SZA_yr,w)
    [struct2]=ncdf_conver(['supermag/raw/',num2str(year(date)),'/',num2str(year(date)),num2str(month(date),'%02.f'),num2str(day(date),'%02.f'),'.raw-supermag.ncdf']);
    [struct2]=SZA_definer(date,SZA_yr,struct2);
    [struct_overlapping]=Overlapper(struct,struct2,w);
    clear('struct')
    
    Norm_Sum=zeros(struct2.nstats,1001);
    Norm_dark_Sum=Norm_Sum;

    Active_Sum=zeros(struct2.nstats,1);
    Active_Dark_Sum=Active_Sum;

    w1=w-1;
    for t=1:w1
        cxc=canonical_matrix(struct_overlapping,t,w);
        [Norm_Dist,active_stats]=Normalizer(cxc);
        Norm_Dist(isnan(Norm_Dist))=0;
        Norm_Sum=Norm_Sum+Norm_Dist;
        Active_Sum=Active_Sum+active_stats;
        
        Dark_stats=repmat(sum(struct_overlapping.daylit(:,t:t+w1),2)/w,1,struct_overlapping.nstats);
        cxc2=cxc.*Dark_stats.*Dark_stats';
        [Norm_Dist,active_stats]=Normalizer(cxc2);
        Norm_Dist(isnan(Norm_Dist))=0;
        Norm_dark_Sum=Norm_dark_Sum+Norm_Dist;
        Active_Dark_Sum=Active_Dark_Sum+active_stats;               
    end

    for t=1:1440-w1
        cxc=canonical_matrix(struct2,t,w);
        [Norm_Dist,active_stats]=Normalizer(cxc);
        Norm_Dist(isnan(Norm_Dist))=0;
        Norm_Sum=Norm_Sum+Norm_Dist;
        Active_Sum=Active_Sum+active_stats;

        Dark_stats=repmat(sum(struct2.daylit(:,t:t+w1),2)/w,1,struct2.nstats);
        cxc2=cxc.*Dark_stats.*Dark_stats';
        [Norm_Dist,active_stats]=Normalizer(cxc2);
        Norm_Dist(isnan(Norm_Dist))=0;
        Norm_dark_Sum=Norm_dark_Sum+Norm_Dist;
        Active_Dark_Sum=Active_Dark_Sum+active_stats;
    end

    struct2.Norm=Norm_Sum./repmat(Active_Sum,1,1001);
    struct2.Norm_Dark=Norm_dark_Sum./repmat(Active_Dark_Sum,1,1001);
    struct2.Active_Sum=Active_Sum;
    struct2.Active_Sum_Dark=Active_Dark_Sum;
    save(['supermag/processed_w',num2str(w),'/',num2str(year(date)),'/',num2str(num),'.processed.mat'],'struct2','-v7.3')
end

