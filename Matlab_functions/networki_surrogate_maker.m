%Load R files and create network files. 

%Total network file
for surr=1:10
%Total network file
load(['R/RSubstorms_surrogate_w128_t5_',num2str(surr),'.mat'])
[networks3_surrogate]=NETSFIND(Substorms_surrogate_w128_t5);
save(['networks3_surrogate_w128_t5_',num2str(surr),'.mat'],'networks3_surrogate','-v6')

%Network with community information
[networki_surrogate]=unspec_no_comms(networks3_surrogate,['w128_t5_',num2str(surr)]);
[networki_surrogate]=bitgiverer(networki_surrogate);
networki_surrogate=modularity_putter(networki_surrogate,['w128_t5_',num2str(surr)]);
networki_surrogate.all=alpha_calculater(Substorms_surrogate_w128_t5,networki_surrogate.all);
save(['networki_surrogate_w128_t5_',num2str(surr),'.mat'],'networki_surrogate','-v7.3')
end

function[NETWORK_all]=alpha_calculater(SUBSTORMS,NETWORK_all)
%Calculate the normalized number of connections
    for num=[1,13,18,29,33]
        NT=1:length(SUBSTORMS{num,4});
        clear('alphaA')
        for tim=NT
           ACTIVE_struct=SUBSTORMS{num,7}(:,:,NT(1,tim));
           for i=1:size(ACTIVE_struct,1)
                ACTIVE_struct(i,i,:)=0;
           end
           AJEC=SUBSTORMS{num,1}(:,:,NT(1,tim)).*ACTIVE_struct;
           alphaA(tim,1)=nansum(nansum(AJEC))./nansum(nansum(ACTIVE_struct));
        end
        NETWORK_all{num,4}=[alphaA];
    end
end


function[netowrki]=modularity_putter(netowrki,Threshold)
%Load R modularity files
    for num=[1,13,18,29,33]
              load(['R/',num2str(num),'modularity_eb_surrogate_',num2str(Threshold),'.mat'])
              netowrki.eb{num,4}=modularity_eb;
    end
end

function[networks3]=NETSFIND(Substorms)
%Calculate total network from upper and lower sections
    for num=[1,13,18,29,33]
        [networks{num,1},networks2{num,1}]=networkfinder2(Substorms,num);
    end
    for i=[1,13,18,29,33]
        for j=1:size(networks{i,1},1)
            networks3{i,1}{j,1}=[networks{i,1}{j,1};networks2{i,1}{j,1}];
        end
        for j=1:size(networks{i,1},1)
            for k=1:size(networks{i,1}{j,1},1)
                if networks{i,1}{j,1}(k,7)<0
                    networks{i,1}{j,1}(k,:)=networks{i,1}{j,1}(k,[2,1,5,6,3,4,7,10,11,8,9]);
                    networks{i,1}{j,1}(k,7)=-networks{i,1}{j,1}(k,7);
                end
            end
            for k=1:size(networks2{i,1}{j,1},1)
                if networks2{i,1}{j,1}(k,7)<0
                    networks2{i,1}{j,1}(k,:)=networks2{i,1}{j,1}(k,[2,1,5,6,3,4,7,10,11,8,9]);
                    networks2{i,1}{j,1}(k,7)=-networks2{i,1}{j,1}(k,7);
                end
            end
            for k=1:size(networks3{i,1}{j,1},1)
                if networks3{i,1}{j,1}(k,7)<0
                    networks3{i,1}{j,1}(k,:)=networks3{i,1}{j,1}(k,[2,1,5,6,3,4,7,10,11,8,9]);
                    networks3{i,1}{j,1}(k,7)=-networks3{i,1}{j,1}(k,7);
                end
            end
        end
    end
end

function[networki]=bitgiverer(networki)
    networki.eb=bitgiver(networki.eb);
    networki.all=bitgiver(networki.all);
end

function[neti]=bitgiver(neti)
%Add useful parameter information to network communities 
    for num=[1,13,18,29,33]
        for i=1:size(neti{num,1},1)
            for j=1:size(neti{num,1},2)
                if ~isempty(neti{num,1}{i,j})
                   MLTs=mod(unique(neti{num,1}{i,j}(:,[9,11]))+12,24);
                   MLATs=unique(neti{num,1}{i,j}(:,[8,10]));
                   neti{num,2}(i,j,1)=size(neti{num,1}{i,j},1); %number of connections
                   neti{num,2}(i,j,2)=length(unique(neti{num,1}{i,j}(:,1:2))); %number of stations in community
                   neti{num,2}(i,j,3)=mean(MLTs,'omitnan'); %mean rotated mlt (mod(mlt+12,24))
                   neti{num,2}(i,j,4)=std(MLTs,'omitnan'); %std rotated mlt 
                   neti{num,2}(i,j,5)=max(MLTs); %max rotated mlt 
                   neti{num,2}(i,j,6)=min(MLTs); %min rotated mlt 
                   neti{num,2}(i,j,7)=mean(MLATs,'omitnan'); %mean mlat
                   neti{num,2}(i,j,8)=std(MLATs,'omitnan'); %std mlat 
                   neti{num,2}(i,j,9)=max(MLATs); %max mlat 
                   neti{num,2}(i,j,10)=min(MLATs); %min mlat 
                    [a,ia,ic]=unique(sort(neti{num,1}{i,j}(:,1:2),2),'rows');
                    LAGS=neti{num,1}{i,j}(ia,7);
                    neti{num,2}(i,j,11)=mean(LAGS,'omitnan'); %mean lags
                    neti{num,2}(i,j,12)=std(LAGS,'omitnan'); %std lags 
                    neti{num,2}(i,j,13)=max(LAGS); %max lags
                    neti{num,2}(i,j,14)=min(LAGS); %min lagc
                    for lg=0:15
                        neti{num,3}(i,j,lg+1)=sum(neti{num,1}{i,j}(:,7)==lg);
                    end
                else
                    neti{num,2}(i,j,1:14)=nan(1,1);
                    neti{num,3}(i,j,1:16)=nan(1,1);
                end
            end
        end
    end
end

function[networks,networks2]=networkfinder2(Substorms_v3,num)
%Calculate network from substorm ajecency matrix
   Sub=Substorms_v3{num,1};
   Subdets=Substorms_v3{num,2};
   glat=squeeze(Subdets(:,4,:));
   glon=squeeze(Subdets(:,3,:));
   mlat=squeeze(Subdets(:,2,:));
   mlt=squeeze(Subdets(:,1,:));
   Max_Lags=Substorms_v3{num,3};
    AJEC=squeeze(Sub);
    w1=127+18;
    durin=size(AJEC,3);
        networks{durin,1}=[];
        networks2{durin,1}=[];
        for t=1:durin
            tnetwork=zeros(nansum(nansum(AJEC(:,:,t))),11);
            tnetwork2=zeros(nansum(nansum(AJEC(:,:,t))),11);
            counter=1;
            counter2=1;
            for i=1:size(AJEC,2)
                for j=i+1:size(AJEC,2)
                    %Upper network
                    if AJEC(i,j,t)==1 
                            tnetwork(counter,:)=[i,j,glat(i,t),glon(i,t),glat(j,t),glon(j,t),Max_Lags(i,j,t),mlat(i,t),mlt(i,t),mlat(j,t),mlt(j,t)];
                            counter=counter+1;
                    end
                    %Lower network 
                    if AJEC(j,i,t)==1 
                            tnetwork2(counter2,:)=[j,i,glat(j,t),glon(j,t),glat(i,t),glon(i,t),Max_Lags(j,i,t),mlat(j,t),mlt(j,t),mlat(i,t),mlt(i,t)];
                            counter2=counter2+1;
                    end

                end
            end
        networks{t,1}=tnetwork(1:counter-1,:);
        networks2{t,1}=tnetwork2(1:counter2-1,:);
        end
end

function[netw]=unspec_no_comms(Networks,Threshold)
%Load community R files and save unique communities as part of the network
for num=[1,13,18,29,33]
    load(['R/',num2str(num),'communities_eb_surrogate_',num2str(Threshold),'.mat'])
    net=Networks{num,1}; 
    duration=size(net,1);
    netw.all=Networks;
    for tm=1:duration
        num_eb_clusts=unique(cluster_list_eb(tm,:));
        net_eb{num,1}{duration,1}=[];

        counter=1;
        for n=num_eb_clusts
            [a1,~]=find((net{tm,1}(:,1)==find(cluster_list_eb(tm,:)==n)));
            [a2,~]=find((net{tm,1}(:,2)==find(cluster_list_eb(tm,:)==n)));
            a3=intersect(a1,a2);
            if ~isempty(a3)
                net_eb{num,1}{tm,counter}=net{tm,1}(a3,:);
                counter=1+counter;
            end
        end

    end
    
end
netw.eb=net_eb;
end