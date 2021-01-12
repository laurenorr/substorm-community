%<Substorm_extracter.m Saves the substorm structure in an R compatible manner.>
%    Copyright (C) <2021>  <Lauren Orr>

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%    Processing data relevent for 116 substorms by loading individual days. Includes zero lag cross-correlation matrix (CCM),
%    max lag CCM and a matrix of the corresponding maximum lags
 [Substorms_v4_w128_t5]=SE(5,127);
 save('R/RSubstorms_v4_w128_t5.mat','Substorms_v4_w128_t5','-v6')

for surr=1:10
    [Substorms_surrogate_w128_t5]=SE_surrogate(5,127,surr);
    save(['R/RSubstorms_surrogate_w128_t5_',num2str(surr),'.mat'],'Substorms_surrogate_w128_t5','-v6')
end

function[Substorms_v4]=SE(THRS,nic)
for num=1:116
    load(['Substorms_v3/',num2str(num),'.substorm_v3.mat'])
    %Save info in cells format
    w1=18+nic;
    Substorms_v4{num,1}=zeros(Substorm.nstats,Substorm.nstats,Substorm.length);
    Substorms_v4{num,2}=zeros(Substorm.nstats,4,Substorm.length);
    Substorms_v4{num,3}=zeros(Substorm.nstats,Substorm.nstats,Substorm.length);
    Substorms_v4{num,4}=zeros(1,Substorm.length);
    Substorms_v4{num,5}=zeros(1,Substorm.length);
    for t=1+w1:1:Substorm.length+w1
        mlt=Substorm.mlt(:,t);
        [~,indexmlt]=sort(mod(mlt+12,24));
        Substorms_v4{num,13}(:,:,t-w1)=Substorm.Geodetic_dist(indexmlt,indexmlt);
        Substorms_v4{num,1}(:,:,t-w1)=Substorm.std.MaxCXC_Dark(indexmlt,indexmlt,t-w1,THRS);
        Substorms_v4{num,2}(:,:,t-w1)=[Substorm.mlt(indexmlt,t),Substorm.mlat(indexmlt,t),Substorm.glon(indexmlt,t),Substorm.glat(indexmlt,t)];
        Substorms_v4{num,3}(:,:,t-w1)=Substorm.Max_Lags(indexmlt,indexmlt,t-w1);
        Substorms_v4{num,10}(:,:,t-w1)=[Substorm.n_nez(indexmlt,t),Substorm.e_nez(indexmlt,t),Substorm.z_nez(indexmlt,t)];
        Substorms_v4{num,11}(:,:,t-w1)=[mean(Substorm.n_nez(indexmlt,t-nic:t),2,'omitnan'),mean(Substorm.e_nez(indexmlt,t-nic:t),2,'omitnan'),mean(Substorm.z_nez(indexmlt,t-nic:t),2,'omitnan')];
        Substorms_v4{num,12}(:,t-w1)=Substorm.names(indexmlt,1);
        q=squeeze(Substorms_v4{num,2}(:,1,t-w1));
        %  Exlude magnetometers not in the nightside (18-6 mlt)
        w=mod(q,18)>6;
        if sum(w)>0
            Substorms_v4{num,1}(:,w,t-w1)=nan;
            Substorms_v4{num,1}(w,:,t-w1)=nan;
            Substorms_v4{num,2}(w,:,t-w1)=nan;
        end
        q=squeeze(Substorms_v4{num,2}(:,2,t-w1));
        %Exclude magnetometers below 60 degrees and above 75 degrees
        w=find(q<60 | q>75);
        if sum(w)>0
            Substorms_v4{num,1}(:,w,t-w1)=nan;
            Substorms_v4{num,1}(w,:,t-w1)=nan;
            Substorms_v4{num,2}(w,:,t-w1)=nan;
        end
         %Use Separation matrix to exclude magnetometers which are too
         %close together (2*10^5 metres)
        A=Substorms_v4{num,13}(:,:,t-w1)<(2*10^5);
        A=A-eye(Substorm.nstats);
        listy=[];
        while nansum(nansum(A))>1
            [a1,~]=find(A==1);
            [c1,c2]=histcounts(a1,[min(a1):1:max(a1)+1]);
            [~,ind]=max(c1);
            listy=[listy,c2(1,ind)];
            A(c2(1,ind),:)=nan;
            A(:,c2(1,ind))=nan;
        end
        Substorms_v4{num,1}(listy,:,t-w1)=nan;
        Substorms_v4{num,1}(:,listy,t-w1)=nan;
        Substorms_v4{num,2}(listy,:,t-w1)=nan;
        Substorms_v4{num,3}(listy,:,t-w1)=nan;
        Substorms_v4{num,3}(:,listy,t-w1)=nan;
        Substorms_v4{num,13}(listy,:,t-w1)=nan;
        Substorms_v4{num,13}(:,listy,t-w1)=nan;
    end

            Substorms_v4{num,4}=[Substorm.timings];
            Substorms_v4{num,5}=[Substorm.datetimes];
            Substorms_v4{num,6}=[Substorm.info];
            Substorms_v4{num,7}=~isnan(Substorms_v4{num,1});
            Active_stats=Substorms_v4{num,7};
            for i=1:size(Active_stats,1)
                Active_stats(i,i,:)=0;
            end
            Substorms_v4{num,8}=squeeze(sum(sum(Active_stats,2),1));
            Substorms_v4{num,9}=Substorm.indicies;
end
end

function[Substorms_v4]=SE_surrogate(THRS,nic,surrno)
for num=[1,13,18,29,33]
    load(['Substorms_v3_surrogates/',num2str(num),'.substorm_v3',num2str(surrno),'.mat'])
    %Save info in cells format
    w1=18+nic;
    Substorms_v4{num,1}=zeros(Substorm.nstats,Substorm.nstats,Substorm.length);
    Substorms_v4{num,2}=zeros(Substorm.nstats,4,Substorm.length);
    Substorms_v4{num,3}=zeros(Substorm.nstats,Substorm.nstats,Substorm.length);
    Substorms_v4{num,4}=zeros(1,Substorm.length);
    Substorms_v4{num,5}=zeros(1,Substorm.length);
    for t=1+w1:1:Substorm.length+w1
        mlt=Substorm.mlt(:,t);
        [~,indexmlt]=sort(mod(mlt+12,24));
        Substorms_v4{num,13}(:,:,t-w1)=Substorm.Geodetic_dist(indexmlt,indexmlt);
        Substorms_v4{num,1}(:,:,t-w1)=Substorm.std.MaxCXC_Dark(indexmlt,indexmlt,t-w1,THRS);
        Substorms_v4{num,2}(:,:,t-w1)=[Substorm.mlt(indexmlt,t),Substorm.mlat(indexmlt,t),Substorm.glon(indexmlt,t),Substorm.glat(indexmlt,t)];
        Substorms_v4{num,3}(:,:,t-w1)=Substorm.Max_Lags(indexmlt,indexmlt,t-w1);
        Substorms_v4{num,10}(:,:,t-w1)=[Substorm.n_nez(indexmlt,t),Substorm.e_nez(indexmlt,t),Substorm.z_nez(indexmlt,t)];
        Substorms_v4{num,11}(:,:,t-w1)=[mean(Substorm.n_nez(indexmlt,t-nic:t),2,'omitnan'),mean(Substorm.e_nez(indexmlt,t-nic:t),2,'omitnan'),mean(Substorm.z_nez(indexmlt,t-nic:t),2,'omitnan')];
        Substorms_v4{num,12}(:,t-w1)=Substorm.names(indexmlt,1);
        q=squeeze(Substorms_v4{num,2}(:,1,t-w1));
        %  Exlude magnetometers not in the nightside (18-6 mlt)
        w=mod(q,18)>6;
        if sum(w)>0
            Substorms_v4{num,1}(:,w,t-w1)=nan;
            Substorms_v4{num,1}(w,:,t-w1)=nan;
            Substorms_v4{num,2}(w,:,t-w1)=nan;
        end
        q=squeeze(Substorms_v4{num,2}(:,2,t-w1));
        %Exclude magnetometers below 60 degrees and above 75 degrees
        w=find(q<60 | q>75);
        if sum(w)>0
            Substorms_v4{num,1}(:,w,t-w1)=nan;
            Substorms_v4{num,1}(w,:,t-w1)=nan;
            Substorms_v4{num,2}(w,:,t-w1)=nan;
        end
         %Use Separation matrix to exclude magnetometers which are too
         %close together (2*10^5 metres)
        A=Substorms_v4{num,13}(:,:,t-w1)<(2*10^5);
        A=A-eye(Substorm.nstats);
        listy=[];
        while nansum(nansum(A))>1
            [a1,~]=find(A==1);
            [c1,c2]=histcounts(a1,[min(a1):1:max(a1)+1]);
            [~,ind]=max(c1);
            listy=[listy,c2(1,ind)];
            A(c2(1,ind),:)=nan;
            A(:,c2(1,ind))=nan;
        end
        Substorms_v4{num,1}(listy,:,t-w1)=nan;
        Substorms_v4{num,1}(:,listy,t-w1)=nan;
        Substorms_v4{num,2}(listy,:,t-w1)=nan;
        Substorms_v4{num,3}(listy,:,t-w1)=nan;
        Substorms_v4{num,3}(:,listy,t-w1)=nan;
        Substorms_v4{num,13}(listy,:,t-w1)=nan;
        Substorms_v4{num,13}(:,listy,t-w1)=nan;
    end

            Substorms_v4{num,4}=[Substorm.timings];
            Substorms_v4{num,5}=[Substorm.datetimes];
            Substorms_v4{num,6}=[Substorm.info];
            Substorms_v4{num,7}=~isnan(Substorms_v4{num,1});
            Active_stats=Substorms_v4{num,7};
            for i=1:size(Active_stats,1)
                Active_stats(i,i,:)=0;
            end
            Substorms_v4{num,8}=squeeze(sum(sum(Active_stats,2),1));
            Substorms_v4{num,9}=Substorm.indicies;
end
end
