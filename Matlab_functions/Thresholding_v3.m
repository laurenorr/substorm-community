%<Thresholding_v3.m Function thresholds the correlation matrix from substorm structure using a
%    months worth of data around the event.>
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

for storm_num=[1:116]
    threshold_substorms(storm_num,'OFF')
end

for storm_num=[1,13,18,29,33]
    for surr=1:10
    threshold_substorms(storm_num,'ON ',surr)
    end
end
        
function[thresholds_Dark]=threshold_substorms(storm_num,SURR,surrno)
        w=128;
        nm=1001;
        if SURR=='OFF'
            load(['Substorms_v3/',num2str(storm_num),'.substorm_v3.mat'])
        else
            load(['Substorms_v3_surrogates/',num2str(storm_num),'.substorm_v3',num2str(surrno),'.mat'])
        end
        Norm=zeros(Substorm.nstats,nm);
        Norm_Dark=Norm;
        Active_Sum=zeros(Substorm.nstats,1);
        Active_Sum_Dark=Active_Sum;
        w1=w-1;
        start_date=Substorm.info.onset;
        surrounding_month=start_date-days(14):days(1):start_date+days(13);
        for i=surrounding_month
            load(['supermag/processed_w',num2str(w),'/',num2str(year(i)),'/',num2str(year(i)),num2str(month(i),'%02.f'),num2str(day(i),'%02.f'),'.processed.mat'])
            index=nan(struct2.nstats,Substorm.nstats);
            for j=1:Substorm.nstats
                index(:,j)=strcmp(Substorm.names(j,:),struct2.names);
            end
            [r,c]=find(index==1);   
            dum_normdark=struct2.Norm_Dark(r,:);
            dum_normdark(isnan(dum_normdark))=0;

            Norm_Dark(c,:)=Norm_Dark(c,:)+dum_normdark.*repmat(struct2.Active_Sum_Dark(r,:),1,nm);
            Active_Sum(c,:)=Active_Sum(c,:)+struct2.Active_Sum(r,:);
            Active_Sum_Dark(c,:)=Active_Sum_Dark(c,:)+struct2.Active_Sum_Dark(r,:);
        end
        
        Norm_Dark=Norm_Dark./repmat(Active_Sum_Dark,1,nm);
        thresholds_Dark=zeros(Substorm.nstats,10);
%Find the thresholds corresponding to a station being connected to 1-10% of the network throughout the surrounding month 
        for k=1:Substorm.nstats
            for l=1:10
                 thresholds_Dark(k,l)=(dsearchn(Norm_Dark(k,:)',l/100))/1001;         
            end
        end
%Exclude stations which were inactive from more than 30% of the month
        Possible_stats=28*1440;
        Stat_Dark_Tresholds=zeros(Substorm.nstats,Substorm.nstats,10);
        active=Active_Sum;
        active(active<0.7.*Possible_stats)=nan;
        active(active>=0.7.*Possible_stats)=1;
        active=repmat(active.*active',1,1,Substorm.length,10);
%Remove A_{ii} (diagonal) 
        A=repmat(squeeze(Substorm.CanonXCorr(:,:,19,:)),1,1,1,10);
        A(~isnan(A))=1;
        fixer=repmat(ones(Substorm.nstats)-eye(Substorm.nstats),1,1,Substorm.length,10);
        active=active.*A.*fixer;
%Calculate the minimum threshold per magnetometer pair.
        for i=1:1:Substorm.nstats   
           for j=1+i:1:Substorm.nstats
                    Stat_Dark_Tresholds(i,j,:)=min(thresholds_Dark(i,:),thresholds_Dark(j,:));
                    Stat_Dark_Tresholds(j,i,:)=Stat_Dark_Tresholds(i,j,:);

           end
        end

        thresholds_Dark=permute(repmat(Stat_Dark_Tresholds,1,1,1,Substorm.length),[1,2,4,3]);

        %Exclude stations that only have low values of magnetic field perturbations (<20)
        noisey=ones(Substorm.nstats,Substorm.nstats,Substorm.length);
        for i=19:Substorm.length+18
            t=i-18:i+w1+18;
            for j=1:Substorm.nstats
                dum=[Substorm.n_nez(j,t);Substorm.e_nez(j,t);Substorm.z_nez(j,t)];
                if nansum(nansum(abs(dum)>20))==0
                    noisey(j,:,i-18)=nan;
                    noisey(:,j,i-18)=nan;
                end
            end
        end
        Substorm.noise=noisey;
%       Zero lagged network
        conXcorr=repmat(squeeze(Substorm.CanonXCorr(:,:,19,:)),1,1,1,10); % nstat x nstat x timeframe x thresholds matrix of canonical cross correlation
        Substorm.std.CXC_Dark=(conXcorr>=thresholds_Dark).*active.*repmat(Substorm.noise,1,1,1,10);
%       Lag of maximum correlation network
        MAX_CXC=repmat(Substorm.Max_CXC,1,1,1,10);
        Substorm.std.MaxCXC_Dark=(MAX_CXC>=thresholds_Dark).*active.*repmat(Substorm.noise,1,1,1,10);
        if SURR=='OFF'
            save(['Substorms_v3/',num2str(storm_num),'.substorm_v3.mat'],'Substorm','-v7.3')
        else
            save(['Substorms_v3_surrogates/',num2str(storm_num),'.substorm_v3',num2str(surrno),'.mat'],'Substorm','-v7.3')
        end
end
