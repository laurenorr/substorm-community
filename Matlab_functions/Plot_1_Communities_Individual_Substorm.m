clear 
load('normalized_times_2.mat')
load('R/RSubstorms_v4_w128_t5.mat')
load('networki_v4_w128_t5.mat')

%Combine all surrogate files
for i=1:10
%     for num=[1,13,18,29,33]%Example substorms from main text & SI
     for num=[13]%Example substorm from main text
        load(['networki_surrogate_w128_t5_',num2str(i),'.mat'])
        MOD{num,i}=networki_surrogate.eb{num,4};
        alphaN{num,i}=networki_surrogate.all{num,4}(:,1);
    end
end

for num=13%[1,13,18,29,33]
    %Find the average surrogate per substorm
    clear('Mod_ave','alphaB_ave','alphaN_ave')
    for i=1:10
           Mod_ave(i,:)=MOD{num,i};
           alphaN_ave(i,:)=alphaN{num,i};
    end
    Mod_ave=mean(Mod_ave,1,'omitnan');
    alphaN_ave=mean(alphaN_ave,1,'omitnan');
%     Plot figure 1 from main text 
    P=community_plotter(Substorms_v4_w128_t5,normalized_times2,networki_v4_w128_t5.eb,networki_v4_w128_t5.all,'w128_t5_eb',num,Mod_ave,alphaN_ave);
end

%Reshape cell output from panel c&d
 P.CD1=[];
P.CD2=[];
P.CD3=[];
P.CD4=[];
P.CD5=[];
P.CD6=[];
for i=1:92
    P.CD1=[P.CD1,P.CD{i,1}];
    P.CD2=[P.CD2,P.CD{i,2}];
    P.CD3=[P.CD3,P.CD{i,3}];
    P.CD4=[P.CD4,P.CD{i,4}];
    P.CD5=[P.CD5,P.CD{i,5}];
    P.CD6=[P.CD6,P.CD{i,6}];
end

function[P]=community_plotter(SUBSTORMS,NORM_TM,NETWORK,NETWORK_all,SNM,num,MOD_ave,alphaN_ave)
%P output is 
subpositions6=[0.065,0.79,0.8,0.16;
               0.065,0.61,0.8,0.16;
               0.065,0.18,0.8,0.05;
               0.065,0.06,0.8,0.1
               0.065,0.43,0.8,0.16;
               0.065,0.25,0.8,0.16;];

clmp=[ 0.9020    0.1882    0.2000
    0.8680    0.2078    0.2305
    0.8340    0.2274    0.2610
    0.8000    0.2471    0.2915
    0.7660    0.2667    0.3220
    0.6233    0.3490    0.4501
    0.4805    0.4314    0.5782
    0.4171    0.4680    0.6351
    0.3536    0.5046    0.6921
    0.2902    0.5412    0.7490
    0.3033    0.5713    0.6843
    0.3163    0.6013    0.6196
    0.3294    0.6314    0.5549
    0.3425    0.6615    0.4902
    0.3555    0.6915    0.4255
    0.3686    0.7216    0.3608
    0.5555    0.7111    0.2837
    0.7425    0.7007    0.2065
    0.9294    0.6902    0.1294
    0.9529    0.6431    0.1203
    0.9765    0.5961    0.1111
    1.0000    0.5490    0.1020
    0.9941    0.5490    0.2270
    0.9883    0.5490    0.3520
    0.9824    0.5490    0.4770
    0.9765    0.5490    0.6020
    0.9745    0.5490    0.6436
    0.9725    0.5490    0.6853
    0.9706    0.5490    0.7269
    0.9686    0.5490    0.7686];

col_lims=6:12/(length(clmp)-1):18;

col=[0.4    0.4     0.4;
        0         0    0.4510;
        0    0.3333    0.6340;
        0    0.6667    0.8170;
        0    1.0000    1.0000;
        0.3333    1.0000    0.6667;
        0.6667    1.0000    0.3333;
        1.0000    1.0000         0;
        1.0000    0.8333         0;
        1.0000    0.6667         0;
        1.0000    0.5000         0;
        1.0000    0.3333         0;
        1.0000    0.1667         0;
        1.0000         0         0;
        0.8750         0         0;
        0.7500         0         0];

    close all
    fig=figure('Outerposition',[1 1 1200 1600]);
%         ONSETLAT=SUBSTORMS{num,6}.OnsetLat;
%         ONSETMLT=mod(SUBSTORMS{num,6}.OnsetMLT+12,24);
%         T05E=mod(SUBSTORMS{num,6}.T05EMLT+12,24);
%         T05W=mod(SUBSTORMS{num,6}.T05WMLT+12,24);
%         T10E=mod(SUBSTORMS{num,6}.T10EMLT+12,24);
%         T10W=mod(SUBSTORMS{num,6}.T10WMLT+12,24);
%         T15E=mod(SUBSTORMS{num,6}.T15EMLT+12,24);
%         T15W=mod(SUBSTORMS{num,6}.T15WMLT+12,24);
%         T20E=mod(SUBSTORMS{num,6}.T20EMLT+12,24);
%         T20W=mod(SUBSTORMS{num,6}.T20WMLT+12,24);
        NT=NORM_TM(num,1):NORM_TM(num,end);
        MOD=NETWORK{num,4}(NT,:);
        MOD_Norm=MOD;
%         MLTS=mod(squeeze(SUBSTORMS{num,2}(:,1,NT))+12,24);
%         BStats=zeros(size(MLTS));
%         BStats(MLTS<=T10E & MLTS>=T10W)=1;
        NET=NETWORK{num,1};
%         alphaB=zeros(length(NT),1);
%         for tim=1:length(NT)
%            temp_B=BStats(:,tim);
%            temp_B=temp_B.*temp_B';
%            ACTIVEB_struct=SUBSTORMS{num,7}(:,:,NT(1,tim)).*temp_B;
%            AJECB=SUBSTORMS{num,1}(:,:,NT(1,tim)).*ACTIVEB_struct;
%            alphaB(tim,1)=nansum(nansum(AJECB))./nansum(nansum(ACTIVEB_struct));
%         end
        modlims=[0,0.85];
        ACTIVE=SUBSTORMS{num,8}(NT,1);
        timings=SUBSTORMS{num,4}(1,NT);
        Indicies=SUBSTORMS{num,9}(NT+128+18,:);
        alphaN=squeeze(NETWORK_all{num,2}(NT,:,1))./ACTIVE;
        NCom=size(NETWORK{num,2}(NT,:,:),2);
        TM=reshape(repmat(timings,NCom,1),1,[]);
        NEB_mean_mlts=NETWORK{num,2}(NT,:,3);
        NEB_mean_mlats=NETWORK{num,2}(NT,:,7);
        NEB_num_lags=NETWORK{num,3}(NT,:,:);
        NEB_sum_lags=zeros(length(timings),NCom,16);
        NEB_max_mlts=NETWORK{num,2}(NT,:,5);
        NEB_max_mlats=NETWORK{num,2}(NT,:,9);
        NEB_min_mlts=NETWORK{num,2}(NT,:,6);
        NEB_min_mlats=NETWORK{num,2}(NT,:,10);
        for lg=16:-1:1
            NEB_sum_lags(:,:,abs(lg-16-1))=sum(NEB_num_lags(:,:,1:lg),3);
        end
        NEB_sum_lags(NEB_sum_lags==0)=nan;
        
        
        %Panel a
        subplot('Position',subpositions6(1,:))
        colormap(col)
        P.A=[];
        for lg=1:16
            hold on;
            scatter(TM,reshape(NEB_mean_mlts',1,[]),reshape((squeeze(NEB_sum_lags(:,:,lg))./repmat(ACTIVE,1,size(6,2)))',1,[])*14000,col(16+1-lg,:),'filled','MarkerFaceAlpha',1);
            P.A=[P.A;TM;reshape(NEB_mean_mlts',1,[]);reshape((squeeze(NEB_sum_lags(:,:,lg))./repmat(ACTIVE,1,size(6,2)))',1,[])];
        end
        ac=gca;
        ylimitz=[6,18];
        set(gca,'fontsize',22)
        ylabel('$\overline{\theta}$ (MLT)','Interpreter','Latex','Fontsize',29)
        yticks([6,9,12,15,18]);
        yticklabels({'18','21','0','3','6'})
        xticklabels([]);
        title([datestr(SUBSTORMS{num,6}.onset)],'FontSize',25,'FontWeight','bold')
        a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        a2=plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3);
        ylim(ylimitz)
        legend([a1,a2],{'Onset','Peak'},'Position',[0.91,0.92,0,0])
        legend('boxoff')
        
        
        %Panel b
        subplot('Position',subpositions6(2,:));
        colormap(col)
        P.B=[];
        for lg=1:16
            hold on;
            scatter(TM,reshape(NEB_mean_mlats',1,[]),reshape((squeeze(NEB_sum_lags(:,:,lg))./repmat(ACTIVE,1,size(6,2)))',1,[])*14000,col(16+1-lg,:),'filled','MarkerFaceAlpha',1);
            P.B=[P.B;TM;reshape(NEB_mean_mlats',1,[]);reshape((squeeze(NEB_sum_lags(:,:,lg))./repmat(ACTIVE,1,size(6,2)))',1,[])];
        end
        ac=gca;
        caxis([0,15])
        cl=colorbar('Position',[0.88,0.7,0.024,0.18]);
        cl.FontSize=20;
        cl.Label.String='Lag, |\tau_c| (mins)';
        cl.Label.FontSize=25;
        xticklabels([]);
        ylimitz=[60,75];
        set(gca,'fontsize',22)
        ylabel('$\mathbf{\overline{\phi}}$ (MLAT)','Interpreter','Latex','Fontsize',29)
        a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        a2=plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3); %0.6365,0.3753,0.6753
        ylim(ylimitz)

        
        %Panel c & d
        MAXY=reshape(NEB_max_mlts',1,[]);MINY=reshape(NEB_min_mlts',1,[]);MEANY=reshape(NEB_mean_mlts',1,[]); 
        MAXYa=reshape(NEB_max_mlats',1,[]);MINYa=reshape(NEB_min_mlats',1,[]);
        td=(timings(1,4)-timings(1,3))/2;
        for i=1:length(TM)
        if ~isnan(MEANY(1,i))
        tt=TM(1,i);
        qs=clmp(find(col_lims<MEANY(1,i),1,'last'),:);
        subplot('Position',subpositions6(5,:));hold on;
        fill([tt-td,tt+td,tt+td,tt-td],[MINY(1,i),MINY(1,i),MAXY(1,i),MAXY(1,i)],qs,'edgecolor','none');
        subplot('Position',subpositions6(6,:));hold on;
        fill([tt-td,tt+td,tt+td,tt-td],[MINYa(1,i),MINYa(1,i),MAXYa(1,i),MAXYa(1,i)],qs,'edgecolor','none');
        end
        end
        subplot('Position',subpositions6(5,:));
        alpha(0.3)
        subplot('Position',subpositions6(6,:));
        alpha(0.3)
        for i=1:1:length(timings)
            for j=1:size(NET,2)
                if ~isempty(NET{NT(1,i),j})
                    UNIQUE_NET=unique(mod(reshape(NET{NT(1,i),j}(:,[9,11]),1,[])+12,24));
                    meann=mean(UNIQUE_NET,'omitnan');
                    qs=clmp(find(col_lims<meann,1,'last'),:);
                    subplot('Position',subpositions6(5,:));hold on;
                    plot(repmat(timings(1,i),1,length(UNIQUE_NET)),UNIQUE_NET,'s','MarkerEdgeColor',qs,'MarkerFaceColor',qs,'MarkerSize',6)
                    P.CD{i,j}=[repmat(timings(1,i),1,length(UNIQUE_NET));UNIQUE_NET];
                    subplot('Position',subpositions6(6,:));
                    hold on;
                    UNIQUE_NET=unique(reshape(NET{NT(1,i),j}(:,[8,10]),1,[]));
                    plot(repmat(timings(1,i),1,length(UNIQUE_NET)),UNIQUE_NET,'s','MarkerEdgeColor',qs,'MarkerFaceColor',qs,'MarkerSize',6)
                    P.CD{i,j}=[P.CD{i,j};UNIQUE_NET];
                end
            end
        end
        
        %Panel c formating
        subplot('Position',subpositions6(5,:));
        ylim([6,18])
        ylimitz=[6,18];
        a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        a2=plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3);
        set(gca,'fontsize',22)
        ylabel('$\theta$ (MLT)','Interpreter','Latex','Fontsize',29)
        yticks([6,9,12,15,18]);
        yticklabels({'18','21','0','3','6'})
        colormap(clmp)
        caxis([6,18])
        c=colorbar('Position',[0.88,0.31,0.024,0.18]);
        c.FontSize=20;
        c.Label.String='$\overline{\theta}$ (MLT)';
        c.Label.Interpreter='Latex';
        c.Label.FontSize=29;
        c.Ticks=[6,9,12,15,18];
        c.TickLabels={'18','21','0','3','6'};
        xticklabels([]);
        xlim([-20,50])
        
        %Panel d formating
        subplot('Position',subpositions6(6,:));
        ylimitz=[60,75];
        set(gca,'fontsize',22)
        ylabel('$\phi$ (MLAT)','Interpreter','Latex','Fontsize',29)
        a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        a2=plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3);
        ylim(ylimitz)
        xticklabels([]);
        xlim([-20,50])


        %Panel e- Modularity
        subplot('Position',subpositions6(3,:));hold on;
        M2=plot(timings,MOD_ave(1,NT),'k','linewidth',4);
        M1=plot(timings,MOD_Norm,'linewidth',5);
        P.E=[timings;MOD_Norm'];
        ylim(modlims);
        xticklabels([]);
        ylimitz=modlims;
        ylim(ylimitz)
        set(gca,'fontsize',22)
        ylabel('Q','Fontsize',24)
        plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3);
        ylim(ylimitz)
        legend([M1,M2],{'Measured','Surrogate'},'Position',[0.913,0.213,0,0])
        legend('boxoff')
        
        %Panel f- Alpha and SML
        subplot('Position',subpositions6(4,:));hold on;
        plot(timings,alphaN_ave(1,NT),'k','linewidth',4);
        A1=plot(timings,alphaN,'Color',[0.00,0.45,0.74],'linewidth',5);
        ylimitz=[0,1];
        ylim(ylimitz)
        set(gca,'fontsize',22)
        xlabel('Normalized time, t^{\prime} (mins^{\prime})','FontSize',22,'FontWeight','bold')
        plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',3);
        plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',3);
        plot([30,30],ylimitz,'--','Color',[0.6500         0         0.2],'linewidth',3);
        ylim(ylimitz)
        yyaxis right 
        plot(timings,-Indicies(:,2)','color',[0.6350    0.0780    0.1840],'linewidth',5);ac=gca;
        P.F=[timings;alphaN';-Indicies(:,2)'];
        ac.YColor=[0.6350    0.0780    0.1840];
        set(gca,'fontsize',22)
        ylimitz=[0,750];
        ylabel({'SML (nT)'},'Fontsize',25,'FontWeight','Bold')
        ylim(ylimitz)
        yyaxis left
        ylabel('\alpha','Fontsize',30)
        
        %set x limits
        for j=1:6
            subplot('Position',subpositions6(j,:));
            xlim([-20,50])
            ac=gca;
            ac.Box='on';
        end
        
        subplot('Position',subpositions6(2,:))
        ax1=gca;
        colormap(ax1,col);
        subplot('Position',subpositions6(5,:))
        ax5=gca;
        colormap(ax5,clmp);

        annotation('textbox',[0.068,0.970+0.006,0,0],'String','a','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')
        annotation('textbox',[0.068,0.790+0.006,0,0],'String','b','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')
        annotation('textbox',[0.068,0.610+0.006,0,0],'String','c','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')
        annotation('textbox',[0.068,0.430+0.006,0,0],'String','d','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')
        annotation('textbox',[0.068,0.250+0.006,0,0],'String','e','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')
        annotation('textbox',[0.068,0.180+0.006,0,0],'String','f','FitBoxToText','on','FontSize',32,'FontWeight','bold','EdgeColor','none')

        set(gcf,'PaperSize',[43,47])
        print(['pngs/Communities_avsurrogate_unnorm_',SNM,'_',num2str(num)],'-dpng');
end
