Mod_flag='normed';
LIM=0.7;
modularity_hist_overlay_med(normalized_times2,networki_v4_w128_t5.eb,intersect([Extreme_quiet],coverage12),'all_12_extreme_quiet_col_w128_t5EB',LIM,Mod_flag)
close all

function[]=modularity_hist_overlay_med(NORM_TM,NETWORK,subs,NAME,probQ,Mod_flag)
close all
clmap =[0         0         0
    0.0345         0    0.0345
    0.0690         0    0.0690
    0.1034         0    0.1034
    0.1379         0    0.1379
    0.1724         0    0.1724
    0.2069         0    0.2069
    0.2414         0    0.2414
    0.2759         0    0.2759
    0.3103         0    0.3103
    0.3448         0    0.3448
    0.3793         0    0.3793
    0.4138         0    0.4138
    0.4483         0    0.4483
    0.4828         0    0.4828
    0.5172         0    0.5172
    0.5517         0    0.5517
    0.5862         0    0.5862
    0.6207         0    0.6207
    0.6552         0    0.6552
    0.6897         0    0.6897
    0.7241         0    0.7241
    0.7586         0    0.7586
    0.7931         0    0.7931
    0.8276         0    0.8276
    0.8621         0    0.8621
    0.8966         0    0.8966
    0.9310         0    0.9310
    0.9655         0    0.9655
    1.0000         0    1.0000
    1.0000    0.0294    1.0000
    1.0000    0.0588    1.0000
    1.0000    0.0882    1.0000
    1.0000    0.1176    1.0000
    1.0000    0.1471    1.0000
    1.0000    0.1765    1.0000
    1.0000    0.2059    1.0000
    1.0000    0.2353    1.0000
    1.0000    0.2647    1.0000
    1.0000    0.2941    1.0000
    1.0000    0.3235    1.0000
    1.0000    0.3529    1.0000
    1.0000    0.3824    1.0000
    1.0000    0.4118    1.0000
    1.0000    0.4412    1.0000
    1.0000    0.4706    1.0000
    1.0000    0.5000    1.0000
    1.0000    0.5294    1.0000
    1.0000    0.5588    1.0000
    1.0000    0.5882    1.0000
    1.0000    0.6176    1.0000
    1.0000    0.6471    1.0000
    1.0000    0.6765    1.0000
    1.0000    0.7059    1.0000
    1.0000    0.7353    1.0000
    1.0000    0.7647    1.0000
    1.0000    0.7941    1.0000
    1.0000    0.8235    1.0000
    1.0000    0.8529    1.0000
    1.0000    0.8824    1.0000
    1.0000    0.9118    1.0000
    1.0000    0.9412    1.0000
    1.0000    0.9706    1.0000
    1.0000    1.0000    1.0000];

A2=0.0075;
A1=(0.8-A2*6)/7;
A3=A1+A2;
sp7=[0.1,0.7,0.8,0.25;
    0.1,0.42,0.8,0.25;
    0.1,0.1,A1,0.21;
    0.1+A3,0.1,A1,0.21;
    0.1+A3*2,0.1,A1,0.21;
    0.1+A3*3,0.1,A1,0.21;
    0.1+A3*4,0.1,A1,0.21;
    0.1+A3*5,0.1,A1,0.21;
    0.1+A3*6,0.1,A1,0.21];
        NT=NORM_TM;  
        for num=1:116
             MOD(num,:)=NETWORK{num,4}(NT(num,:),:);
             NO_CONNECTIONS(num,:)=nansum(NETWORK{num,2}(NORM_TM(num,:),:,1),2);
        end
        NO_CONNECTIONS=(NO_CONNECTIONS<5);

        if Mod_flag=='unnorm'
            MOD_Norm=MOD;
            upperlimm=0.8;
        else
            MOD_Norm=MOD./max(MOD,[],2);
                upperlimm=1;
        end
        MOD_Norm(NO_CONNECTIONS)=nan;
        for i=1:71
            mod_counts2(:,i)=histcounts(MOD_Norm(subs,i),[0:0.05:1]);
        end
        count=1;
        for i=size(mod_counts2,1):-1:1
            mod_counts(count,:)=mod_counts2(i,:);
            count=count+1;
        end
        
    figure('Outerposition',[1 1 1600 1100]);
    subplot('Position',sp7(1,:)); hold on;
    colormap(clmap)
    imagesc([-20:50],[1:-0.05:-0],(mod_counts./sum(mod_counts,1)));
    ylimitz=[0,upperlimm];
    ylim([-0.025,upperlimm+0.025])
    xlim([-20,50])
    a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',4);
    plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',4);
    plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',4);
    a2=plot([30,30],ylimitz,'--','Color',[0.6500,0,0.2],'linewidth',4);
    legend([a1,a2],{'Onset','Peak'},'NumColumns',1,'Position',[0.91,0.6,0.047,0.04])
    xticklabels([]);
    ylabel('Modularity, Q_N','Fontsize',22,'FontWeight','Bold')
    caxis([0,probQ])
    cl=colorbar('Position',[0.91,0.7,0.024,0.25]);
    cl.FontSize=20;
    cl.Label.String='Probability';
    cl.Label.FontSize=22;
    ac=gca;
    ac.Box='on';
    set(gca,'fontsize',22)
    
    subplot('Position',sp7(2,:)); hold on;
    for i=subs
        plot(-21:1:50,MOD_Norm(i,:),'color',[0.7,0.73,0.7,0.15],'linewidth',4.5);
        hold on;
    end
    plot(-21:1:50,median(MOD_Norm(subs,:),'omitnan'),'color',[0,0,0],'linewidth',5)
    plot(-21:1:50,quantile(MOD_Norm(subs,:),0.25),'color',[0.4,0.4,0.4],'linewidth',4)
    plot(-21:1:50,quantile(MOD_Norm(subs,:),0.75),'color',[0.4,0.4,0.4],'linewidth',4)
    ylim([0,upperlimm])
    ylimitz=[0,upperlimm];
    a1=plot([0,0],ylimitz,'--','Color',[0.1000,0.6000,0.2000],'linewidth',4);
    plot([10,10],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',4);
    plot([20,20],ylimitz,'--','Color',[0.5,0.5,0.5],'linewidth',4);
    a2=plot([30,30],ylimitz,'--','Color',[0.6500,0,0.2],'linewidth',4);
    xlabel('Normalized time, t^{\prime} (mins^{\prime})','FontSize',22,'FontWeight','bold')
    ylabel('Modularity, Q_N','Fontsize',22,'FontWeight','Bold')
    ac=gca;
    ac.Box='on';
    xlim([-20,50])
    set(gca,'fontsize',22)

for sp=1:7
    if mod(sp,8)==1
        TM=1:10;
    elseif mod(sp,8)==2
        TM=11:20;
    elseif mod(sp,8)==3
        TM=21:30;
    elseif mod(sp,8)==4
        TM=31:40;
    elseif mod(sp,8)==5
        TM=41:50;
    elseif mod(sp,8)==6
        TM=51:60;
    elseif mod(sp,8)==7
        TM=61:71;
    end
    subplot('Position',sp7(sp+2,:));
    histogram(MOD_Norm(subs,TM),[0:0.05:1],'Normalization','probability')
    hold on;aq=plot([median(median(MOD_Norm(subs,TM),'omitnan'),'omitnan'),median(median(MOD_Norm(subs,TM),'omitnan'),'omitnan')],[0,1],'k','linewidth',3);
    ac=gca;
    ac.Box='on';
    xlim([0,upperlimm])
    xticks([0:0.5:1])
      if sp==1
            ylabel({'Probability'},'FontSize',22,'FontWeight','bold')
            yticks([0:0.1:probQ])
            
            title('-20\leq t^{\prime}\leq -11','FontSize',30)
        elseif sp==2
            yticklabels([])
            title('-10\leq t^{\prime}\leq -1','FontSize',30)
        elseif sp==3
                        yticklabels([])
            title('0\leq t^{\prime}\leq 9','FontSize',30)
        elseif sp==4
                        yticklabels([])
            title('10\leq t^{\prime}\leq 19','FontSize',30)
        elseif sp==5
                        yticklabels([])
            title('20\leq t^{\prime}\leq 29','FontSize',30)
        elseif sp==6
                        yticklabels([])
            title('30\leq t^{\prime}\leq 39','FontSize',30)
        elseif sp==7
                        yticklabels([])
            title('40\leq t^{\prime}\leq 50','FontSize',30)
      end
    xlabel('Q_N','FontSize',22,'FontWeight','bold')
    ylim([0,probQ])
    set(gca,'fontsize',22)
end
legend([aq],{'Median'},'NumColumns',1,'Position',[0.915,0.2,0.047,0.04])
set(gcf,'PaperSize',[55,35])
if Mod_flag=='unnorm'
        print(['/Figures/Communities_Paper_Figures/Mod_op/pdfs/Mod_ophist_unnormd_MED_ge5',NAME],'-dpdf','-r600');
        savefig(['/Figures/Communities_Paper_Figures/Mod_op/figs/Mod_ophist__unnormd_MED_ge5',NAME,'.fig']);
else
        print(['/Figures/Communities_Paper_Figures/Mod_op/pdfs/Mod_ophist_normd_MED_ge5',NAME],'-dpdf','-r600');
        savefig(['/Figures/Communities_Paper_Figures/Mod_op/figs/Mod_ophist__normd_MED_ge5',NAME,'.fig']);
end

end