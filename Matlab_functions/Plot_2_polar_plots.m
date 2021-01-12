P=polar_plot_printer(networki_v4_w128_t5.eb,'_v4_w128_t5_eb_',Substorms_v4_w128_t5,normalized_times2);
function[P]=polar_plot_printer(cluster_net,TYPE,SUBSTORMS,normalized_times)
close all
%Set colour schemes
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
c6=[0,0,0];
%Begin plot
figure('Outerposition',[1 1 1600 1600]);
for num=[13]
    for tim=16:5:51%1:size(cluster_net{num,1},1)
        clf
        a2=[];
        tm=normalized_times(num,tim);
        activestats=sum(SUBSTORMS{num,7}(:,:,tm))>0;
        q1=[SUBSTORMS{num,10}(:,1,tm),SUBSTORMS{num,10}(:,2,tm)];
        %Mark in black all magnetometers included and their vectors
        polarplot(deg2rad(SUBSTORMS{num,2}(activestats,1,tm)*15),SUBSTORMS{num,2}(activestats,2,tm),'o','color',c6,'MarkerSize',13,'MarkerFaceColor',c6)
        hold on;
        for ssn=1:length(activestats)
            if activestats(1,ssn)==1
                polarplot(deg2rad([SUBSTORMS{num,2}(ssn,1,tm)*15,SUBSTORMS{num,2}(ssn,1,tm)*15+q1(ssn,2)/18]),[SUBSTORMS{num,2}(ssn,2,tm),SUBSTORMS{num,2}(ssn,2,tm)+q1(ssn,1)/18],'color',c6,'linewidth',5)
            end
        end
        Parm_out=[q1,SUBSTORMS{num,2}(:,1,tm),SUBSTORMS{num,2}(:,2,tm)];
        %Group magnetometers into communities with the same colour scheme
        %as figure 1 panels c-d
        for groupo=1:size(cluster_net{num,1},2)
            if ~isempty(cluster_net{num,1}{tm,groupo})
                a1=unique(cluster_net{num,1}{tm,groupo}(:,[1,2]))';
                a2=[a2,reshape(a1,1,[])];
                Parm_out(a1,5)=groupo;
                meann=mean(mod(SUBSTORMS{num,2}(a1,1,tm)+12,24),'omitnan');
                qs=clmp(find(col_lims<meann,1,'last'),:);
                %Plot magnetometers with community colour
                polarplot(deg2rad(SUBSTORMS{num,2}(a1,1,tm)*15),SUBSTORMS{num,2}(a1,2,tm),'o','color',qs,'MarkerSize',13,'MarkerFaceColor',qs)
                LON=deg2rad(SUBSTORMS{num,2}(a1,1,tm)*15); hold on;
                LAT=SUBSTORMS{num,2}(a1,2,tm);
                x=mod(LON+pi/2,2*pi);
                y=LAT;
                DT=delaunayTriangulation(x,y);
                    
                if ~isempty(DT.ConnectivityList) 
                    C=convexHull(DT);
                    %Manual positioning of dashed lines surrounding
                    %communities for times in paper
                    if num==13
                        if tim==21 && groupo==1
                            C=[1;3;4;2;1];
                        elseif tim==31 && groupo==2
                            C=[1;3;4;7;10;5;2;1];
                        elseif tim==31 && groupo==3
                            C=[1;3;5;8;7;6;4;2;1];
                        elseif tim==41 
                            C=[1;8;13;17;23;25;24;22;20;16;11;6;2;1];
                        elseif tim==51
                            C=[1;2;9;14;18;24;26;25;23;22;17;12;7;1];
                        elseif tim==26 && groupo==3
                            C=[1;3;4;2;1];
                         elseif tim==36 && groupo==1   
                            C=[1;6;11;17;16;13;9;4;1];
                         elseif tim==36 && groupo==2
                             C=[1;5;4;2;1];
                         elseif tim==46   
                             C=[1;2;9;18;24;26;25;23;21;17;12;7;1];
                         elseif tim==56    
                             C=[1;2;8;9;12;18;19;17;15;11;6;1];
                        end
                    end
                    if num==1
                        if tim==16 && groupo==1
                            C=[1;2;6;11;12;14;5;1];
                        elseif tim==31 && groupo==1
                            C=[1;3;7;12;17;16;13;6;2;1];
                            elseif tim==36 && groupo==1
                            C=[1;3;4;7;1];
                             elseif tim==51 && groupo==1
                            C=[1;2;6;11;16;17;18;12;5;1];
                        end
                    end
                    if num==29
                        if tim==16 && groupo==1
                        C=[1;2;3;6;10;11;9;7;5;1];  
                        elseif tim==21 && groupo==1
                            C=[1;4;10;13;14;15;11;9;3;1];
                        elseif tim==26 && groupo==1    
                            C=[1;2;6;12;17;19;20;16;13;11;5;1];
                        elseif tim==31 && groupo==1      
                            C=[1;2;6;12;17;19;18;16;13;11;5;1];
                        elseif tim==36 && groupo==1      
                            C=[1;2;6;12;17;19;20;21;16;13;11;5;1];
                        elseif tim==41 && groupo==1      
                            C=[1;3;7;13;17;18;20;21;19;16;14;12;6;1];
                        elseif tim==46 && groupo==1     
                            C=[1;3;7;13;18;20;22;23;24;21;17;14;12;6;1];
                        elseif tim==51 && groupo==1     
                            C=[1;3;7;13;18;20;22;23;21;17;14;12;6;1];
                        end
                    end
                    if num==33
                        if tim==16 && groupo==1 
                            C=[1;2;3;4;5;1];                            
                        elseif tim==21 && groupo==1 
                            C=[1;4;7;10;11;13;12;9;2;1];
                        elseif tim==26 && groupo==1 
                            C=[1;5;7;9;11;8;2;1];
                         elseif tim==31 && groupo==1    
                            C=[1;5;8;9;10;12;7;2;1];
                         elseif tim==36 && groupo==2    
                            C=[1;2;5;11;16;18;17;14;10;4;1];
                          elseif tim==41 && groupo==1     
                            C=[1;5;11;17;19;18;16;14;10;4;2;1];
                         elseif tim==46 && groupo==1     
                            C=[1;2;6;12;18;20;19;17;15;11;5;1];
                         elseif tim==51 && groupo==1       
                            C=[1;3;7;13;19;21;22;16;12;6;1];
                        end
                    end
                   if num==18
                      if tim==26 && groupo==2     
                        C=[1;2;4;6;5;1];
                       elseif tim==36 && groupo==1
                        C=[1;5;9;14;16;15;8;3;1];
                        elseif tim==41 && groupo==2
                        C=[1;2;5;7;6;4;1];
                        elseif tim==51 && groupo==1
                        C=[1;5;8;12;15;16;14;3;1];
                      end
                   end
                    xhull=DT.Points(C,1);
                    yhull=DT.Points(C,2);  
                    thull=mod(xhull-pi/2,2*pi);
                    polarplot(thull,yhull,':','color',qs,'linewidth',6)
                else
                    polarplot([LON(1,1),LON(2,1)],[LAT(1,1),LAT(2,1)],':','color',qs,'linewidth',6)
                end
                %Rotate map to correct orientation
                ax=gca;
                ax.ThetaTick=[0:45:360];ax.ThetaTickLabel={'0','3','6','9','12','15','18','21'};
                ax.RLim=[58,90];
                ax.RDir='reverse';
                ax.ThetaZeroLocation='bottom';
                ax.Color='none';
                ax.RAxisLocation=90;
                
                %Plot magnetometer vectors with community colour
                Q1=q1(a1,:);
                LON=SUBSTORMS{num,2}(a1,1,tm);
                LONG=deg2rad([LON*15,LON*15+Q1(:,2)/18]);
                LATI=[LAT,LAT+Q1(:,1)/18];
                hold on;
                for ssn=1:length(a1)
                        polarplot(LONG(ssn,:),LATI(ssn,:),'color',qs,'linewidth',5)
                end
            end
        end
        Parm_out(activestats==0,:)=[];
        title(datestr(SUBSTORMS{num,5}(1,tm)),'fontweight','bold')
        set(gca,'fontsize',30)
        set(gcf,'PaperSize',[43,47])
        print(['pngs/polar_plot_',TYPE,'_',num2str(num),'_',num2str(tim)],'-dpng');
        P{tim,1}=Parm_out;
    end
end
end