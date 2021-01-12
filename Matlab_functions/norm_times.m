load(['R/RSubstorms_v4_w128_t5.mat'])
for num=1:116
    times2=[find( Substorms_v4_w128_t5{num,4}<=-21,1,'last'):1:find( Substorms_v4_w128_t5{num,4}>=50,1,'first')];
    timings2=Substorms_v4_w128_t5{num,4}(1,times2);
    for tm=-20:1:50
        normalized_times2(num,tm+21)=times2(1,find(timings2>=tm,1,'first'));
    end
end
save('normalized_times_2.mat','normalized_times2','-v6')