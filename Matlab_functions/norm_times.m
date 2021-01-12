%<norm_times.m Calculate normalized times.>
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
load(['R/RSubstorms_v4_w128_t5.mat'])
for num=1:116
    times2=[find( Substorms_v4_w128_t5{num,4}<=-21,1,'last'):1:find( Substorms_v4_w128_t5{num,4}>=50,1,'first')];
    timings2=Substorms_v4_w128_t5{num,4}(1,times2);
    for tm=-20:1:50
        normalized_times2(num,tm+21)=times2(1,find(timings2>=tm,1,'first'));
    end
end
save('normalized_times_2.mat','normalized_times2','-v6')
