%<SZA_sorter.m Read in a year of SZA data txt downloaded from superMAG.>
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
%    Read in a year of SZA data txt downloaded from superMAG. Write as a structure
%    with data matrix- no. stations x minutes in year
%    and names vector- no.stations
Q=dlmread('20180116-09-34-supermag.txt','\t',73,5);
fileID=fopen('20180116-09-34-supermag.txt');
names=textscan(fileID,'%s %*[^\n]');
Names=names{1,1}(54:end,:);
SZA=Q(:,1);
ordered_names=unique(Names);
no_names=length(ordered_names);
names_indexed=nan(no_names-1,1440*366);
for i=2:no_names
    indexes=strcmp(ordered_names(i,1),Names);
    r=find(indexes==1);
    names_indexed(i-1,:)=r';
end

SZA_ordered=nan(no_names-1,1440*366);
for q=2:no_names
    for w=1:1440*366
        SZA_ordered(q-1,w)=SZA(names_indexed(q-1,w),1);
    end
end
SZA2000.data=SZA_ordered;
SZA2000.names=ordered_names(2:end,1);

save('SZA_data/SZA2000.mat','SZA2000','-v7.3')
