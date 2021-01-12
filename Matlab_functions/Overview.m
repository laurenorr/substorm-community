%<Overview.m Running instructions for substorm community code.>
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

%Make Communities_paper the main directory and add all folders within it to
%the path

%Preprocessing of data per day- the most recent data should be downloaded
%from superMAG and functions may need editing depending on format. Processing takes approx 10 minutes per day and >1000 days
%were used.
% .mat files used in orr.et al are included in code folder.

 %Read in a year of SZA data txt downloaded from superMAG. Write as a structure
 %with data matrix- no. stations x minutes in year
 %and names vector- no.stations
 %Code will need editing depending on format downloaded and what year is being used.
 SZA_sorter.m
 %Code is provided as an example of the shape of the structures required for later code.
% 
%Code to read ncdf files into a matlab structure. Data is saved as year month day in separate year folders between 1996 to
%2001. Each day contains a normalization matrix of what proportion of the network would be connected if the correlation
 %threshold was t where t=0:0.001:1. This matrix will later be used for a
%station and event specific threshold.
%Should be edited to reflected most recent version of SuperMAG dataset.
 %Takes approximated 10 minutes per day
 ncdf_converter


% ----------------------------------------------------------------------------------------------------------------------------------

%Processing data relevent for 116 substorms by loading individual days. Includes zero lag cross-correlation matrix (CCM),
% max lag CCM and a matrix of the corresponding maximum lags
%Takes approx 10 minutes per substorm
Substorm_Lags_pm15

%Function thresholds the correlation matrix from substorm structure using a
%months worth of data around the event.
%Takes approx 20 minutes total
Thresholding_v3

%Saves the substorm structure in an R compatible manner
%5 minutes
Substorm_extracter

% ----------------------------------------------------------------------------------------------------------------------------------
% Run R functions, ~2 minutes. igraph package is used.
% Run in R, Membership_save.R
% Run in R, Membership_save_surrogates.R
% ----------------------------------------------------------------------------------------------------------------------------------

%Load R files and make network with community information, 2 minutes
networki_function_maker
networki_surrogate_maker

%Calculate normalized times, 5 seconds
norm_times

% ----------------------------------------------------------------------------------------------------------------------------------

%Plots from main text
Plot_1_Communities_Individual_Substorm
Plot_2_polar_plots
Plot_3_Modularity_histogram
