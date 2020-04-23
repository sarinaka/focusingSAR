%Nicole Bienert, Sean Peters
%Purpose: plot aPRES measurements

%Instructions: Change the variables file, burst, datafolder, figureFolder.
%If you want the plots to save automatically, uncomment the saveas

clc
close all
clear all


%% vars
%The pRES data file. Do not include the .dat file extension
file = 'Survey_2018-07-05_133859_nearMetal' 
%The folder where the data file exists
dataFolder = 'C:\Users\jeany\OneDrive - Leland Stanford Junior University\Documents\School\Research\Greenland 2018\measurements\point measurements\7_05_DrillSite\Data';
%If plots are saving automatically, where do you want them saved to?
figureFolder = 'C:\Users\jeany\Downloads';
%include the pRES processing codes (fmcw_Catherine) in matlab's path
addpath(genpath(pwd))


fileType='.dat'; %file extension
maxH=2000; %Limit the max ice thickness on the plot
permittivity = 3.18; %relative permittivity of ice


filename= [dataFolder,'\',file,fileType]

% calculations
%read the header and extract the voltage data from the PRES file
data = fmcw_load(filename,permittivity);
%calculate range data
[Rcoarse,Rfine,spec_cor,spec] =fmcw_range(data,1,maxH,@blackman);
%average all bursts
R = abs(mean(spec_cor));
%find the max amplitude
maxA = max(max(abs(spec_cor(:,1000:4000))));

%Plot
figure()
plot(Rcoarse,abs(spec_cor)),xlim([0 maxH]),ylim([0 1.5*maxA])
hTitle = title(regexprep(file,"_"," "))
hXlabel = xlabel('Range')
hYlabel = ylabel('|Magnitude|')
Aesthetics_Script
% saveas(gcf, [figureFolder,'\',file,'-ptMeasurement'], 'fig')
% saveas(gcf, [figureFolder,'\',file,'-ptMeasurement'], 'png')
    
figure()
plot(Rcoarse,10*log10(abs(R)),'b-'),xlim([0 maxH])
hTitle = title(regexprep(file,"_"," "))
hXlabel=xlabel('Range')
hYlabel = ylabel('|Magnitude| (dB)')
Aesthetics_Script
% saveas(gcf, [figureFolder,'\',file,'-ptMeasurement-coherentSum'], 'fig')
% saveas(gcf, [figureFolder,'\',file,'-ptMeasurement-coherentSum'], 'png')

