function []=PolarPlotdB(theta,data,range,Ntick,Line,option)
%   function [] = PolarPlotdBtheta,data,range,Ntick,Line,option)
%   Polar_dB plots the data in dB 'data' vs 'theta' is polar fashion. 
%
%  	theta 	: in radian
%   data   	: in dB
%   range   : e.g. [ min max ], default [-40 0] 
%   Ntick   : No. of radial ticks, default 5
% 	Line    : Line properties e.g. '--b', defaul 'k'
% 	option  : 1) Full circle 2) Half circle display, defaul 1 
%%
switch nargin
    case 2
        range =	[-40 round(max(data))]; % default
        Ntick =	5;
		Line = 'k';
		option = 1;
	case 3
        Ntick =	5;
		Line = 'k';
		option = 1;
	case 4
		Line = 'k';
		option = 1;
    case 5
		option = 1;
end
data(isnan(data))       =   min(range);
data(data < min(range))	=   min(range);
polarplot(theta,data,Line,'LineWidth',1)
rlim(range)
%% 
if option == 1
thetaticks(0:15:345)
thetaticklabels({'0^\circ',' ','30^\circ',' ','60^\circ',' ','90^\circ',' ','120^\circ',' ','150^\circ',' ','180^\circ',' ','-150^\circ',' ','-120^\circ',' ','-90^\circ',' ','-60^\circ',' ','-30^\circ'})
end
%% 
if option == 2
thetalim([-90 90])
thetaticks(-90:15:90)
thetaticklabels({'-90^\circ',' ','-60^\circ',' ','-30^\circ',' ','0^\circ',' ','30^\circ',' ','60^\circ',' ','90^\circ'})
end
%%
range_                  =   linspace(min(range),max(range),Ntick); 
rticks(range_)
Labels                  =   cell(1,Ntick);
for i=1:Ntick
    if i==Ntick
        Labels{1,i}  	=   strcat(num2str(round(range_(i),2)),' dB');
    else
        Labels{1,i}    	=   num2str(round(range_(i),2));
    end
end
rticklabels(Labels)
%%
pax = gca;
pax.RAxisLocation   =   0;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
end