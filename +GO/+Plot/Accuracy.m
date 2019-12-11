function h = Accuracy(stat,varargin)
% GO.Plot.Accuracy
% 
% Description:	plot an accuracy figure using info from a stat struct returned
%				by GO.Analyze.ROIMVPA
% 
% Syntax:	h = GO.Plot.Accuracy(stat,<options>)
% 
% Updated: 2014-03-06
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.

acc	= 100*stat.accuracy.mean.allway;
err	= 100*stat.accuracy.se.allway;
p	= stat.accuracy.pfdr.allway;

cLabelGroup	= stat.label{1};
cLabelBar	= stat.label{2};

%plot
	h	= GO.Plot.Bar(acc,err,p,...
			'name'		, 'accuracy'		, ...
			'groups'	, cLabelGroup		, ...
			'bars'		, cLabelBar			, ...
			'ylabel'	, 'Accuracy (%)'	, ...
			'ymin'		, 15				, ...
			'ymax'		, 50				, ...
			'thresh'	, 25				, ...
			varargin{:});
