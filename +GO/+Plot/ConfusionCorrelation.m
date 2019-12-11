function h = ConfusionCorrelation(stat,varargin)
% GO.Plot.ConfusionCorrelation
% 
% Description:	plot confusion correlation figures using info from a stat
%				struct returned by GO.Analyze.ROIMVPA
% 
% Syntax:	h = GO.Plot.ConfusionCorrelation(stat,<options>)
% 
% In:
%	stat	- the stat struct returned from GO.Analyze.ROIMVPA
%	<options>:
%		jackknife:	(false) true to construct the plot based on the jackknife
%					statistics
%		models:		(<all>) the indices of the models to plot
%		ymin:		(0) the minimum y value
%		ymax:		(3.5) the maximum y value
%		<other>:	extra options to GO.Plot.Bar
% 
% Updated: 2014-07-28
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'jackknife'	, false	, ...
		'models'	, []	, ...
		'ymin'		, 0		, ...
		'ymax'		, 3.5	  ...
		);

if isempty(opt.models)
	nModel	= size(stat.confusion.corrcompare.group.allway.r,3);;
	kModel	= (1:nModel)';
else
	kModel	= opt.models;
	nModel	= numel(kModel);
end

cLabelGroup	= stat.label{1};
cLabelBar	= stat.label{2};

for kM=1:nModel
	kMCur	= kModel(kM);
	
	if opt.jackknife
		z		= fisherz(stat.confusion.corrcompare.subjectJK.allway.group.r(:,:,kMCur));
		err		= fisherz(stat.confusion.corrcompare.subjectJK.allway.group.se(:,:,kMCur));
		p		= stat.confusion.corrcompare.subjectJK.allway.group.p(:,:,kMCur);
		zThresh	= [];
	else
		z	= fisherz(stat.confusion.corrcompare.group.allway.r(:,:,kMCur));
		err	= [];
		p	= stat.confusion.corrcompare.group.allway.p(:,:,kMCur);
		
		zThresh	= fisherz(stat.confusion.corrcompare.group.allway.cutoff);
	end
	
	strName	= sprintf('confusioncorrelation%d',kMCur);
	
	h	= GO.Plot.Bar(z,err,p,...
			'name'		, strName			, ...
			'groups'	, cLabelGroup		, ...
			'bars'		, cLabelBar			, ...
			'ylabel'	, 'Fisher''s Z(r)'	, ...
			'ymin'		, opt.ymin			, ...
			'ymax'		, opt.ymax			, ...
			'thresh'	, zThresh			, ...
			varargin{:});
end
