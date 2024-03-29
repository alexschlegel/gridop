function Organize(varargin)
% GO.Preprocess.Organize
% 
% Description:	organize the raw data
% 
% Syntax:	GO.Preprocess.Organize(<options>)
% 
% In:
% 	<options>:
%		cores:	(12)
%		force:	(false)
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirData

opt	= ParseArgs(varargin,...
		'cores'	, 12	, ...
		'force'	, false	  ...
		);

strDirRaw	= DirAppend(strDirData,'raw');

%organize the data
	[b,cPathRaw,cPathOut]	= PARRECOrganize(strDirRaw,...
								'cores'	, opt.cores	, ...
								'force'	, opt.force	  ...
								);
