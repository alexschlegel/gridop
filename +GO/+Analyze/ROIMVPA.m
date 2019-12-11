function [res,stat] = ROIMVPA(varargin)
% GO.Analyze.ROIMVPA
% 
% Description:	perform an ROI classification of shapes and operations
% 
% Syntax:	[res,stat] = GO.Analyze.ROIMVPA(<options>)
% 
% In:
% 	<options>:
%		tag:		('') an optional tag to identify the analysis
%		subject:	(<all>) the subjects to include
%		mask:		(<core>) the names of the masks to use
%		dim:		(<auto>) the PCA dimensionality
%		ifo:		(<load>) the result of a call to GO.SubjectInfo
%		cores:		(12) the number of processor cores to use
%		load:		(true) true to load the results if we previously saved them
%		force_mvpa:	(true) true to force classification
%		force_each:	(false) true to force each mask computation
%		force_pre:	(false) true to force preprocessing steps
%		silent:		(false) true to suppress status messages
% 
% Out:
% 	res		- the MVPA results
%	stat	- extra stats on the MVPA results
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'tag'			, ''	, ...
		'subject'		, {}	, ...
		'mask'			, {}	, ...
		'dim'			, []	, ...
		'ifo'			, []	, ...
		'cores'			, 12	, ...
		'load'			, true	, ...
		'force_mvpa'	, true	, ...
		'force_each'	, false	, ...
		'force_pre'		, false	, ...
		'silent'		, false	  ...
		);

strDirOut	= GO.Data.Directory(['roimvpa' opt.tag]);

%subject codes
	cSubject	= GO.Subject('subject',opt.subject);
	nSubject	= numel(cSubject);
%masks
	[cPathMask,cMask]	= GO.Path.Mask('subject',cSubject,'mask',opt.mask);
	nMask				= numel(cMask);

%have we done this already?
	param	= {opt.tag cSubject cMask opt.dim};
	
	if opt.load
		sData	= GO.Data.Load('roimvpa',param);
		
		if ~isempty(sData)
			res		= sData.res;
			stat	= sData.stat;
			
			return;
		end
	end

%data paths
	[C,cPathData,cPathMask]	= GO.Analyze.PCA(...
								'subject'	, cSubject		, ...
								'mask'		, cMask			, ...
								'dim'		, opt.dim		, ...
								'cores'		, opt.cores		, ...
								'load'		, opt.load		, ...
								'force'		, opt.force_pre	, ...
								'silent'	, opt.silent	  ...
								);

%get the subject info
	if isempty(opt.ifo)
		ifo	= GO.SubjectInfo('subject',cSubject);
	else
		ifo	= opt.ifo;
	end

%run each classification
	cScheme	= GO.Param('scheme');
	nScheme	= numel(cScheme);
	
	durRun	= GO.Param('trrun');
	nRun	= size(ifo.shape,2);
	kRun	= reshape(repmat(1:nRun,[durRun 1]),[],1);
	
	res	= struct;
	for kS=1:nScheme
		strScheme	= cScheme{kS};
		
		cTarget	= ifo.label.mvpa.target.(strScheme).correct;
		kChunk	= ifo.label.mvpa.chunk.correct;
		
		cOutPrefix	= cellfun(@(s) [s '-' strScheme],cSubject,'uni',false);
		
		res.(strScheme)	= MVPAClassify(cPathData,cTarget,kChunk,...
							'path_mask'			, cPathMask			, ...
							'mask_name'			, cMask				, ...
							'spatiotemporal'	, true				, ...
							'target_blank'		, 'Blank'			, ...
							'zscore'			, kRun				, ...
							'output_dir'		, strDirOut			, ...
							'output_prefix'		, cOutPrefix		, ...
							'cores'				, opt.cores			, ...
							'debug'				, 'all'				, ...
							'force'				, opt.force_mvpa	, ...
							'force_each'		, opt.force_each	, ...
							'silent'			, opt.silent		  ...
							);
	end

%calculate some extra stats
	conf	= GO.ConfusionModels;
	
	stat	= MVPAClassifyExtraStats(res,...
				'confusion_model'	, conf			, ...
				'silent'			, opt.silent	  ...
				);
	
%save the result
	sData.res			= res;
	sData.stat			= stat;
	
	GO.Data.Save(sData,'roimvpa',param);
