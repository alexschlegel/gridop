function [res,stat,cMaskPair] = CCMVPA(varargin)
% GO.Analyze.CCMVPA
% 
% Description:	perform an ROI cross-classification (between each pair of ROIs)
%				of shapes and operations
% 
% Syntax:	[res,stat,cMaskPair] = GO.Analyze.CCMVPA(<options>)
% 
% In:
% 	<options>:
%		tag:		('') an optional tag to identify the analysis
%		subject:	(<all>) the subjects to include
%		mask:		(<core>) the names of the masks to use
%		dim:		(50) the number of PCA/ICA components to use
%		ica:		(false) true to use ICA
%		ifo:		(<load>) the result of a call to GO.SubjectInfo
%		cores:		(12) the number of processor cores to use
%		load:		(true) true to load the results if we previously saved them
%		force_mvpa:	(true) true to force classification
%		force_each:	(false) true to force each mask computation
%		force_pre:	(false) true to force preprocessing steps
%		silent:		(false) true to suppress status messages
% 
% Out:
% 	res			- the MVPA results
%	stat		- extra stats on the MVPA results
%	cMaskPair	- an nMaskPair x 2 cell of mask pairs for each classification
%				  result
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'tag'			, ''	, ...
		'subject'		, {}	, ...
		'mask'			, {}	, ...
		'dim'			, 50	, ...
		'ica'			, false	, ...
		'ifo'			, []	, ...
		'cores'			, 12	, ...
		'load'			, true	, ...
		'force_mvpa'	, true	, ...
		'force_each'	, false	, ...
		'force_pre'		, false	, ...
		'silent'		, false	  ...
		);

strDirOut	= GO.Data.Directory(['ccmvpa' opt.tag]);

%subject codes
	cSubject	= GO.Subject('subject',opt.subject);
	nSubject	= numel(cSubject);
%masks
	[cPathMask,cMask]	= GO.Path.Mask('subject',cSubject,'mask',opt.mask);
	nMask				= numel(cMask);

%have we done this already?
	param	= {opt.tag cSubject cMask opt.dim opt.ica};
	
	if opt.load
		sData	= GO.Data.Load('ccmvpa',param);
		
		if ~isempty(sData)
			res			= sData.res;
			stat		= sData.stat;
			cMaskPair	= sData.cMaskPair;
			
			return;
		end
	end

%get the subject info
	if isempty(opt.ifo)
		ifo	= GO.SubjectInfo('subject',cSubject);
	else
		ifo	= opt.ifo;
	end

%align the ROIs to a common space
	[cPathData,cMaskPair]	= GO.Analyze.ROIAlign(...
								'subject'	, cSubject		, ...
								'mask'		, cMask			, ...
								'dim'		, opt.dim		, ...
								'ica'		, opt.ica		, ...
								'ifo'		, ifo			, ...
								'cores'		, opt.cores		, ...
								'load'		, opt.load		, ...
								'force'		, opt.force_pre	, ...
								'silent'	, opt.silent	  ...
								);
	nMaskPair				= size(cMaskPair,1);

%run each cross-classification
	cScheme	= GO.Param('scheme');
	nScheme	= numel(cScheme);
	
	durRun	= GO.Param('trrun');
	nRun	= size(ifo.shape,2);
	kRun	= reshape(repmat(1:nRun,[durRun 1]),[],1);
	
	conf	= GO.ConfusionModels;
	nModel	= numel(conf);
	
	[res,stat]	= deal(struct);
	for kS=1:nScheme
		strScheme	= cScheme{kS};
		
		cTarget	= repmat(ifo.label.mvpa.target.(strScheme).correct,[1 nMaskPair]);
		kChunk	= repmat(ifo.label.mvpa.chunk.correct,[1 nMaskPair]);
		
		cOutPrefix	= cellfun(@(f) [getfield(regexp(PathGetFilePre(f,'favor','nii.gz'),'(?<prefix>.+)-[^-]+$','names'),'prefix') '-' strScheme],cPathData(:,:,1),'uni',false);
		
		res.(strScheme)	= MVPAClassify(cPathData,cTarget,kChunk,...
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
		
		stat.(strScheme).acc	= reshape(res.(strScheme).allway.accuracy.mean,[nSubject nMaskPair]);
		stat.(strScheme).mAcc	= mean(stat.(strScheme).acc)';
		stat.(strScheme).seAcc	= stderr(stat.(strScheme).acc)';
		
		stat.(strScheme).conf	= reshape(res.(strScheme).allway.confusion,[4 4 nSubject nMaskPair]);
		stat.(strScheme).mConf	= squeeze(mean(stat.(strScheme).conf,3));
		stat.(strScheme).seConf	= squeeze(stderr(stat.(strScheme).conf,[],3));
		
		[h,p,ci,stats]				= ttest(stat.(strScheme).acc,0.25,'tail','right');
		[pThresh,pFDR]				= fdr(p,0.05);
		stat.(strScheme).pAcc		= p';
		stat.(strScheme).pfdrAcc	= pFDR';
		stat.(strScheme).tAcc		= stats.tstat';
		stat.(strScheme).dfAcc		= stats.df';
		
		[r,stats]					= corrcoef2(reshape(conf{1},[],1),reshape(permute(stat.(strScheme).mConf,[3 1 2]),nMaskPair,[]));
		[pThresh,pFDR]				= fdr(stats.p,0.05);
		stat.(strScheme).pConf		= stats.p;
		stat.(strScheme).pfdrConf	= pFDR;
		stat.(strScheme).rConf		= stats.r;
		stat.(strScheme).dfConf		= stats.df;
		
		for kM=1:nModel
			[r,stats]									= corrcoef2(reshape(conf{kM},[],1),reshape(permute(stat.(strScheme).mConf,[3 1 2]),nMaskPair,[]));
			[pThresh,pFDR]								= fdr(stats.p,0.05);
			stat.(strScheme).modelcompare(kM).pConf		= stats.p;
			stat.(strScheme).modelcompare(kM).pfdrConf	= pFDR;
			stat.(strScheme).modelcompare(kM).rConf		= stats.r;
			stat.(strScheme).modelcompare(kM).dfConf	= stats.df;
		end
	end

%save the result
	sData.res		= res;
	sData.stat		= stat;
	sData.cMaskPair	= cMaskPair;
	
	GO.Data.Save(sData,'ccmvpa',param);
