%Analysis_20150307_NewCCTest.m
%test the new MVPAROICrossClassify script
nCore			= 12;
strNameAnalysis	= '20150307_newcctest';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;
cSubject	= ifo.code.fmri;
cMask		= GO.UnionMasks;
cMask		= [cMask.core];
nMask		= numel(cMask);
cMaskCore	= [cMask; 'core'];

cMaskLabel		= upper(cMask);
cMaskLabel{5}	= 'LOC';

dimPCAGC	= 10;
dimPCACC	= 50;

durRun	= GO.Param('trrun');
nRun	= size(ifo.shape,2);
kRun	= reshape(repmat(1:nRun,[durRun 1]),[],1);

cTarget	= ifo.label.mvpa.target.shape.correct;
kChunk	= ifo.label.mvpa.chunk.correct;

%***
%cMask						= {'occ','pcu'};
%[cSubject,cTarget,kChunk]	= varfun(@(x) x(1:10),cSubject,cTarget,kChunk);

res	= MVPAROICrossClassify(...
		'output_dir'		, strDirOut		, ...
		'dir_data'			, strDirData	, ...
		'subject'			, cSubject		, ...
		'mask'				, cMask			, ...
		'targets'			, cTarget		, ...
		'chunks'			, kChunk		, ...
		'spatiotemporal'	, true			, ...
		'target_blank'		, 'Blank'		, ...
		'zscore'			, kRun			, ...
		'debug'				, 'all'			, ...
		'cores'				, nCore			, ...
		'force'				, false			  ...
		);

conf	= GO.ConfusionModels;
conf	= conf{1};
	
stat	= MVPAClassifyExtraStats(res,...
			'confusion_model'	, conf	  ...
			);

strPathOut	= PathUnsplit(strDirOut,'result','mat');
save(strPathOut,'res','stat');
