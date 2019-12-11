% Analysis_20140723_Control_ROI
% address reviewer concerns from the Science submission:
%	1)	are there areas where we can't classify? as a control analysis, try to
%		classify in THAL, the ROI that showed the lowest classification
%		accuracies in our previous study, and the ventricles.
nCore			= 12;
strNameAnalysis	= '20140723_control_roi';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;
cSubject	= ifo.code.fmri;

dimPCAGC	= 10;
dimPCACC	= 50;

%ROI control
	cMask	= {
				'thal'
				'ventricle'
				};
	nMask		= numel(cMask);
	
	cMaskLabel		= cellfun(@(m) upper(m(1:min(4,numel(m)))),cMask,'uni',false);
	
	[res,stat]	= GO.Analyze.ROIMVPA(...
					'subject'		, cSubject	, ...
					'mask'			, cMask		, ...
					'mindim'		, dimPCAGC	, ...
					'ifo'			, ifo		, ...
					'cores'			, nCore		  ...
					);
	
	%figures
		h	= GO.Plot.Accuracy(stat,'outdir',strDirOut);
		h	= GO.Plot.ConfusionCorrelation(stat,'outdir',strDirOut,'jackknife',true);
		h	= GO.Plot.Confusion(stat,'outdir',strDirOut);
