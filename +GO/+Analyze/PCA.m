function [C,cPathPCA,cPathMaskPCA] = PCA(varargin)
% GO.Analyze.PCA
% 
% Description:	run FSL's MELODIC tool on gridop functional data
% 
% Syntax:	[C,cPathPCA,cPathMaskPCA] = GO.Analyze.PCA(<options>)
% 
% In:
% 	<options>:
%		subject:(<all>) the subjects to include
%		mask:	(<core>) the names of the masks to use
%		mindim:	(10) the minimum number of PCA dimensions
%		dim:	([]) manually set the number of PCA dimensions
%		ica:	(false) true to return ICA components rather than PCA
%		cores:	(12) the number of processor cores to use
%		load:	(true) true to load the results if we previously saved them
%		force:	(false) true to force computation
%		silent:	(false) true to suppress status messages
% 
% Out:
% 	C				- an nSubject x nMask cell of PCA/ICA component signals
%	cPathPCA		- an nSubject x 1 cell of PCA/ICA NIfTI files
%	cPathMaskPCA	- an nSubject x 1 cell of nMask x 1 cells of PCA/ICA mask
%					  files
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'subject'	, {}	, ...
		'mask'		, {}	, ...
		'mindim'	, 10	, ...
		'dim'		, []	, ...
		'ica'		, false	, ...
		'cores'		, 12	, ...
		'load'		, true	, ...
		'force'		, false	, ...
		'silent'	, false	  ...
		);

strTx		= conditional(opt.ica,'ica','pca');
strDim		= tostring(unless(opt.dim,'auto'));
strDataName	= sprintf('pca-%s',strDim);

status(sprintf('Using %s data (%s, mindim=%d)',strDataName,strTx,opt.mindim),'silent',opt.silent);
strDirOut	= GO.Data.Directory(strDataName);

%subject codes
	cSubject	= GO.Subject('subject',opt.subject);
	nSubject	= numel(cSubject);
%mask paths
	[cPathMask,cMask]	= GO.Path.Mask('subject',opt.subject,'mask',opt.mask);
	nMask				= numel(cMask);

%have we done this already?
	param	= {cSubject cMask opt.mindim opt.ica};
	
	if ~opt.force && opt.load
		sData	= GO.Data.Load(strDataName,param);
		
		if ~isempty(sData)
			C				= sData.C;
			cPathPCA		= sData.cPathPCA;
			cPathMaskPCA	= sData.cPathMaskPCA;
			
			return;
		end
	end

%get the functional data paths
	cPathData	= GO.Path.Functional('subject',cSubject);
	
%call MELODIC
	cPathData	= repmat(cPathData,[1 nMask]);
	cPathMask	= cat(2,cPathMask{:})';
	
	cSubjectRep		= repmat(cSubject,[1 nMask]);
	cMaskRep	 	= repmat(cMask',[nSubject 1]);
	
	cDirMELODIC	= cellfun(@(s,m) DirAppend(strDirOut,[s '-' m]),cSubjectRep,cMaskRep,'uni',false);
	
	s	= FSLMELODIC(cPathData,...
			'path_functional'	, cPathData		, ...
			'path_mask'			, cPathMask		, ...
			'dir_out'			, cDirMELODIC	, ...
			'mindim'			, opt.mindim	, ...
			'dim'				, opt.dim		, ...
			'comptype'			, strTx			, ...
			'load_weight'		, false			, ...
			'cores'				, opt.cores		, ...
			'force'				, opt.force		, ...
			'silent'			, opt.silent	  ...
			);
	C	= s.comp;

%save the NIfTI files
	strSuffix	= join(cMask,'');
	
	cPathPCA		= cellfun(@(s) PathUnsplit(strDirOut,[s '-' strSuffix],'nii.gz'),cSubject,'uni',false);
	cPathMaskPCA	= cellfun(@(s) cellfun(@(m) PathUnsplit(strDirOut,[s '-' m],'nii.gz'),cMask,'uni',false),cSubject,'uni',false);
	
	progress('action','init','total',nSubject,'label','saving PCA/ICA data');
	for kS=1:nSubject
		%make the data file
			data	= cat(2,C{kS,:});
			nData	= size(data,2);
			
			data	= permute(data,[2 3 4 1]);
			nii		= NIfTI.Create(data);
			
			NIfTI.Write(nii,cPathPCA{kS});
		%make the masks
			kMaskPre	= 0;
			for kM=1:nMask
				nCompMask	= size(C{kS,kM},2);
				kMaskStart	= kMaskPre+1;
				kMaskEnd	= kMaskStart + nCompMask - 1;
				
				msk							= zeros(nData,1);
				msk(kMaskStart:kMaskEnd)	= 1;
				
				niiMask	= NIfTI.Create(msk);
				
				NIfTI.Write(niiMask,cPathMaskPCA{kS}{kM});
				
				kMaskPre	= kMaskEnd;
			end
		
		progress;
	end

%save the result
	sData.C				= C;
	sData.cPathPCA		= cPathPCA;
	sData.cPathMaskPCA	= cPathMaskPCA;
	
	GO.Data.Save(sData,strDataName,param);
