function [cPathData,cMaskPair] = ROIAlign(varargin)
% GO.Analyze.ROIAlign
% 
% Description:	align PCA/ICA components between ROIs
% 
% Syntax:	[cPathData,cMaskPair] = ROIAlign(<options>)
% 
% In:
% 	<options>:
%		subject:	(<all>) the subjects to include
%		mask:		(<core>) the names of the masks to use
%		dim:		(50) the number of PCA/ICA components to use
%		ica:		(false) true to use ICA
%		ifo:		(<load>) the result of a call to GO.SubjectInfo
%		cores:		(12) the number of processor cores to use
%		load:		(true) true to load the results if we previously saved them
%		force:		(false) true to force computation
%		force_pre:	(false) true to force preprocessing steps
%		silent:		(false) true to suppress status messages
% 
% Out:
%	cPathData	- an nSubject x nMaskPair x 2 cell of aligned PCA/ICA NIfTI
%				  files
%	cMaskPair	- an nMaskPair x 2 cell of mask pairs
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'subject'	, {}	, ...
		'mask'		, {}	, ...
		'dim'		, 50	, ...
		'ica'		, false	, ...
		'ifo'		, []	, ...
		'cores'		, 12	, ...
		'load'		, true	, ...
		'force'		, false	, ...
		'force_pre'	, false	, ...
		'silent'	, false	  ...
		);

strDirOut	= GO.Data.Directory('roialign');

%subject codes
	cSubject	= GO.Subject('subject',opt.subject);
	nSubject	= numel(cSubject);
%mask paths
	[cPathMask,cMask]	= GO.Path.Mask('subject',opt.subject,'mask',opt.mask);
	nMask				= numel(cMask);

%have we done this already?
	param	= {cSubject cMask opt.dim opt.ica};
	
	if ~opt.force && opt.load
		sData	= GO.Data.Load('roialign',param);
		
		if ~isempty(sData)
			cPathData	= sData.cPathData;
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

%data paths
	[C,cPathPCA,cPathMask]	= GO.Analyze.PCA(...
								'subject'	, cSubject		, ...
								'mask'		, cMask			, ...
								'dim'		, opt.dim		, ...
								'ica'		, opt.ica		, ...
								'cores'		, opt.cores		, ...
								'load'		, opt.load		, ...
								'force'		, opt.force_pre	, ...
								'silent'	, opt.silent	  ...
								);

%align the PCA/ICA components for each mask pair. do the alignment for each
%trial, leaving data for that trial out of the alignment calculation
	[cMaskPair,kMaskPair]	= handshakes(cMask);
	cMaskPair				= [cMaskPair; repmat(cMask,[1 2])];
	kMaskPair				= [kMaskPair; repmat((1:nMask)',[1 2])];
	nMaskPair				= size(cMaskPair,1);
	
	%get the output file paths
		cPathData	= cell(nSubject,nMaskPair,2);
		
		for kS=1:nSubject
			strSubject	= cSubject{kS};
			
			for kM=1:nMaskPair
				strMask1	= cMaskPair{kM,1};
				strMask2	= cMaskPair{kM,2};
				
				cPathData{kS,kM,1}	= PathUnsplit(strDirOut,sprintf('%s-%s_%s-%s',strSubject,strMask1,strMask2,strMask1),'nii.gz');
				cPathData{kS,kM,2}	= PathUnsplit(strDirOut,sprintf('%s-%s_%s-%s',strSubject,strMask1,strMask2,strMask2),'nii.gz');
			end
		end
	
	b	= MATLABPoolOpen(opt.cores);
	
	progress('action','init','total',nSubject,'label','aligning PCA/ICA components for each subject','silent',opt.silent);
	for kS=1:nSubject
		kChunk	= ifo.label.mvpa.chunk.correct{kS};
		
		%do this to make debugging easier
		if opt.cores > 1
			parfor kM=1:nMaskPair
				if opt.force || ~FileExists(cPathData{kS,kM,1}) || ~FileExists(cPathData{kS,kM,2})
					kM1	= kMaskPair(kM,1);
					kM2	= kMaskPair(kM,2);
					
					%PCA/ICA subsets
						C1	= C{kS,kM1}(:,1:opt.dim);
						C2	= C{kS,kM2}(:,1:opt.dim);
					%match up the two sets of components, chunk by chunk (to make sure
					%we aren't artifically inflating the similarity between signals)
						[C1,C2]	= SignalMatch(C1,C2,'chunk',kChunk);
					
					%save the components to a NIfTI file
						nii1	= NIfTI.Create(permute(C1,[2 3 4 1]));
						nii2	= NIfTI.Create(permute(C2,[2 3 4 1]));
						
						NIfTI.Write(nii1,cPathData{kS,kM,1});
						NIfTI.Write(nii2,cPathData{kS,kM,2});
				end
			end
		else
			progress('action','init','total',nMaskPair,'name','maskpair','label','aligning PCA/ICA components for each mask pair','silent',opt.silent);
			for kM=1:nMaskPair
				if opt.force || ~FileExists(cPathData{kS,kM,1}) || ~FileExists(cPathData{kS,kM,2})
					kM1	= kMaskPair(kM,1);
					kM2	= kMaskPair(kM,2);
					
					%PCA/ICA subsets
						C1	= C{kS,kM1}(:,1:opt.dim);
						C2	= C{kS,kM2}(:,1:opt.dim);
					%match up the two sets of components, chunk by chunk (to make sure
					%we aren't artifically inflating the similarity between signals)
						try
							[C1,C2]	= SignalMatch(C1,C2,'chunk',kChunk);
						catch me
							rethrow(me);
						end
					
					%save the components to a NIfTI file
						nii1	= NIfTI.Create(permute(C1,[2 3 4 1]));
						nii2	= NIfTI.Create(permute(C2,[2 3 4 1]));
						
						NIfTI.Write(nii1,cPathData{kS,kM,1});
						NIfTI.Write(nii2,cPathData{kS,kM,2});
				end
				
				progress('name','maskpair');
			end
		end
		
		progress;
	end
	
	b	= MATLABPoolClose;

%save the result
	sData.cPathData	= cPathData;
	sData.cMaskPair	= cMaskPair;
	
	GO.Data.Save(sData,'roialign',param);
