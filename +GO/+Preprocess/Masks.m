function b = Masks(varargin)
% GO.Masks
% 
% Description:	prepare masks from CI, an anatomical occipital mask, and a
%				ventricle mask for control analyses
% 
% Syntax:	b = GO.Masks(<options>)
% 
% In:
% 	<options>:
%		subject:	(<all>) the subjects to process
%		force:		(false)
%		cores:		(12)
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirData;

opt	= ParseArgs(varargin,...
		'subject'	, {}	, ...
		'force'		, false	, ...
		'cores'		, 12	  ...
		);

cSubject	= GO.Subject('subject',opt.subject,'state','fmri');

strDirMaskOut	= DirAppend(strDirData,'mask');

%CI masks
	strDirMaskCI	= DirAppend(strDirData,'mni-mask');
	cPathMaskMNI	= FindFilesByExtension(strDirMaskCI,'nii.gz');
	
	%remove the occ masks
		bRemove					= cellfun(@(f) ~isempty(strfind(PathGetFilePre(f,'favor','nii.gz'),'occ')),cPathMaskMNI);
		cPathMaskMNI(bRemove)	= [];
	
	nMask			= numel(cPathMaskMNI);
	
	cDirReg		= cellfun(@(s) DirAppend(strDirData,'functional',s,'feat_cat','reg'),cSubject,'uni',false);
	bDo			= cellfun(@isdir,cDirReg);
	
	[cSubject,cDirReg]	= varfun(@(x) x(bDo),cSubject,cDirReg);
	nSubject			= numel(cSubject);
	
	cDirMask	= cellfun(@(s) DirAppend(strDirMaskOut,s),cSubject,'uni',false);
	cellfun(@CreateDirPath,cDirMask);
	
	cPathXFMMNI2Func	= cellfun(@(d) PathUnsplit(d,'standard2example_func','mat'),cDirReg,'uni',false);
	cPathFunc			= cellfun(@(d) PathUnsplit(d,'example_func','nii.gz'),cDirReg,'uni',false);
	
	cPathMaskMNIRep		= repmat(cPathMaskMNI,[1 nSubject]);
	cPathXFMMNI2FuncRep	= repmat(cPathXFMMNI2Func',[nMask 1]);
	cPathFuncRep		= repmat(cPathFunc',[nMask 1]);
	cDirMaskRep			= repmat(cDirMask',[nMask 1]);
	
	cPathMaskRep	= cellfun(@(dm,fm) PathUnsplit(dm,PathGetFileName(fm)),cDirMaskRep,cPathMaskMNIRep,'uni',false);
	
	b	= FSLRegisterFLIRT(cPathMaskMNIRep,cPathFuncRep,...
			'output'	, cPathMaskRep			, ...
			'xfm'		, cPathXFMMNI2FuncRep	, ...
			'interp'	, 'nearestneighbour'	, ...
			'force'		, opt.force				, ...
			'cores'		, opt.cores				  ...
			);

%occipital mask
	cMaskLabel	=	{
						{'ctx_lh_G_and_S_occipital_inf' 'ctx_lh_G_occipital_middle' 'ctx_lh_G_occipital_sup' 'ctx_lh_G_cuneus' 'ctx_lh_Pole_occipital' 'ctx_lh_S_oc_middle_and_Lunatus' 'ctx_lh_S_oc_sup_and_transversal' 'ctx_lh_S_occipital_ant'}
						{'ctx_rh_G_and_S_occipital_inf' 'ctx_rh_G_occipital_middle' 'ctx_rh_G_occipital_sup' 'ctx_rh_G_cuneus' 'ctx_rh_Pole_occipital' 'ctx_rh_S_oc_middle_and_Lunatus' 'ctx_rh_S_oc_sup_and_transversal' 'ctx_rh_S_occipital_ant'}
						{'ctx_lh_G_and_S_occipital_inf' 'ctx_lh_G_occipital_middle' 'ctx_lh_G_occipital_sup' 'ctx_lh_G_cuneus' 'ctx_lh_Pole_occipital' 'ctx_lh_S_oc_middle_and_Lunatus' 'ctx_lh_S_oc_sup_and_transversal' 'ctx_lh_S_occipital_ant' 'ctx_rh_G_and_S_occipital_inf' 'ctx_rh_G_occipital_middle' 'ctx_rh_G_occipital_sup' 'ctx_rh_G_cuneus' 'ctx_rh_Pole_occipital' 'ctx_rh_S_oc_middle_and_Lunatus' 'ctx_rh_S_oc_sup_and_transversal' 'ctx_rh_S_occipital_ant'}
					};
	cCrop		=	{
						[]
						[]
						[]
					};
	cMaskName	=	{
						'occ-left'
						'occ-right'
						'occ'
					};
	nMask		= numel(cMaskLabel);
	
	cDirFEAT		= reshape(cellfun(@(s) DirAppend(strDirData,'functional',s,'feat_cat'),cSubject,'uni',false),1,[]);
	cDirFreeSurfer	= reshape(cellfun(@(s) DirAppend(strDirData,'structural',s,'freesurfer'),cSubject,'uni',false),1,[]);
	
	%compute the Freesurfer to functional space transforms
		b	= FreeSurfer2FEAT(cDirFreeSurfer,cDirFEAT,...
				'force'	, opt.force	, ...
				'cores'	, opt.cores	  ...
				);
		
		cPathFS2F	= cellfun(@(d) PathUnsplit(DirAppend(d,'reg'),'freesurfer2example_func','mat'),cDirFEAT,'UniformOutput',false);
		cPathF		= cellfun(@(d) PathUnsplit(DirAppend(d,'reg'),'example_func','nii.gz'),cDirFEAT,'UniformOutput',false);
	%extract the masks
		cDirFreeSurferR	= repmat(cDirFreeSurfer,[nMask 1]);
		cDirMaskR		= repmat(cDirMask',[nMask 1]);
		cMaskLabelR		= repmat(cMaskLabel,[1 nSubject]);
		cMaskNameR		= repmat(cMaskName,[1 nSubject]);
		cCropR			= repmat(cCrop,[1 nSubject]);
		cPathFS2FR		= repmat(cPathFS2F,[nMask 1]);
		cPathFR			= repmat(cPathF,[nMask 1]);
		
		cPathMaskOcc	= cellfun(@(d,m) PathUnsplit(d,m,'nii.gz'),cDirMaskR,cMaskNameR,'uni',false);
		
		cInput	=	{...
						cDirFreeSurferR						, ...
						cMaskLabelR							, ...
						'crop'			, cCropR			, ...
						'xfm'			, cPathFS2FR		, ...
						'ref'			, cPathFR			, ...
						'output'		, cPathMaskOcc		, ...
						'force'			, opt.force			  ...
					};
		
		b	= MultiTask(@FreeSurferMask,cInput,...
					'description'	, 'extracting masks'	, ...
					'cores'			, opt.cores				  ...
					);

%ventricle masks
	%construct the ventricle masks in freesurfer space
		[bSuccess,cPathVentricle]	= FreeSurferMaskVentricle(cDirFreeSurfer,...
										'force'	, opt.force	, ...
										'cores'	, opt.cores	  ...
										);
		cPathVentricle				= reshape(cPathVentricle,[],1);
	
	%convert to functional space
		cPathVFunc	= cellfun(@(d) PathUnsplit(d,'ventricle','nii.gz'),cDirMask,'uni',false);
		
		b	= FSLRegisterFLIRT(cPathVentricle,cPathFunc,...
			'output'	, cPathVFunc				, ...
			'xfm'		, reshape(cPathFS2F,[],1)	, ...
			'interp'	, 'nearestneighbour'		, ...
			'force'		, opt.force					, ...
			'cores'		, opt.cores					  ...
			);

%union masks
	sUnion	= GO.UnionMasks;
	
	cUnionName	= fieldnames(sUnion);
	nUnion		= numel(cUnionName);
	
	MultiTask(@(n) UnionMask(sUnion.(n),n),{cUnionName},...
		'description'	, 'constructing union masks'	, ...
		'cores'			, opt.cores						  ...
		);
	
%------------------------------------------------------------------------------%
function UnionMask(cMask,strName)
	cHemi	= {'-left';'-right';''};
	nHemi	= numel(cHemi);
	
	cMask	= reshape(cMask,[],1);
	
	for kH=1:nHemi
		strHemi	= cHemi{kH};
		
		for kS=1:nSubject
			strSubject	= cSubject{kS};
			
			strDirMaskCur	= DirAppend(strDirMaskOut,strSubject);
			cPathMask		= cellfun(@(m) PathUnsplit(strDirMaskCur,[m strHemi],'nii.gz'),cMask,'uni',false);
			strPathUnion	= PathUnsplit(strDirMaskCur,[strName strHemi],'nii.gz');
			
			MRIMaskMerge(cPathMask,strPathUnion,'force',opt.force,'silent',true);
		end
	end
end
%------------------------------------------------------------------------------%

end
