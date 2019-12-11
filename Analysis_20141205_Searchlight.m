%PNAS reviewer 2 wants a searchlight done to explore how widespread the effects were
%prepare the parameters for searchlight shape and operation classifications
PrepGO

strName	= '20141205_searchlight';

bForce	= false;
nCore	= 11;

ifo = GO.SubjectInfo;

param	= struct;

%subject info
	param.subject	= ifo.code.fmri;
	nSubject		= numel(param.subject);

%output directory
	param.output_dir	= DirAppend(strDirAnalysis,strName);
	CreateDirPath(param.output_dir);
	
%run for each sample
	durRun	= GO.Param('trrun');
	nRun	= size(ifo.shape,2);
	kRun	= reshape(repmat(1:nRun,[durRun 1]),[],1);

%attributes
	strDirAttr	= DirAppend(param.output_dir,'attr');
	CreateDirPath(strDirAttr);
	
	param.scheme	= {'shape'; 'operation'};
	nScheme			= numel(param.scheme);
	
	for kS=1:nScheme
		scheme	= param.scheme{kS};
		param.(scheme).attr	= cell(nSubject,1);
		
		for kU=1:nSubject
			subject	= param.subject{kU};
			
			attr.target	= ifo.label.mvpa.target.(scheme).all{kU};
			attr.chunk	= kRun;
			
			strAttr	= struct2table(attr,'heading',false);
			
			strPathAttr				= PathUnsplit(strDirAttr,sprintf('%s_%s',scheme,subject),'attr');
			param.(scheme).attr{kU}	= strPathAttr;
			
			if bForce || ~FileExists(strPathAttr)
				fput(strAttr,strPathAttr);
			end
		end
	end

%transform functional data to mni 3mm space
	cPathFunctional		= GO.Path.Functional('subject',param.subject);
	[b,param.data_path]	= FSLFunc2Standard3(cPathFunctional,...
							'force'	, bForce	, ...
							'cores'	, nCore		  ...
							);

%MNI gray matter mask
	strDirFreeSurfer	= DirAppend(strDirData,'mni-freesurfer');
	
	[b,strPathGM]	= FreeSurferMaskGM(strDirFreeSurfer,...
						'force'	, bForce	  ...
						);
	
	%transform freesurfer brain to standard 3mm space
		strDirMRI		= DirAppend(strDirFreeSurfer,'mri');
		strPathBrain	= PathUnsplit(strDirMRI,'brain','nii.gz');
		
		strDirReg	= DirAppend(PathGetDir(param.data_path{1}),'feat_cat','reg');
		strPathRef	= PathUnsplit(strDirReg,'standard-3mm','nii.gz');
		
		[b,strPathBrain3,strPathXFM3]	= FSLRegisterFLIRT(strPathBrain,strPathRef,...
											'force'		, bForce				  ...
											);
	
	%transform gm mask to FSL standard 3mm space
		strPathGM3	= PathAddSuffix(strPathGM,'-2mni-3mm','favor','nii.gz');
		
		b	= FSLRegisterFLIRT(strPathGM,strPathRef,...
				'output'	, strPathGM3			, ...
				'xfm'		, strPathXFM3			, ...
				'interp'	, 'nearestneighbour'	, ...
				'force'		, bForce				  ...
				);
	
	param.gm_path	= strPathGM3;

%save the parameters
	strPathParam	= PathUnsplit(strDirCode,'Analysis_20141205_Searchlight','json');
	
	param	= structtreefun(@(x) reshape(x,1,[]),param);
	fput(json.to(param),strPathParam);

%%%RUN THE PYTHON SCRIPT!!!%%%

%analyze the results
	%load the searchlight results
	cPathSL	= cellfun(@(sc) cellfun(@(su) PathUnsplit(param.output_dir,sprintf('%s_%s',su,sc),'nii'),param.subject','uni',false),param.scheme','uni',false);
	dSL		= cellfun(@(sc) cellfun(@(fsu) double(NIfTI.Read(fsu,'return','data')),sc,'uni',false),cPathSL,'uni',false);
	
	%load the mask
	msk	= logical(NIfTI.Read(param.gm_path,'return','data'));
	
	%perform a t-test on each voxel compared to chance
		dAnalyze	= cellfun(@(sc) cellfun(@(d) d(msk),sc,'uni',false),dSL,'uni',false);
		dAnalyze	= cellfun(@(sc) cat(2,sc{:}),dAnalyze,'uni',false);
		
		[h,p,ci,stats]	= cellfun(@(d) ttest(d,0.25,'dim',2,'tail','right'),dAnalyze,'uni',false);
		
		[dummy,pfdr]	= cellfun(@(p) fdr(p,0.05),p,'uni',false);
	
	%save the results
		strDirOut	= DirAppend(param.output_dir,'result');
		CreateDirPath(strDirOut);
		
		nii			= NIfTI.Read(cPathSL{1}{1});
		nii.data	= NaN(size(nii.data));
		
		for kS=1:nScheme
			strScheme	= param.scheme{kS};
			
			%acc
				nii.data(msk)	= mean(dAnalyze{kS},2);
				
				strPathOut	= PathUnsplit(strDirOut,sprintf('%s_%s',strScheme,'acc'),'nii.gz');
				NIfTI.Write(nii,strPathOut);
			
			%tstat
				nii.data(msk)	= stats{kS}.tstat;
				
				strPathOut	= PathUnsplit(strDirOut,sprintf('%s_%s',strScheme,'tstat'),'nii.gz');
				NIfTI.Write(nii,strPathOut);
			
			%1 - p
				nii.data(msk)	= 1 - p{kS};
				
				strPathOut	= PathUnsplit(strDirOut,sprintf('%s_%s',strScheme,'p'),'nii.gz');
				NIfTI.Write(nii,strPathOut);
			
			%1 - pfdr
				nii.data(msk)	= 1 - pfdr{kS};
				
				strPathOut	= PathUnsplit(strDirOut,sprintf('%s_%s',strScheme,'pfdr'),'nii.gz');
				NIfTI.Write(nii,strPathOut);
		end
	
	%copy over some helper files
		strPathGM	= PathUnsplit(strDirOut,'gm','nii.gz');
		FileCopy(param.gm_path,strPathGM);
		
		strPathMNIFrom	= PathUnsplit(DirAppend(PathGetDir(param.data_path{1}),'feat_cat','reg'),'standard-3mm','nii.gz');
		strPathMNITo	= PathUnsplit(strDirOut,'standard','nii.gz');
		FileCopy(strPathMNIFrom,strPathMNITo);
