%PNAS reviewer 2 wants us to verify that no univariate effects exist.
%1) perform a first level analysis fitting betas to each condition for each experiment
%2) perform a higher-level analysis, 1-way ANOVA with four levels and repeated measures
%		(see http://fsl.fmrib.ox.ac.uk/fsl/fsl4.0/feat5/detail.html#ANOVA1factor4levelsRepeatedMeasures)
PrepGO

strName		= '20141209_univariate';
strDirOut	= DirAppend(strDirAnalysis,strName);
CreateDirPath(strDirOut);

nCore	= 12;

ifo			= GO.SubjectInfo;
nSubject	= numel(ifo.id);
cScheme		= {'shape';'operation'};
nScheme		= numel(cScheme);

%first level analysis
	cPathFunctional	= GO.Path.Functional;
	cDirFEAT		= cellfun(@FSLDirFEAT,cPathFunctional,'uni',false);
	
	for kS=1:nScheme
		strScheme	= cScheme{kS};
		
		cDirFEATScheme	= cellfun(@(d) regexprep(d,'/([^/]+)/$',sprintf('/$1-%s/',strScheme)),cDirFEAT,'uni',false);
		
		cNameEV	= [ifo.condition.(strScheme); 'bad'];
		
		cLabel		= cell(nSubject,1);
		for kU=1:nSubject
			cLabelCorrect	= ifo.label.mvpa.target.(strScheme).correct{kU};
			cLabelAll		= ifo.label.mvpa.target.(strScheme).all{kU};
			
			bBad				= ~strcmp(cLabelCorrect,cLabelAll);
			cLabel{kU}			= cLabelCorrect;
			cLabel{kU}(bBad)	= {'bad'};
		end
		
		d	= cellfun(@(L) cellfun(@(ev) strcmp(L,ev),cNameEV,'uni',false),cLabel,'uni',false);
		d	= cellfun(@(ds) cat(2,ds{:}),d,'uni',false);
		
		[b,cDirFEATScheme]	= FSLFEATFirst(cPathFunctional,d,...
								'output'			, cDirFEATScheme	, ...
								'ev_name'			, cNameEV			, ...
								'convolve'			, true				, ...
								'tfilter'			, true				, ...
								'tcontrast_name'	, cNameEV			, ...
								'cores'				, nCore				, ...
								'force'				, false				  ...
								);
		
		%copy the reg folders over
			cDirRegFrom	= cellfun(@(d) DirAppend(d,'reg'),cDirFEAT,'uni',false);
			cDirRegTo	= cellfun(@(d) DirAppend(d,'reg'),cDirFEATScheme,'uni',false);
			b			= cellfunprogress(@FileCopy,cDirRegFrom,cDirRegTo);
	end

%higher level analysis
	cPathFunctional	= GO.Path.Functional;
	cSubject		= ifo.code.fmri;
	
	%remove xw, which for some reason that i am too lazy to determine is not
	%making it past the first level analysis
	cPathFunctional(18)	= [];
	cSubject(18)		= [];
	nSubject			= numel(cPathFunctional);
	
	cDirFEAT		= cellfun(@FSLDirFEAT,cPathFunctional,'uni',false);
	
	cDirOut	= cellfun(@(s) DirAppend(strDirOut,sprintf('anova_%s',s)),cScheme,'uni',false);
	
	cCondition		= cellfun(@(s) ifo.condition.(s),cScheme,'uni',false);
	nCondition		= cellfun(@numel,cCondition);
	
	cNameEVSubject		= cSubject;
	cNameEVCondition	= cellfun(@(cc) cc(1:end-1),cCondition,'uni',false);
	cNameEV				= cellfun(@(cnc) [cNameEVSubject; cnc],cNameEVCondition,'uni',false);
	
	cTContrast		= arrayfun(@(n) [zeros(n-1,nSubject) eye(n-1)],nCondition,'uni',false);
	cTContrastName	= cellfun(@(cc) cellfun(@(c) sprintf('%s-%s',c,cc{end}),cc(1:end-1),'uni',false),cCondition,'uni',false);
	
	fTest	= [1 1 1];
	
	[cPathCOPE,cD]	= deal(cell(nScheme,1));
	
	for kS=1:nScheme
		strScheme	= cScheme{kS};
		
		cDirFEATScheme	= cellfun(@(d) regexprep(d,'/([^/]+)/$',sprintf('/$1-%s/',strScheme)),cDirFEAT,'uni',false);
		cDirStats		= cellfun(@(d) DirAppend(d,'stats'),cDirFEATScheme,'uni',false);
		
		cPathCOPE{kS}	= arrayfun(@(k) cellfun(@(d) PathUnsplit(d,sprintf('cope%d',k),'nii.gz'),cDirStats,'uni',false),(1:nCondition(kS))','uni',false); 
		cPathCOPE{kS}	= cat(1,cPathCOPE{kS}{:});
		
		evSubject	= arrayfun(@(k) repmat([zeros(k-1,1); 1; zeros(nSubject-k,1)],[nCondition(kS) 1]),(1:nSubject)','uni',false);
		evSubject	= cat(2,evSubject{:});
		
		evCondition	= arrayfun(@(k) [zeros((k-1)*nSubject,1); ones(nSubject,1); zeros((nCondition(kS)-k)*nSubject,1)],(1:nCondition(kS)-1)','uni',false);
		evCondition	= cat(2,evCondition{:});
		
		cD{kS}	= [evSubject evCondition];
	end
	
	[b,cDirOut]	= FSLFEATHigher(cPathCOPE,cD,...
					'output'			, cDirOut			, ...
					'ev_name'			, cNameEV			, ...
					'tcontrast'			, cTContrast			, ...
					'tcontrast_name'	, cTContrastName	, ...
					'ftest'				, fTest				, ...
					'cores'				, nCore				, ...
					'force'				, false				  ...
					);
