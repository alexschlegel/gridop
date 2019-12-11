% Analysis_20150910_Control_RTMatch
% address JoCN reviewer concern about RTs. perform RT-matched shape ROI
% classifications.
nCore			= 12;
strNameAnalysis	= '20150910_control_rtmatch';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;

cSubject	= ifo.code.fmri;
nSubject	= numel(cSubject);

cMask		= GO.UnionMasks;
cMask		= [cMask.core];
nMask		= numel(cMask);

dimPCAGC	= 10;

%construct the targets/chunks so we have RT-matched conditions
	nRun	= GO.Param('exp','runs');
	
	%taken from GO.SubjectInfo
		sTargetBlank	= dealstruct('all','correct',cell(nSubject,1));
		sLabelBlank		= struct(...
							'chunk'		, sTargetBlank									, ...
							'target'	, structfun2(@(x) sTargetBlank, ifo.condition)	  ...
							);
		
		%mvpa:
		%	consider the middle three TRs of the trial. this gives us an effective
		%	HRF shift of 1 TR, and we don't consider any TR at or after the test
		%	screen.
		%te:
		%	get some more data in there for TE calculations. hopefully this gets the
		%	tail end of any information transfer. 
		sLabelScheme	= struct(...
							'mvpa'	, struct('HRF',1,'BlockOffset',0,'BlockSub',3)						, ...
							'te'	, struct('HRF',1,'BlockOffset',0,'BlockSub',GO.Param('te','block'))	  ...
							);
		cLabelScheme	= fieldnames(sLabelScheme);
		nLabelScheme	= numel(cLabelScheme);
		
		durBlock	= GO.Param('trtrial');
		durRest		= GO.Param('trrest');
		durPre		= GO.Param('trrestpre') - durRest;
		durPost		= GO.Param('trrestpost') - durRest;
		durRun		= GO.Param('trrun');
	
	strScheme	= 'shape';
	cCondition	= reshape(ifo.condition.(strScheme),[],1);
	nCondition	= numel(cCondition);
	kCondition	= (1:nCondition)';
	
	cConditionCI	= [repmat({'Blank'},[nCondition 1]); cCondition];
	
	sLabel	= sLabelScheme.mvpa;
	
	nTrialRecord	= NaN(nSubject,1);
	rtTrialRecord	= cell(nSubject,1);
	
	[cTargetRT,cChunkRT]			= deal(cell(nSubject,1));
	for kS=1:nSubject
		bCorrect	= reshape(ifo.correct(kS,:,:),[],1);
		rt			= reshape(ifo.rt(kS,:,:),[],1);
		
		kConditionTrial	= reshape(ifo.(strScheme)(kS,:,:),[],1);
		
		cKTrialByCondition			= arrayfun(@(kc) find(kConditionTrial==kc),kCondition,'uni',false);
		cKTrialByConditionCorrect	= cellfun(@(kt) kt(bCorrect(kt)==1),cKTrialByCondition,'uni',false);
		
		nTrialCorrectByCondition	= cellfun(@numel,cKTrialByConditionCorrect);
		nTrialUseMax				= min(nTrialCorrectByCondition);
		
		%find out how many trials we can use so that the lowest RT trials of
		%the highest mean RT condition match the highest RT trials of the lowest
		%mean RT condition
			cRTByCondition		= cellfun(@(kt) rt(kt),cKTrialByConditionCorrect,'uni',false);
			rtMeanByCondition	= cellfun(@mean,cRTByCondition);
			
			kConditionRTMin		= find(rtMeanByCondition==min(rtMeanByCondition));
			kConditionRTMax		= find(rtMeanByCondition==max(rtMeanByCondition));
			
			[cRTByConditionSort,cKSort]		= cellfun(@sort,cRTByCondition,'uni',false);
			cKTrialByConditionCorrectSort	= cellfun(@(kt,ks) kt(ks),cKTrialByConditionCorrect,cKSort,'uni',false);
			
			rtDiff		= arrayfun(@(n) abs(mean(cRTByConditionSort{kConditionRTMin}(end-n+1:end)) - mean(cRTByConditionSort{kConditionRTMax}(1:n))),(1:nTrialUseMax)');
			nTrialUse	= find(rtDiff==min(rtDiff));
		
		%fairly dumb way of choosing the trials to use
			rtMeanTarget	= mean([mean(cRTByConditionSort{kConditionRTMin}(end-nTrialUse+1:end)) mean(cRTByConditionSort{kConditionRTMax}(1:nTrialUse))]);
			
			rtMeanDiff		= cellfun(@(rt) arrayfun(@(k) abs(mean(rt(k:k+nTrialUse-1))-rtMeanTarget),(1:numel(rt)-nTrialUse+1)'),cRTByConditionSort,'uni',false);
			kMeanClosest	= cellfun(@(rtm) find(rtm==min(rtm),1),rtMeanDiff);
			
			cKTrialUse	= cellfun(@(kt,kstart) sort(kt(kstart:kstart+nTrialUse-1)),cKTrialByConditionCorrectSort,num2cell(kMeanClosest),'uni',false);
		
		%construct the target and chunk arrays
			kTrialUseAll	= cat(1,cKTrialUse{:});
			
			bUse				= zeros(size(squeeze(ifo.correct(kS,:,:))));
			bUse(kTrialUseAll)	= 1;
			
			[cTarget,cEvent]	= deal(cell(nRun,1));
			for kR=1:nRun
				block	= squeeze(ifo.(strScheme)(kS,kR,:));
				bUseRun	= reshape(bUse(kR,:),[],1);
				blockCI	= block + nCondition*bUseRun;
				
				cTarget{kR}	= block2target(blockCI,durBlock,durRest,cConditionCI,durPre,durPost,...
								'hrf'			, sLabel.HRF			, ...
								'block_offset'	, sLabel.BlockOffset	, ...
								'block_sub'		, sLabel.BlockSub		  ...
								);
				
				cEvent{kR}	= block2event(blockCI,durBlock,durRest,durPre,durPost);
			end
			
			cTargetRT{kS}	= cat(1,cTarget{:});
			
			event		= eventcat(cEvent,durRun);
			nEvent		= size(event,1);
			durRunTotal	= durRun*nRun;
			
			event(:,1)	= 1:nEvent;
			event(:,2)	= event(:,2) + sLabel.HRF + sLabel.BlockOffset;
			event(:,3)	= sLabel.BlockSub;
			ev			= event2ev(event,durRunTotal);
			
			cChunkRT{kS}	= sum(ev.*repmat(1:nEvent,[durRunTotal 1]),2);
		
		nTrialRecord(kS)	= nTrialUse;
		rtTrialRecord{kS}	= cellfun(@(k) rt(k),cKTrialUse,'uni',false);
	end

%some descriptive statistics about the trials
	disp(sprintf('trials per subject: %.3f (SEM %.3f)',mean(nTrialRecord),stderr(nTrialRecord)));
	
	mRT					= cellfun(@(rt) cellfun(@mean,rt)',rtTrialRecord,'uni',false);
	mRT					= cat(1,mRT{:});
	[pRT,atRT,statRT]	= anova1(mRT,[],'off');
	disp(sprintf('matched RT difference: F(%d,%d)=%e, p=%e',atRT{2,3},atRT{3,3},atRT{2,5},atRT{2,6}));

%analyze
	%data paths
		[C,cPathData,cPathMask]	= GO.Analyze.PCA(...
									'subject'	, cSubject		, ...
									'mask'		, cMask			, ...
									'dim'		, dimPCAGC		, ...
									'cores'		, nCore			, ...
									'load'		, true			, ...
									'force'		, false			, ...
									'silent'	, false			  ...
									);
	
	res	= struct;
	
	cOutPrefix	= cellfun(@(s) [s '-' strScheme '-rtmatch'],cSubject,'uni',false);
	
	kRun	= reshape(repmat(1:nRun,[durRun 1]),[],1);
	
	res.(strScheme)	= MVPAClassify(cPathData,cTargetRT,cChunkRT,...
						'path_mask'			, cPathMask			, ...
						'mask_name'			, cMask				, ...
						'spatiotemporal'	, true				, ...
						'target_blank'		, 'Blank'			, ...
						'zscore'			, kRun				, ...
						'output_dir'		, strDirOut			, ...
						'output_prefix'		, cOutPrefix		, ...
						'cores'				, nCore				, ...
						'debug'				, 'all'				, ...
						'force'				, false				, ...
						'force_each'		, false				, ...
						'silent'			, false				  ...
						);
	
	conf	= GO.ConfusionModels;
	
	%do this manually since we changed the code in between gridop and mwlearn
		for kM=1:nMask
			strMask	= cMask{kM};
			
			cm		= res.shape.result.(strMask).allway.confusion;
			cCM		= squeeze(mat2cell(cm,4,4,ones(nSubject,1)));
			cCMJK	= jackknife(@(x) {nanmean(cat(3,x{:}),3)},cCM);
			
			[r,ccstat]	= cellfun(@(c) corrcoef2(reshape(conf{1},[],1),reshape(c,1,[])),cCMJK);
			
			stat.(strMask)	= restruct(ccstat);
			stat.(strMask)	= rmfield(stat.(strMask),{'tails','df','t','p','cutoff','m','b'});
			
			stat.(strMask).mr	= nanmean(stat.(strMask).r);
			stat.(strMask).ser	= nanstderrJK(stat.(strMask).r);
			stat.(strMask).mz	= nanmean(stat.(strMask).z);
			stat.(strMask).sez	= nanstderrJK(stat.(strMask).z);
			
			[h,p,ci,stats]	= ttestJK(stat.(strMask).z,0,0.05,'right');
			
			stat.(strMask).p	= p;
			stat.(strMask).t	= stats.tstat;
			stat.(strMask).df	= stats.df;
		end
		
		p			= structfun(@(msk) msk.p,stat);
		[~,pfdr]	= fdr(p,0.05);
	
	disp([cMask cellfun(@(m) stat.(m).mz,cMask,'uni',false) cellfun(@(m) stat.(m).t,cMask,'uni',false) cellfun(@(m) stat.(m).p,cMask,'uni',false) num2cell(pfdr)]);
	
%save the results
	strPathOut	= PathUnsplit(strDirOut,'result','mat');
	save(strPathOut,'res','stat');
