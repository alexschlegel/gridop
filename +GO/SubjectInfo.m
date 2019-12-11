function s = SubjectInfo(varargin)
% GO.SubjectInfo
% 
% Description:	compile a struct of subject info
% 
% Syntax:	s = GO.SubjectInfo(<options>)
% 
% In:
% 	<options>:
%		subject:		(<all>) a cell of subject ids to include
%		exclude:		(<none>) a cell of subject ids to exclude
%		state:			('preprocess') the state of subjects to return. one of
%						the following:
%							all:		all subjects
%							fmri:		subjects with fmri sessions
%							preprocess:	preprocessed subjects
%		shuffle:		(false) true to shuffle condition labels
%		shuffle_seed:	(19811011) the seed for the label shuffling
% 
% Updated: 2014-07-26
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirBase strDirData;

warning('off','MATLAB:indeterminateFields');

opt	= ParseArgs(varargin,...
		'subject'		, {}			, ...
		'exclude'		, {}			, ...
		'state'			, 'preprocess'	, ...
		'shuffle'		, false			, ...
		'shuffle_seed'	, 19811011		  ...
		);

opt.exclude	= ForceCell(opt.exclude);
opt.state	= CheckInput(opt.state,'state',{'all','fmri','preprocess'});

status(sprintf('selected subject state: %s',opt.state));

%condition information
	s.condition	= struct(...
					'shape'		, {{'R1';'R2';'P1';'P2'}}	, ...
					'operation'	, {{'CW';'CCW';'H';'V'}}	  ...
					);
	
	nShape		= numel(s.condition.shape);
	nOperation	= numel(s.condition.operation);
	shpRep		= repmat(s.condition.shape,[1 nOperation]);
	opRep		= repmat(s.condition.operation',[nShape 1]);
	s.condition.shapeop	= cellfun(@(s,o) [s '_' o],shpRep,opRep,'uni',false);
	
	cScheme	= fieldnames(s.condition);
	nScheme	= numel(cScheme);
%get the subject ids
	if isempty(opt.subject)
		cPathSubject	= FindFiles(strDirData,'^\w\w\w?\.mat$');
		
		cID	= cellfun(@PathGetFilePre,cPathSubject,'uni',false);
	else
		cID	= reshape(cellfun(@(id) regexprep(id,'\d{2}\w{3}\d{2}(\w{2,3})$','$1'),ForceCell(opt.subject),'uni',false),[],1);
	end
	
	%exclude
		cID	= setdiff(cID,opt.exclude);
	
	nSubject	= numel(cID);

	s.id	= cID;
%get some subject info
	s.subject		= dealstruct('age','gender','handedness',NaN(nSubject,1));
	s.subject.map	= dealstruct('shape','operation',NaN(nSubject,4));
	
	for kS=1:nSubject
		strPathSubject	= PathUnsplit(strDirData,cID{kS},'mat');
		x				= load(strPathSubject);
		
		s.subject.gender(kS)		= switch2(x.ifoSubject.gender,'f',0,'m',1,NaN);
		s.subject.handedness(kS)	= switch2(x.ifoSubject.handedness,'r',0,'l',1,NaN);
		
		s.subject.map.shape(kS,:)		= x.ifoSubject.map_stim;
		s.subject.map.operation(kS,:)	= x.ifoSubject.map_op';
	end

%get the practice and fmri session paths
	[s.code.practice,s.code.fmri,s.code.eeg,cPathPractice,cPathFMRI,cPathEEG]	= deal(cell(nSubject,1));
	
	for kS=1:nSubject
		cPathSubject		= FindFiles(strDirData,['\d' s.id{kS} '\.mat$']);
		nSession			= numel(cPathSubject);
		tSession			= cellfun(@(f) ParseSessionCode(PathGetFilePre(f)),cPathSubject);
		[tSessionS,kSort]	= sort(tSession);
		
		if nSession>0
			cPathPractice{kS}	= cPathSubject{kSort(1)};
			s.code.practice{kS}	= PathGetFilePre(cPathPractice{kS});
			
			if nSession>1
				cPathFMRI{kS}	= cPathSubject{kSort(2)};
				s.code.fmri{kS}	= PathGetFilePre(cPathFMRI{kS});
				
				if nSession>2
					cPathEEG{kS}	= cPathSubject{kSort(3)};
					s.code.eeg{kS}	= PathGetFilePre(cPathEEG{kS});
					
					if nSession>3
						error([num2str(nSession) ' sessions found for subject ' s.id{kS} '.']);
					end
				end
			end
		end
	end

%is the subject data preprocessed? the presence of the occ mask is a good indication.
	s.subject.preprocess	= cellfun(@(s) FileExists(PathUnsplit(DirAppend(strDirData,'mask',s),'occ','nii.gz')),s.code.fmri);

%read the practice data (+ subject age)
	s.practice.history	= cell(nSubject,1);
	s.practice.ntrial	= NaN(nSubject,1);
	
	for kS=1:nSubject
		if ~isempty(cPathPractice{kS})
			x	= load(cPathPractice{kS});
			
			s.subject.age(kS)	= ConvertUnit(x.PTBIFO.experiment.start - x.PTBIFO.subject.dob,'ms','day')/365.25;
			
			if isfield(x.PTBIFO.go,'practice_record')
				s.practice.history{kS}	= logical(x.PTBIFO.go.practice_record);
				s.practice.ntrial(kS)	= numel(x.PTBIFO.go.practice_record);
			else
				status(['no practice record for ' s.code.practice{kS} '!'],'warning',true);
				s.practice.history{kS}	= logical([]);
				s.practice.ntrial(kS)	= 0;
			end
		end
	end
%read the fmri data
	nRun	= GO.Param('exp','runs');
	nTrial	= GO.Param('trialperrun');
	tTest	= GO.Param('time','prompt') + GO.Param('time','operation');
	TR		= GO.Param('time','tr');
	
	rBase	= GO.Param('reward','base');
	dRight	= GO.Param('rewardpertrial');
	dWrong	= GO.Param('penaltypertrial');
	
	[s.shape,s.operation,s.correct,s.rt]	= deal(NaN(nSubject,nRun,nTrial));
	s.reward								= NaN(nSubject,1);
	
	for kS=1:nSubject
		if ~isempty(cPathFMRI{kS})
			x	= load(cPathFMRI{kS});
			
			for kR=1:nRun
				if ~isempty(x.PTBIFO.go.result{kR})
					s.shape(kS,kR,:)		= s.subject.map.shape(kS,[x.PTBIFO.go.result{kR}.input]);
					s.operation(kS,kR,:)	= s.subject.map.operation(kS,[x.PTBIFO.go.result{kR}.operation]);
					s.correct(kS,kR,:)		= reshape([x.PTBIFO.go.result{kR}.correct],[],1);
					
					%response time
						tResponse	= NaN(nTrial,1);
						for kT=1:nTrial
							if ~isempty(x.PTBIFO.go.result{kR}(kT).tresponse)
								tResponse(kT)	= x.PTBIFO.go.result{kR}(kT).tresponse{end};
							end
						end
						
						s.rt(kS,kR,:)	= TR*(tResponse - tTest);
				end
			end
			
			%reward
				bCorrect	= reshape(s.correct(kS,:,:),[],1);
				bCorrect	= logical(bCorrect(~isnan(bCorrect)));
				
				s.reward(kS)	= rBase;
				
				nTrialTotal	= numel(bCorrect);
				for kT=1:nTrialTotal
					if bCorrect(kT)
						s.reward(kS)	= s.reward(kS) + dRight;
					else
						s.reward(kS)	= max(rBase,s.reward(kS) - dWrong);
					end
				end 
		end
	end
	
	%should we shuffle the labels?
		if opt.shuffle
			strm	= RandStream.create('mt19937ar','seed',opt.shuffle_seed);
			RandStream.setGlobalStream(strm);
			
			s.shape		= randomize(s.shape,3,'seed',floor(rand*1000000));
			s.operation	= randomize(s.operation,3,'seed',floor(rand*1000000));
		end
	
	%combined shape+operation (corresponds to indices in s.condition.shapeop)
		s.shapeop	= s.shape + nShape*(s.operation-1);
%differences between conditions?
	s.mean	= structfun2(@(x) dealstruct('correct','rt',NaN(nSubject,numel(x))),s.condition);
	
	for kS=1:nSubject
		bCorrect	= squeeze(s.correct(kS,:,:));
		rt			= squeeze(s.rt(kS,:,:));
		
		for kC=1:nScheme
			strScheme	= cScheme{kC};
			nCondition	= numel(s.condition.(strScheme));
			
			kCondition		= squeeze(s.(strScheme)(kS,:,:));
			
			s.mean.(strScheme).correct(kS,:)	= arrayfun(@(k) nanmean(bCorrect(kCondition==k)),(1:nCondition)');
			s.mean.(strScheme).rt(kS,:)			= arrayfun(@(k) nanmean(rt(kCondition==k)),(1:nCondition)');
			s.mean.(strScheme).rtcorrect(kS,:)	= arrayfun(@(k) nanmean(rt(kCondition==k & bCorrect==1)),(1:nCondition)');
		end
	end

%construct the target and chunk arrays
	sTargetBlank	= dealstruct('all','correct',cell(nSubject,1));
	sLabelBlank		= struct(...
						'chunk'		, sTargetBlank									, ...
						'target'	, structfun2(@(x) sTargetBlank, s.condition)	  ...
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
	
	for kL=1:nLabelScheme
		strLabel	= cLabelScheme{kL};
		sLabel		= sLabelScheme.(strLabel);
		
		s.label.(strLabel)	= sLabelBlank;
		
		for kS=1:nSubject
			if any(isnan(reshape(s.shape(kS,:,:),[],1)))
				continue;
			end
			
			for kC=1:nScheme
				strScheme	= cScheme{kC};
				
				cCondition		= reshape(s.condition.(strScheme),[],1);
				nCondition		= numel(cCondition);
				cConditionCI	= [repmat({'Blank'},[nCondition 1]); cCondition];
				
				%all
					[cTarget,cEvent]	= deal(cell(nRun,1));
					for kR=1:nRun
						block		= squeeze(s.(strScheme)(kS,kR,:));
						cTarget{kR}	= block2target(block,durBlock,durRest,cCondition,durPre,durPost,...
										'hrf'			, sLabel.HRF			, ...
										'block_offset'	, sLabel.BlockOffset	, ...
										'block_sub'		, sLabel.BlockSub		  ...
										);
						
						if kC==1
							cEvent{kR}	= block2event(block,durBlock,durRest,durPre,durPost);
						end
					end
					
					s.label.(strLabel).target.(strScheme).all{kS}	= cat(1,cTarget{:});
					
					if kC==1
						event		= eventcat(cEvent,durRun);
						nEvent		= size(event,1);
						durRunTotal	= durRun*nRun;
						
						event(:,1)	= 1:nEvent;
						event(:,2)	= event(:,2) + sLabel.HRF + sLabel.BlockOffset;
						event(:,3)	= sLabel.BlockSub;
						ev			= event2ev(event,durRunTotal);
						
						s.label.(strLabel).chunk.all{kS}	= sum(ev.*repmat(1:nEvent,[durRunTotal 1]),2);
					end
					
				%just correct
					[cTarget,cEvent]	= deal(cell(nRun,1));
					for kR=1:nRun
						block		= squeeze(s.(strScheme)(kS,kR,:));
						bCorrect	= squeeze(s.correct(kS,kR,:));
						blockCI		= block + nCondition*bCorrect;
						
						cTarget{kR}	= block2target(blockCI,durBlock,durRest,cConditionCI,durPre,durPost,...
										'hrf'			, sLabel.HRF			, ...
										'block_offset'	, sLabel.BlockOffset	, ...
										'block_sub'		, sLabel.BlockSub		  ...
										);
						
						if kC==1
							cEvent{kR}	= block2event(blockCI,durBlock,durRest,durPre,durPost);
						end
					end
					
					s.label.(strLabel).target.(strScheme).correct{kS}	= cat(1,cTarget{:});
					
					if kC==1
						event							= eventcat(cEvent,durRun);
						event(event(:,1)<=nCondition,:)	= [];
						
						nEvent		= size(event,1);
						durRunTotal	= durRun*nRun;
						
						event(:,1)	= 1:nEvent;
						event(:,2)	= event(:,2) + sLabel.HRF + sLabel.BlockOffset;
						event(:,3)	= sLabel.BlockSub;
						ev			= event2ev(event,durRunTotal);
						
						s.label.(strLabel).chunk.correct{kS}	= sum(ev.*repmat(1:nEvent,[durRunTotal 1]),2);
					end
			end
		end
	end

%only keep the specified subjects
	switch opt.state
		case 'all'
		%nothing to do
			return;
		case 'fmri'
			bKeep	= ~cellfun(@isempty,s.code.fmri);
		case 'preprocess'
			bKeep	= s.subject.preprocess;
	end
	
	s	= KeepSubject(s,bKeep);

%------------------------------------------------------------------------------%
function s = KeepSubject(s,b)
	sKeep			= structtreefun(@(x) KeepOne(x,b),s);
	sKeep.condition	= s.condition;
	s				= sKeep;
end
%------------------------------------------------------------------------------%
function x = KeepOne(x,b)
	if size(x,1)==nSubject
		cSub	= subsall(x);
		cSub{1}	= b;
		x		= x(cSub{:});
	end
end
%------------------------------------------------------------------------------%

end
