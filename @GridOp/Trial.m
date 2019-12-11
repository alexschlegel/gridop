function res = Trial(go,kRun,kTrial,varargin)
% GridOp.Trial
% 
% Description:	run a GridOp trial
% 
% Syntax:	res = go.Trial(kRun,kTrial,[tStart]=<now>,<options>)
% 
% In:
%	kRun	- the run number
% 	kTrial	- the trial number
%	tStart	- the start time, in TRs. if unspecified, the trial starts
%			  immediately and uses PTB.Now time to advance
%	<options:
%		prompttexture:		(<none>) the name of a texture holding the prompt
%		promptshape:		(<none>) if the prompt was rendered already, the
%							shape location array
%		promptoperation:	(<none>) if the prompt was rendered already, the
%							operation location array
% 
% Out:
% 	res	- a struct of results
% 
% Updated: 2013-09-24
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
[tStart,opt]	= ParseArgs(varargin,[],...
					'prompttexture'		, []	, ...
					'promptshape'		, []	, ...
					'promptoperation'	, []	  ...
					);

bPractice	= isempty(tStart);
bShowPrompt	= isempty(opt.prompttexture);

res	= struct;

%get the mapping between button indices and test choices
	kButtonCorrect		= cell2mat(go.Experiment.Input.Get('correct'));
	kButtonIncorrect	= cell2mat(go.Experiment.Input.Get('incorrect'));
%get the stimulus and operation mappings for the current subject
	mapStim	= go.Experiment.Subject.Get('map_stim');
	mapOp	= go.Experiment.Subject.Get('map_op');

%get/generate the trial parameters
	trial	= go.Experiment.Info.Get('go','trial');
	
	%target input stimulus and operation
		res.input		= trial.input(kRun,kTrial);
		res.operation	= trial.op(kRun,kTrial);
	%target stimulus and operation location in the prompt
		res.input_location	= trial.prompt.location(kRun,kTrial);
	%output stimulus and operation
		res.test_correct		= trial.test.correct(kRun,kTrial);
		res.output_operation	= conditional(res.test_correct,...
									res.operation							, ...
									randFrom(1:4,'exclude',res.operation)	  ...
									);
	%test stimuli and operations
		res.output_location		= trial.test.location(kRun,kTrial);
		kDistractor				= setdiff(1:4,res.output_location);
		
		[res.test_shape,res.test_operation]	= deal(NaN(4,1));
		
		res.test_shape(res.output_location)	= res.input;
		res.test_shape(kDistractor)			= randomize(setdiff(1:4,res.input));
		
		res.test_operation(res.output_location)	= res.output_operation;
		res.test_operation(kDistractor)			= randomize(setdiff(1:4,res.output_operation));
%set up the sequence
	t	= GO.Param('time');
	
	cSequence	=	{
						@ShowPrompt
						@ShowOperation
						@ShowTest
					};
	tSequence	=	cumsum([
						t.prompt
						t.operation
						t.test
					]);
	fWait			=	{
							@Wait_Default
							@Wait_Default
							@Wait_Response
						};
	
	if bPractice
		tSequence	= tSequence*t.tr;
		tStart		= PTB.Now;
		strTUnit	= 'ms';
	else
		strTUnit	= 'tr';
	end
%run the sequence
	bFirstWait	= true;
	
	if bPractice
		go.Experiment.Scanner.StartScan;
	end
	
	kLastResponse	= NaN;
	[tStart,tEnd,tShow,bAbort,kResponse,tResponse]	= go.Experiment.Show.Sequence(cSequence,tSequence,...
															'tunit'			, strTUnit		, ...
															'tstart'		, tStart		, ...
															'tbase'			, 'sequence'	, ...
															'fwait'			, fWait			, ...
															'fixation'		, false			  ...
															);
	
	if bPractice
		go.Experiment.Show.Blank('fixation',true);
		go.Experiment.Window.Flip;
		
		go.Experiment.Scanner.StopScan;
	end
	
	res.tstart		= tStart;
	res.tend		= tEnd;
	res.tshow		= tShow;
	res.abort		= bAbort;
	res.kresponse	= kResponse;
	res.tresponse	= tResponse;
	
	if ~isempty(res.kresponse)
		res.correct	= 	(res.test_correct  && ismember(res.kresponse{end},kButtonCorrect)) || ...
						(~res.test_correct && ismember(res.kresponse{end},kButtonIncorrect));
	else
		res.correct	= false;
	end

%------------------------------------------------------------------------------%
function ShowPrompt(varargin)
	if bShowPrompt
	%construct the prompt now
		[res.prompt_shape,res.prompt_operation]	= go.ShowPrompt(kRun,kTrial,varargin{:});
	else
	%someone else already did it, just transfer the texture
		res.prompt_shape		= opt.promptshape;
		res.prompt_operation	= opt.promptoperation;
		
		go.Experiment.Show.Texture(opt.prompttexture,varargin{:});
	end
end
%------------------------------------------------------------------------------%
function ShowOperation(varargin)
	go.Experiment.Show.Blank('fixation',true,varargin{:});
end
%------------------------------------------------------------------------------%
function ShowTest(varargin)
	%show each test stimulus
		dOffset	= GO.Param('size','offset');
		sStim	= GO.Param('size','stimva');
		xStim	= dOffset*[-1 1 -1 1];
		yStim	= dOffset*[-1 -1 1 1];
		
		for kS=1:4
			[rot,flip]	= GO.Operate(res.test_operation(kS),'map',mapOp);
			
			im	= GO.Stim.Stimulus(res.test_shape(kS),...
					'map'		, mapStim	, ...
					'rotation'	, rot		, ...
					'flip'		, flip		  ...
					);
			
			pStim	= [xStim(kS) yStim(kS)];
			
			go.Experiment.Show.Image(im,pStim,sStim,varargin{:});
		end
	%show the fixation dot
		go.Experiment.Show.Fixation(varargin{:});
end
%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%
function [bAbort,kResponse,tResponse] = Wait_Default(tNow,tNext)
	bAbort		= false;
	kResponse	= [];
	tResponse	= [];
	
	if bFirstWait
	%wait longer during the first call, so log events, etc. can post
		go.Experiment.Scheduler.Wait(PTB.Scheduler.PRIORITY_LOW);
		
		bFirstWait	= false;
	else
		go.Experiment.Scheduler.Wait(PTB.Scheduler.PRIORITY_CRITICAL);
	end
end
%------------------------------------------------------------------------------%
function [bAbort,kResponse,tResponse] = Wait_Response(tNow,tNext)
	bAbort							= false;
	[dummy,dummy,dummy,kResponse]	= go.Experiment.Input.DownOnce('response');
	
	tResponse	= conditional(isempty(kResponse),[],tNow);
	
	go.Experiment.Scheduler.Wait(PTB.Scheduler.PRIORITY_CRITICAL);
end
%------------------------------------------------------------------------------%

end
