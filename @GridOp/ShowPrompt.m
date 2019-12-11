function [shp,op] = ShowPrompt(go,kRun,kTrial,varargin)
% GridOp.ShowPrompt
% 
% Description:	show the trial prompt screen
% 
% Syntax:	go.ShowPrompt(kRun,kTrial,<options>)
% 
% In:
%	kRun	- the run number
% 	kTrial	- the trial number
%	<options>:
%		window:	('main') the name of the window on which to show the prompt
% 
% Out:
%	shp	- the shape display order
%	op	- the operation display order
% 
% Updated: 2013-09-25
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'window'	, 'main'	  ...
		);

%get the shape and operation order
	trial	= go.Experiment.Info.Get('go','trial');
	
	kInput		= trial.input(kRun,kTrial);
	kOperation	= trial.op(kRun,kTrial);
	
	locInput	= trial.prompt.location(kRun,kTrial);
	kDistractor	= setdiff(1:4,locInput);
	
	[shp,op]		= deal(NaN(4,1));
	
	shp(locInput)		= kInput;
	shp(kDistractor)	= randomize(setdiff(1:4,kInput));
	
	op(locInput)	= kOperation;
	op(kDistractor)	= randomize(setdiff(1:4,kOperation));
%blank the screen
	go.Experiment.Show.Blank('fixation',false,'window',opt.window);
%show the prompts
	mStimulus	= GO.Param('prompt','stimulus');
	mOperation	= GO.Param('prompt','operation');
	
	chrPrompt		= mStimulus(shp);
	chrOperation	= mOperation(op);

	dPrompt	= GO.Param('prompt','distance');
	
	xPrompt	= dPrompt*[-1 0 1 0];
	yPrompt	= dPrompt*[0 -1 0 1] + 0.25;
	
	strSize	= num2str(GO.Param('text','size'));
	
	for kP=1:4
		go.Experiment.Show.Text(['<size:' strSize '>' chrPrompt(kP) chrOperation(kP) '</size>'],[xPrompt(kP) yPrompt(kP)],'window',opt.window);
	end
%show the arrow
	im	= imrotate(go.arrow,(1-locInput)*90);
	
	go.Experiment.Show.Image(im,[],GO.Param('prompt','arrow'),'window',opt.window);
