function p = Param(varargin)
% GO.Param
% 
% Description:	get a gridop parameter
% 
% Syntax:	p = GO.Param(f1,...,fN)
% 
% In:
% 	fK	- the the Kth parameter field
% 
% Out:
% 	p	- the parameter value
%
% Example:
%	p = GO.Param('color','back');
% 
% Updated: 2013-08-30
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global SIZE_MULTIPLIER;
persistent P;

if isempty(SIZE_MULTIPLIER)
	SIZE_MULTIPLIER	= 1;
end


if isempty(P)
	%stimulus parameters
		P.color	= struct(...
					'back'	, [128 128 128]	, ...
					'fore'	, [0 0 0]		, ...
					'text'	, [0 0 0]		  ...
					);
		P.shape	= struct(...
					'rect'	, {{[0 1 1 1; 0 1 0 0; 0 1 0 0; 1 1 0 0],[1 1 1 1; 0 1 0 0; 0 1 0 0; 0 1 0 0]}}	, ...
					'polar'	, {{[1 1 0 1; 1 1 0 1; 0 1 0 0; 0 1 0 0],[0 0 1 0; 1 1 1 0; 1 1 0 0; 1 1 0 0]}}	, ...
					...%'rect'	, {{[1 0 0 0; 1 1 1 1; 0 0 0 1; 0 0 0 1],[1 1 1 1; 0 1 0 0; 0 1 0 0; 0 1 0 0]}}	, ...
					...%'polar'	, {{[1 1 0 1; 1 1 0 1; 0 1 0 0; 0 1 0 0],[0 1 0 0; 0 1 1 1; 0 1 1 1; 0 1 1 1]}}	, ...
					'polar_exp'	, 1		, ... %makes each ring a different radius
					'rect_f'	, 0.94	  ... %make the rect stimuli take up a smaller space, to mac
					);
		P.size	= struct(...
					'stim'		, 200					, ...
					'stimva'	, 2.5*SIZE_MULTIPLIER	, ...
					'offset'	, 2*SIZE_MULTIPLIER		  ... %offset of test screen, in d.v.a.
					);
	%timing
		P.time	= struct(...
					'tr'		, 2000	, ...
					'prompt'	, 1		, ...
					'operation'	, 3		, ...
					'test'		, 1		, ...
					'feedback'	, 1		, ...
					'rest'		, 4		, ...
					'prepost'	, 4		  ...
					);
	%experiment design
		P.exp	= struct(...
					'runs'	, 15	, ...
					'reps'	, 1		  ...
					);
	%text
		P.text	= struct(...
					'font'	, 'Helvetica'			, ...
					'size'	, 0.75*SIZE_MULTIPLIER	  ...
					);
	%reward
		P.reward	= struct(...
						'base'		, 20	, ...
						'max'		, 40	, ...
						'penalty'	, 5		  ... %penalty is <- times the reward
						);
	%prompt mapping
		P.prompt	= struct(...
						'stimulus'		, 'ABCD'				, ... 
						'operation'		, '1234'				, ... 
						'distance'		, 1*SIZE_MULTIPLIER		, ... %distance from fixation to prompt
						'arrow'			, 0.5*SIZE_MULTIPLIER	  ...
						);
	%response buttons
		P.response	= struct(...
						'correct'	, 'left'	, ...
						'incorrect'	, 'right'	  ...
						);
	%scheme
		P.scheme	= {'shape';'operation'};
	%stuff for transfer entropy calculations
		P.te		= struct(...
						'block'	, 5	  ...	%number of samples to use per trial
						);
end

p	= P;

for k=1:nargin
	v	= varargin{k};
	
	switch class(v)
		case 'char'
			switch v
				case 'sizemultiplier'
					p	= SIZE_MULTIPLIER;
				case 'trialperrun'
					p	= 4*4*P.exp.reps;
				case 'trtrial'
					p	= P.time.prompt + P.time.operation + P.time.test;
				case 'trrest'
					p	= P.time.feedback + P.time.rest;
				case 'trrestpre'
					p	= P.time.prepost;
				case 'trrestpost'
					p	= P.time.feedback + P.time.prepost;
				case 'trrun'
					nTrial	= GO.Param('trialperrun');
					
					p	= GO.Param('trrestpre') + nTrial*GO.Param('trtrial') + (nTrial-1)*GO.Param('trrest') + GO.Param('trrestpost');
				case 'trtotal'
					p	= P.exp.runs*GO.Param('trrun');
				case 'trun'
					p	= GO.Param('trrun')*GO.Param('time','tr')/1000/60;
				case 'ttotal'
					p	= GO.Param('trtotal')*GO.Param('time','tr')/1000/60;
				case 'trialtotal'
					p	= GO.Param('trialperrun')*P.exp.runs;
				case 'rewardpertrial'
					p	= (P.reward.max - P.reward.base)/GO.Param('trialtotal');
				case 'penaltypertrial'
					p	= GO.Param('rewardpertrial')*P.reward.penalty;
				otherwise
					if isfield(p,v)
						p	= p.(v);
					else
						p	= [];
						return
					end
			end
		otherwise
			if iscell(p)
				p	= p{v};
			else
				p	= [];
				return;
			end
	end
end
