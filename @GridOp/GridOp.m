classdef GridOp < PTB.Object
% GridOp
%
% Description:	the gridop experiment object 
%
% Syntax: go = GridOp(<options>)
%
%			subfunctions:
%				Start(<options>):	start the object
%				End:				end the object
%               Prepare:            prepare necessary info
%               Run:                execute a gridop run
%
% In:
% 	<options>:
%       debug:		(0) the debug level
%
% Out: 
%
% Updated: 2013-08-23
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.

	%PUBLIC PROPERTIES---------------------------------------------------------%
	properties
		Experiment;
		
		%running reward total
			reward;
		%images
			arrow	= [];
			op		= {};
	end
	%PUBLIC PROPERTIES---------------------------------------------------------%
	
	
	%PRIVATE PROPERTIES--------------------------------------------------------%
	properties (SetAccess=private, GetAccess=private)
		argin;
	end
	%PRIVATE PROPERTIES--------------------------------------------------------%
	
	
	%PROPERTY GET/SET----------------------------------------------------------%
	methods
		
	end
	%PROPERTY GET/SET----------------------------------------------------------%
	
	
	%PUBLIC METHODS------------------------------------------------------------%
	methods
		function go = GridOp(varargin)
			go	= go@PTB.Object([],'gridop');
			
			go.argin	= varargin;
			
			%parse the inputs
			opt = ParseArgs(varargin,...
				'input'			, []	, ...
				'session'		, []	, ...
				'disable_key'	, false	, ...
                'fullscreen'    , []    , ...
				'debug'			, 0 	  ...
				);
			if isempty(opt.session)
				opt.session	= conditional(opt.debug==2,1,2);
			end
			
			opt.name	= 'gridop';
			opt.context	= switch2(opt.session,1,'psychophysics',2,'fmri');
			opt.tr		= GO.Param('time','tr');
			
			%window
				opt.background	= GO.Param('color','back');
				opt.text_color	= GO.Param('color','text');
				opt.text_size	= GO.Param('text','size');
				opt.text_family	= GO.Param('text','font');
			
			opt.input_scheme	= 'lr';
			
			%options for PTB.Experiment object
			cOpt = opt2cell(opt);
			
			%initialize the experiment
			go.Experiment	= PTB.Experiment(cOpt{:});
			
			%set the session
			go.Experiment.Info.Set('go','session',opt.session);
			
			%start
			go.Start;
		end
		%----------------------------------------------------------------------%
		function Start(go,varargin)
		%start the gridop object
			go.argin	= append(go.argin,varargin);
			
			if ~notfalse(go.Experiment.Info.Get('go','prepared'))
				%prepare info
				go.Prepare(varargin{:});
			end
		end
		%----------------------------------------------------------------------%
		function End(go,varargin)
		%end the gridop object
			v	= varargin;
            
			go.Experiment.End(v{:});
		end
		%----------------------------------------------------------------------%
	end
	%PUBLIC METHODS------------------------------------------------------------%
end
