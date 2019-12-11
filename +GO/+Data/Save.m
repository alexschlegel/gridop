function Save(d,strName,varargin)
% GO.Data.Save
% 
% Description:	save data for later retrieval
% 
% Syntax:	GO.Data.Save(d,strName,[param]=<none>)
% 
% In:
% 	d		- the data
%	strName	- the data name
%	[param]	- a variable storing parameters for the data, in case different
%			  versions of the data exist with different parameters 
% 
% Updated: 2014-03-03
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirData;

param	= ParseArgs(varargin,[]);

save(GO.Data.Path(strName,param),'d','param');
