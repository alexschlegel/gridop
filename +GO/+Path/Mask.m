function [cPathMask,cMask] = Mask(varargin)
% GO.Path.Mask
% 
% Description:	get the paths to mask files
% 
% Syntax:	[cPathMask,cMask] = GO.Path.Mask(<options>)
% 
% In:
% 	<options>:
%		subject:	(<all>) the subjects to include
%		mask:		(<core>) the mask names, or an nSubject x 1 cell of mask
%					paths
%		state:		('preprocess') see GO.SubjectInfo
% 
% Out:
% 	cPathMask	- an nSubject x 1 cell of nMask x 1 cells of mask paths
%	cMask		- an nMask x 1 cell of mask names
% 
% Updated: 2014-04-20
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirData;

opt	= ParseArgs(varargin,...
		'subject'	, {}			, ...
		'mask'		, {}			, ...
		'state'		, 'preprocess'	  ...
		);

cSubject	= GO.Subject('subject',opt.subject,'state',opt.state);
nSubject	= numel(cSubject);

if iscell(opt.mask) && numel(opt.mask)==nSubject
	cPathMask	= opt.mask;
	cMask		= cellfun(@(c) cellfun(@(m) getfield(regexp(PathGetFilePre(m,'favor','nii.gz'),'-(?<name>.+)$','names'),'name'),c,'uni',false),cPathMask,'uni',false);
	return;
end

if isempty(opt.mask)
	sMask		= GO.UnionMasks;
	opt.mask	= sMask.core;
end

strDirMask	= DirAppend(strDirData,'mask');

cMask		= opt.mask;
cPathMask	= cellfun(@(s) cellfun(@(m) PathUnsplit(DirAppend(strDirMask,s),m,'nii.gz'),cMask,'uni',false),cSubject,'uni',false);
