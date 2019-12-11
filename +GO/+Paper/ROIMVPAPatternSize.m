function [s,cMask] = ROIMVPAPatternSize(ifo)
% GO.Paper.ROIMVPAPatternSize
% 
% Description:	calculate the size of patterns used in the ROIMVPA analysis
% 
% Syntax:	[s,cMask] = GO.Paper.ROIMVPAPatternSize(ifo)
% 
% Updated: 2015-04-13
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirBase

strDirPCA	= DirAppend(strDirBase,'data','store','pca-auto_1');

cSubject	= ifo.code.fmri;
nSubject	= numel(cSubject);
cMask		= [GO.UnionMasks().core];
nMask		= numel(cMask);

s	= NaN(nSubject,nMask);

progress('action','init','total',nSubject,'name','subject','label','processing subjects');
for kS=1:nSubject
	strSubject	= cSubject{kS};
	
	progress('action','init','total',nMask,'name','mask','label','processing masks');
	for kM=1:nMask
		strMask	= cMask{kM};
		
		strPathPCA	= PathUnsplit(strDirPCA,sprintf('%s-%s',strSubject,strMask),'nii.gz');
		
		msk	= double(NIfTI.Read(strPathPCA,'return','data'));
		
		s(kS,kM)	= sum(msk);
		
		progress('name','mask');
	end
	
	progress('name','subject');
end
