% get the MNI coordinates of the masks
strDirMask	= DirAppend(strDirData,'mni-mask');

cPathMask	= FindFiles(strDirMask,'(left|right)\.nii\.gz');
cMask		= cellfun(@(f) PathGetFilePre(f,'favor','nii.gz'),cPathMask,'uni',false);
nMask		= numel(cPathMask);

[m,p]	= cellfunprogress(@MRIMaskPosition,cPathMask,'uni',false);

for kM=1:nMask
	s		= regexp(cMask{kM},'(?<name>.+)-(?<hemi>.+)','names');
	strHemi	= switch2(s.hemi,'left','LH','right','RH');
	
	disp(sprintf('% 5s/%s: x=%6.2f, y=%6.2f, z=%6.2f',s.name,strHemi,m{kM}(1),m{kM}(2),m{kM}(3)));
end
