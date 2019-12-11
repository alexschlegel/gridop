% Analysis_20140503_GCMVPA_Control
% control of the GC MVPA analysis, with split halves of the data. hopefully we
% get qualitatively similar results, but if not it could just be a power issue
nCore			= 12;
strNameAnalysis	= '20140503_gcmvpa_control';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;
cSubject	= ifo.code.fmri;
nSubject	= numel(cSubject);

cMask		= GO.UnionMasks().core;
nMask		= numel(cMask);

cMaskLabel		= upper(cMask);
cMaskLabel{5}	= 'LOC';

dimICAGC		= 10;

colMask	=	[
				0	222	222	%dlpfc
				204	0	204	%fef
				255	0	0	%occ
				179	75	0	%pcu
				0	51	204	%loc
				0	96	0	%ppc
			]/255;

%split the subjects in half randomly
	kOrder	= randomize(1:nSubject)';
	kSplit	= floor(nSubject/2);
	kGroup1	= kOrder(1:kSplit);
	kGroup2	= kOrder(kSplit+1:end);
	
	cSubjectGroup	= {cSubject(kGroup1); cSubject(kGroup2)};
	nGroup			= numel(cSubjectGroup);

for kG=1:nGroup
	cSubjectCur	= cSubjectGroup{kG};
	
	ifo			= GO.SubjectInfo('subject',cSubjectCur);
	
	%GC Connectivity Pattern Classification!
		[res,stat] = GO.Analyze.GCMVPA(...
						'subject'	, cSubjectCur	, ...
						'mask'		, cMask			, ...
						'dim'		, dimICAGC		, ...
						'ifo'		, ifo			, ...
						'cores'		, nCore			  ...
						);
		
		%what did we get?
			cScheme	= stat.label{1};
			nScheme	= numel(cScheme);
			
			[conn,pConn,pFDRConn]	= deal(NaN(nMask,nMask,nScheme));
			
			for kS=1:nScheme
				strScheme	= stat.label{1}{kS};
				
				disp(sprintf('Group %d, %s:',kG,strScheme));
				
				cMaskSrc	= stat.label{2};
				nMaskSrc	= numel(cMaskSrc);
				
				for kMS=1:nMaskSrc
					strMaskSrc	= cMaskSrc{kMS};
					
					cMaskDst		= cMaskSrc;
					bSelf			= strcmp(strMaskSrc,cMaskDst);
					cMaskDst(bSelf)	= [];
					nMaskDst		= numel(cMaskDst);
					
					for kMD=1:nMaskDst
						strMaskDst	= cMaskDst{kMD};
						
						r		= stat.confusion.corr.group.allway.r(kS,kMS,kMD);
						z		= fisherz(r);
						
						p		= stat.confusion.corr.group.allway.p(kS,kMS,kMD);
						pfdr	= stat.confusion.corr.group.allway.pfdr(kS,kMS,kMD);
						
						pacc	= stat.accuracy.p.allway(kS,kMS,kMD);
						paccfdr	= stat.accuracy.pfdr.allway(kS,kMS,kMD);
						
						kSrc	= find(strcmp(strMaskSrc,cMask));
						kDst	= find(strcmp(strMaskDst,cMask));
						
						conn(kSrc,kDst,kS)		= z;
						pConn(kSrc,kDst,kS)		= p;
						pFDRConn(kSrc,kDst,kS)	= pfdr;
						
						if pfdr <= sqrt(0.1)
							disp(sprintf('   %6s to %6s: r=%f, p=%f, pfdr=%f, pacc=%f, paccfdr=%f',strMaskSrc,strMaskDst,r,p,pfdr,pacc,paccfdr));
						end
					end
				end
			end
		
		%connectivity figures
			%custom anterior to posterior order
				cLabel	=		upper({
									'dlpfc'
									'fef'
									'pcu'
									'ppc'
									'loc'
									'occ'
								});
				
				cLabel		= cMaskLabel;
				[b,kOrder]	= ismember(cLabel,cMaskLabel);
				colLabel	= colMask(kOrder,:);
			
			col	= GetPlotColors(2);
			col	= col(end:-1:1,:);
			
			for kS=1:nScheme
				strScheme	= cScheme{kS};
				
				C		= ReorderConfusion(conn(:,:,kS),cMaskLabel,cLabel);
				p		= ReorderConfusion(pConn(:,:,kS),cMaskLabel,cLabel);
				pfdr	= ReorderConfusion(pFDRConn(:,:,kS),cMaskLabel,cLabel);
				
				%hierarchical network plot
					%get the line with bounds
						rAll	= stat.confusion.corr.group.allway.r;
						zAll	= fisherz(rAll);
						pAll	= stat.confusion.corr.group.allway.p;
						
						zShow	= zAll(pAll<=0.05);
						zMin	= min(zShow);
						zMax	= max(zShow);
						
						tLineMin	= 1;
						tLineMax	= 8;
					
					x			= C;
					x(p>0.05)	= 0;
					x(isnan(x))	= 0;
					
					b		= biograph(x,cLabel);
					
					n		= b.Nodes;
					nNode	= numel(n);
					for kN=1:nNode
						n(kN).Color		= colMask(kN,:);
						n(kN).LineColor	= colMask(kN,:);
						n(kN).TextColor	= GetGoodTextColor(n(kN).Color);
					end
					
					e		= b.Edges;
					nEdge	= numel(e);
					for kE=1:nEdge
						sNode	= regexp(e(kE).ID,'(?<src>.+) -> (?<dst>.+)','names');
						kSrc	= find(strcmp(cMask,sNode.src));
						kDst	= find(strcmp(cMask,sNode.dst));
						
						e(kE).LineWidth	= MapValue(e(kE).Weight,zMin,zMax,tLineMin,tLineMax);
						e(kE).LineColor	= col(kS,:);
					end
					
					g	= biograph.bggui(b);
					f	= get(g.biograph.hgAxes,'Parent');
					
					strPathOut	= PathUnsplit(strDirOut,sprintf('gcnetwork-%s-group%d',strScheme,kG),'png');
					saveas(f,PathAddSuffix(strPathOut,'','eps'));
					fig2png(f,strPathOut);
			end
end
