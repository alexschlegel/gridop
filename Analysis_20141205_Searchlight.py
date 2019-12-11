print 'preparing script'

import os
import json
from mvpa2.suite import *

#__file__ = '/mnt/tsestudies/wertheimer/gridop/code/Analysis_20141205_Searchlight.py'

# enable debug output for searchlight call
if __debug__:
	debug.active += ["SLC"]

#parameters
param_path = os.path.splitext(os.path.realpath(__file__))[0] + '.json'
with open(param_path,'r') as f:
	param = json.load(f)

num_subject = len(param['subject'])

#classifier
clf = LinearCSVMC()

#partitioner
partitioner = NFoldPartitioner(cvtype=1)
#we're taking all targets, not just correct ones, so no balancing needed
#balancer = Balancer(attr='targets', count=10, limit='partitions', apply_selection=True)
#partitioner = ChainNode([partitioner,balancer],space='partitions')

#cross validator
cv_kwargs = {
		'errorfx':   lambda p,t: np.mean(p == t),
		'enable_ca': ['stats'],
	}
cv = CrossValidation(clf, partitioner, **cv_kwargs)

#searchlight
sl = sphere_searchlight(cv, radius=3, space='voxel_indices', postproc=mean_sample())

#perform each classification
for idx_subject in range(0,num_subject):
	subject = str(param['subject'][idx_subject])
	
	#get some initial attributes
	scheme = param['scheme'][0]
	attr_path = param[scheme]['attr'][idx_subject]
	attr = SampleAttributes(attr_path)
	
	#load the data
	print 'loading data for %s' % (subject)
	data_path = param['data_path'][idx_subject]
	dataset = fmri_dataset(
		samples=data_path,
		targets=attr.targets,
		chunks=attr.chunks,
		mask=param['gm_path'],
		)
	
	#zscore by run
	zscore(dataset, chunks_attr='chunks', param_est=('targets', ['Blank']), dtype='float32')
	
	#perform each searchlight
	for scheme in param['scheme']:
		print 'searchlight: %s / %s' % (subject,scheme)
		
		#load the scheme parameters
		attr_path = param[scheme]['attr'][idx_subject]
		attr = SampleAttributes(attr_path)
		
		#copy a stripped down version of the data and assign targets
		ds = dataset.copy(deep=False,
							sa=['targets', 'chunks'],
							fa=['voxel_indices'],
							a=['mapper'])
		ds.sa['targets'] = attr.targets
		
		#average samples within each chunk
		avg = mean_group_sample(['targets','chunks'])
		ds  = ds.get_mapped(avg)
		
		#remove blanks
		ds = ds[ds.sa.targets != 'Blank']
		
		#perform the searchlight
		sl_map = sl(ds)
		
		#save the results
		output_path = os.path.join(param['output_dir'],'%s_%s.nii' % (subject,scheme))
		nii = map2nifti(sl_map, imghdr=dataset.a.imghdr)
		nii.to_filename(output_path)
