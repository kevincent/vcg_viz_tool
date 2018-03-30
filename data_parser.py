import os
import glob as gb
import pandas as pd
import numpy as np
import re

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def parse(initial_path = "/Users/kevin/Continuity/Simulations/MachineLearningPaper/LBBB_LVdyssync", vcg_which = False, verbose = False, vcg_select = '*', patient_number=0):

    """
    Get all data from patient number 'patient_number' in folder LBBB_LVdyssync stored into the 'patient'variable
	Returns:
	--------
	patient: 
        pandas.Series containing all data of the patient. See comments in __main__ for indexing.
        
        
    Warning: python3 is __next__, python2 is next
    Warning: in folder patient 1 you need to copy the sync file out of the vcg folder, otherwise it will be missing
    Credit: Yngve Moe wrote an original version of this script
    """


    vcg_type = 'dVCG' if vcg_which else 'kVCG'

    patients_folders = [os.path.join(initial_path, sub_path) for sub_path in os.walk(initial_path).next()[1]] 
    #patients = [None]*len(patients_folders)

    # for patient_folder in patients_folders[patient_number:]:
    for patient_folder in [patients_folders[patient_number]]:

        # include patient identifier
        pt_id = patient_folder.rsplit('/',1)[1]

        # load the optimal dysynchrony
        opt_desync = pd.read_csv(gb.glob(patient_folder + '/*sync_opt*')[0], header = None)

        # load the dysynchrony times of each simulation
        desync = pd.read_csv(gb.glob(patient_folder + '/*sync.*')[0], header = None)

	   # load the parameters of each simulation
        eval_values = pd.read_csv(gb.glob(patient_folder + '/*eval*')[0],
                                  header = None, names=['x','y','z','cond','L','R'],delim_whitespace=True, index_col=False)
        # search for nodes file
        if os.path.isfile(patient_folder+'/'+pt_id+'_nodes.txt'):
            nodes = np.loadtxt(patient_folder+'/'+pt_id+'_nodes.txt',skiprows=1,usecols=(0,8,16))
            #find unique pacing locations
            eval_values['pacing_nodes'] = zip(eval_values['x'],eval_values['y'],eval_values['z'])
            unique_nodes = eval_values['pacing_nodes'].unique()
            #match unique nodes to their node numbers
            node_list = []
            for n in unique_nodes:
                node_list.append(np.linalg.norm(np.asarray(n)-nodes,axis=1).argmin()+1)
            
            #map nodes back into dataframe
            xmap = dict(zip(eval_values['pacing_nodes'].unique(),node_list))
            #make these int?
            eval_values['pacing_nodes'] = eval_values['pacing_nodes'].map(xmap)
        
        # load the experimental VCG
        vcg_real = pd.read_csv(gb.glob(patient_folder + '/*measured*/*' + vcg_type + '*')[0],
                               header = None, names=['VCGx','VCGy','VCGz'], sep=' ', index_col=False)    

        # load the VCG of each simulation
        vcg_reading = gb.glob(patient_folder + '/*model*/' + vcg_select + '.txt')
        vcg_reading = sorted(vcg_reading, key=numericalSort)
        vcg_sims = [int(filename.split('_')[-2]) for filename in vcg_reading]
        vcg_model = [np.nan]*max(vcg_sims) #may be better to fill DataFrame with nan or zeros for access during plotting
        for i in range(len(vcg_reading)):
            vcg_model[vcg_sims[i]-1] = pd.read_csv(vcg_reading[i], header = None, names=['VCGx','VCGy','VCGz'], sep='\t', index_col=False)
        
        #Print some info about data
        if False:
            print(pt_id + '\n'+ ''.join(str(node_list)))
        if verbose:
            print(pt_id +
                  '\n\tSimulations in pts_to_eval: %i \
                  \n\tVCG files: %i \
                  \n\tMax VCG number: %i \
                  \n\tNumber of pacing nodes %i'
                  %(len(eval_values),len(vcg_reading),max(vcg_sims),len(unique_nodes)))

        # store all into a panda series
        patient = pd.Series([pt_id, opt_desync[0], desync[0], eval_values, vcg_real,
                             vcg_model], index=['pt_id', 'opt_desync', 'desync', 'eval_values', 'vcg_real', 'vcg_model'])

    return patient

if __name__ == "__main__":

    patient = parse(initial_path = "/Users/kevin/Continuity/Simulations/MachineLearningPaper/LBBB_LVdyssync", vcg_which = False, patient_number=0)

    #mean_plot(patient)
            
    # example calls
    print patient['opt_desync'][0]
    #print patient['desync']
    #print patient['eval_values']
    #print patient['vcg_real'] # ['px'],['py'],['pz']
    #print patient['vcg_model'] # [i] # ['px'],['py'],['pz']


