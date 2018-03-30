import numpy as np
import pandas as pd
from copy import deepcopy

from numpy.linalg import eig, inv, norm
from numpy import mat
from sklearn.decomposition import PCA
from scipy import interpolate, ndimage


def project_PCA(patient_vcg, by_coords=False, plotting=False):
    """Transform the VCG signal  with PCA. The new x, y, z axis are ordered
    by decreasing variance.
    Returns:
    --------
    by_coords = False
     vcg_projected : 
        Numpy arrays of same sizes as the vcg signal. With transformed coordinates.
    
    by_coords = True
    x(s), y(s), z(s) : np.ndarray, np.ndarray, np.ndarray
        Numpy arrays of same sizes as the vcg signal. x, y, z coordinates
        for the transformed points. In this way it is easier to project on the 2D plane
        by grabbing the first two components.
        
    Credit: Yngve Moe wrote an original version of this script
    """    
    #http://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_3d.html
    
    vcg = patient_vcg.as_matrix()
    center(vcg)
    pca = PCA(n_components=3)
    pca.fit(vcg)
    vcg_projected = pca.transform(vcg)
    
    if plotting:
        plt.scatter(vcg_projected[:,0], vcg_projected[:,1])
        plt.show()

    if by_coords: return vcg_projected[:,0], vcg_projected[:,1], vcg_projected[:,2]
    return vcg_projected


def center(vcg):
    vcg -= vcg.mean(0, keepdims=True) 


def normalize_patient(patient):
    """Normalises the VCG by dividing by the maximum heart vector.
        Creates a new 
    """
#    mag = np.sqrt(x.dot(x)) is much faster, consider refactoring
#   https://stackoverflow.com/questions/9171158/how-do-you-get-the-magnitude-of-a-vector-in-numpy

    patient['vcg_real_norm'] = patient['vcg_real']  / norm(patient['vcg_real'],axis=1).max()
    patient['vcg_model_norm'] = [None]*len(patient['vcg_model'])
    for i, vcg_model in enumerate(patient['vcg_model']):

        # Standardise VCG signal:
        if isinstance(patient['vcg_model'][i], pd.DataFrame):
            patient['vcg_model_norm'][i] = patient['vcg_model'][i] / norm(patient['vcg_model'][i],axis=1).max()
        else:
            patient['vcg_model_norm'][i] = np.nan

    return patient


def center_patient(patient):
    """Centers the heart vector for each patient such that its average is 0.
    This requires the heart vector to be in cartesian coordinates!
    """

    patient = deepcopy(patient)
    patient['vcg_model'] = deepcopy(patient['vcg_model'])

    patient['vcg_real'] -= patient['vcg_real'].mean()

    for i, vcg_model in enumerate(patient['vcg_model']):
        patient['vcg_model'][i] -= vcg_model.mean()

    return patient


def resample_patient(patient, length=None):
    """Resamples the heart vector for each patient by velocity so that it is uniformly sampled in space (instead of in time).
    This requires the heart vector to be in cartesian coordinates!
    """
    
    patient = deepcopy(patient)
    patient['vcg_model'] = deepcopy(patient['vcg_model'])

    for i, vcg in enumerate(patient['vcg_model']):
        resampled_vcg = np.array(resample_by_velocity(vcg, length=length)).T
        patient['vcg_model'][i] = pd.DataFrame(resampled_vcg, columns=vcg.columns)
        #patient['vcg_model'][i]['VCGx'], patient['vcg_model'][i]['VCGy'], patient['vcg_model'][i]['VCGz'] = resample_by_velocity(patient['vcg_model'][i], length=length)

    resampled_vcg = np.array(resample_by_velocity(vcg, length=length)).T
    patient['vcg_real'] = pd.DataFrame(resampled_vcg, columns = patient['vcg_real'].columns)
    #patient['vcg_real']['VCGx'], patient['vcg_real']['VCGy'], patient['vcg_real']['VCGz'] = resample_by_velocity(patient['vcg_real'], length=length)

    return patient


def project_patient(patient):
    """Uses PCA on the heart vector and transforms it so that the z axis contains the least variance.
    This requires the heart vector to be in cartesian coordinates!
    """
    
    patient = deepcopy(patient)
    patient['vcg_model'] = deepcopy(patient['vcg_model'])
    
    for i in range(len(patient['vcg_model'])):
        patient['vcg_model'][i]['VCGx'], patient['vcg_model'][i]['VCGy'], patient['vcg_model'][i]['VCGz'] = project_PCA(patient['vcg_model'][i], by_coords=True)   

    patient['vcg_real']['VCGx'], patient['vcg_real']['VCGy'], patient['vcg_real']['VCGz'] = project_PCA(patient['vcg_real'], by_coords=True)

    return patient

def create_data_matrix(patient, transforms=None):
    """Returns a datamatrix for the patients where each column is a simulation.
    Arguments
    ---------
    patient : pandas.DataFrame
        The dataframe containing all information about the patient
    transforms : Array like
        List containing three transformations to be applied to the different
        axis of the VCG signal.
    
    Returns
    -------
    data_matrix : np.ndarray
        A matrix where each column is a VCG signal for each simulation (the
        three dimensions are just concatenated to one vector).
    """
    data_matrix = np.zeros((len(patient['vcg_model']), np.prod(patient['vcg_model'][0].shape)))
    for i, simulation in enumerate(patient['vcg_model']):
        sim = simulation.values.copy()

        if transforms is not None:
            for j in range(3):
                if transforms[j] is not None:
                    sim[j, :] = transforms[j](sim[j, :])

        data_matrix[i] = sim.reshape(np.prod(sim.shape))

    return data_matrix
