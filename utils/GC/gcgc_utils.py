
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import scipy.integrate as integrate
import matplotlib.colors




def parse_chromatogram(path)->pd.DataFrame:
    '''Parses a FID chromatoram from a .csv file into a dataframe with the colymns 'Ret.Time[s]' and 'Absolute Intensity'''
    df = pd.read_csv(path, sep = ',',  skiprows = 1)
    df.columns = ['Time(ms)', 'Time(min)', 'unknown','Absolute Intensity']
    df['Ret.Time'] = df['Time(min)']
    df['Ret.Time']-=df['Ret.Time'][0]
    df['Ret.Time[s]']=df['Ret.Time']*60
    df.drop(['Time(ms)','Time(min)','unknown'], axis=1, inplace=True)
    return df


def split_solvent(df: pd.DataFrame, solvent_time: float)->pd.DataFrame:
    '''sets absolute intensity to 0 when time is < solvent_time'''
    df.loc[df['Ret.Time[s]'] <= solvent_time, 'Absolute Intensity'] = 0
    return df

def add_split(df,split_time,sampling_interval):#split time in s, sampling interval in ms
    '''Adds a new column to the Dataframe that indicates the number of the modulation. e.g., all rows with no. 10 will belong to the same injection on the second column.
    The function will add rows to the dataframe if the last modulation is not complete.'''
    rows_splitting = split_time/sampling_interval*1000
    df['split_no_fromindex'] = df.index//rows_splitting 
    while len(df[df['split_no_fromindex']==df['split_no_fromindex'].max()]) < len(df[df['split_no_fromindex']==df['split_no_fromindex'].min()]):
        df.loc[len(df)+1] = df.iloc[-1]
    return df

def min_correct(df): # 'global minimum' as baseline
    df['Absolute Intensity'] = df['Absolute Intensity'] - df['Absolute Intensity'].min()
    return df

def baseline_stridewise(df_array): #minimum in each stride as baseline
    df_array = df_array - df_array.min(axis=0)
    return df_array

def convert_to2D(df,split_time):
    df_short = df[['split_no_fromindex','Absolute Intensity']]
    array_list = []
    
    for i in range(0,int(df_short['split_no_fromindex'].max())):
        array_list.append(df_short[df_short['split_no_fromindex']==i]['Absolute Intensity'].values)

    #turn arraylist into an 2D array
    array = np.zeros((len(array_list),len(array_list[0])))
    for i in range(0,len(array_list)):
        array[i,:] = array_list[i]

    index_list_retention_time = []
    for i in range(len(array_list)):
        index_list_retention_time.append(i*split_time)

    columns_splittime = []
    for i in range(len(array_list[0])):
        columns_splittime.append(round(i*split_time/len(array_list[0]),3))

    df_array = pd.DataFrame(array, index = index_list_retention_time, columns = columns_splittime)
    df_array= df_array.T
    df_array = df_array.iloc[::-1]
    return df_array

def shift_phase(df_array, shift):
    '''shifts the phase of the 2D chromatogram'''
    df_array_shifted = np.roll(df_array, shift, axis=0)
    return df_array_shifted

import tifffile





def normalize_array(df_array):
    '''Normalizes the Volume of a 2D chromatogram to 1'''
    arrray_for_integral = np.array(df_array)
    integral_non_norm = integrate.trapz(integrate.trapz(arrray_for_integral, axis=0), axis=0)
    df_array_norm = df_array/integral_non_norm
    return df_array_norm
    
def integrate_masked(df_norm_array, maskpath):
    mask =  tifffile.imread(maskpath)/255 # if the mask is binary no need to divide by 255
    df_norm_array_masked  = df_norm_array*mask
    array_norm_mask_diarom = np.array(df_norm_array_masked)
    row_integrated = integrate.trapz(array_norm_mask_diarom, axis=0)
    column_integrated = integrate.trapz(row_integrated, axis=0)
    return column_integrated
    


def mask_integrate(df_array_norm, mask_dir) -> pd.DataFrame:
    '''
    Returns a frame consisting of a the integrals for each mask.
    Masks are a 2D array of the same size as the 2D chromatogram. The mask is a binary image where the areas of interest is 255 and the rest is 0.
    The chromatogram is multiplied by each mask individually, followed by an integration.
    To adjust, change the masks in the mask_dir folder.
    '''
    mask_list = glob.glob(mask_dir + '*.tif') #
    mask_names = []
    for i in range(len(mask_list)):
        mask_name = mask_list[i].split('\\')[-1].split('.')[0]
        mask_names.append(mask_name.split('Mask_')[-1])

    integral_list = []
    for i in range(len(mask_list)):
        integral_list.append(integrate_masked(df_array_norm, mask_list[i]))

    df_integral = pd.DataFrame([integral_list], columns = mask_names)
    df_integral['unassigned'] = 1-df_integral[mask_names].sum(axis=1) # the remaining chromatogram volume is left as unassigned.
    return df_integral

def process_chromatogram(filepath, split_time: float, sampling_interval: float, mask_dir,shift: float =0, solvent_time: float =0)->pd.DataFrame:
    df = parse_chromatogram(filepath)
    df = split_solvent(df, solvent_time)
    df = add_split(df,split_time,sampling_interval)
    df_array = convert_to2D(df,split_time)
    df_array = baseline_stridewise(df_array)
    df_array = shift_phase(df_array, shift)
    df_array_norm = normalize_array(df_array)
    df_integral = mask_integrate(df_array_norm, mask_dir)
    return df_integral, df_array_norm



def plot_2Dchromatogram(df_array_norm,maskdir,savedir,plotmask=True,title = '2D Chromatogram',split_time = 20):
    mask_list = glob.glob(maskdir + '*.tif')
    plt.figure(figsize=(8,8/1.615))
    plt.imshow(np.sqrt(df_array_norm), cmap='viridis', interpolation='nearest', extent=[0, 106, 0, split_time], aspect='auto')
    plt.colorbar(label='$\sqrt{\mathrm{intensity}}$')

    colormaplist = [matplotlib.colors.ListedColormap(['none', 'C'+str(i)]) for i in range(6) ]
    annotations = [['Alkanes/Alkenes',28,7],['Monoaromatics',45,12],['Diaromatics',30,19],['Triaromatics',64,4.5],['Pyrenes',93,8]]
    if plotmask:
        for i in range(len(mask_list)):
            mask =  tifffile.imread(mask_list[i])/255 # if the mask is binary no need to divide by 255
            plt.imshow(mask, cmap=colormaplist[i], interpolation='nearest', extent=[0, 106, 0, split_time], aspect='auto', alpha=0.2)
            plt.text(annotations[i][1], annotations[i][2], annotations[i][0], fontsize=8, color='white')

    plt.xlabel('Retention time 1 (min)')
    plt.ylabel('Retention time 2 (s)')
    plt.title(title)
    plt.savefig(savedir, dpi=1000, bbox_inches='tight')

