import pandas as pd
import numpy as np
from scipy import integrate
import glob
def get_meta(dir_path):
    meta_frame = pd.read_csv(filepath_or_buffer=dir_path+'meta.csv', sep=',',index_col=0,
                            skip_blank_lines=True, names=['item','value'],
                            usecols=[7, 8])
    return meta_frame

def get_calib_factors(meta_frame):
    N2_calibs = []
    H2_calibs = []
    TCD_HC_calibs = []
    FID_HC_calibs = []
    injection = []
    for i in range(14): # 14 sample loops in total
        injection.append(i+2)

        cellname_N2="Calibration factor TCD_N2_Inj{}"
        cellname_N2.format(i+2) #we start with injection no.2, counting from 1
        N2_calibs.append(float(meta_frame['value'][cellname_N2.format(i+2)]))

        cellname_H2="Calibration factor TCD_H2_Inj{}"
        cellname_H2.format(i+2)
        H2_calibs.append(float(meta_frame['value'][cellname_H2.format(i+2)]))

        cellname_HC_TCD="Calibration factor TCD_HC_Inj{}"
        cellname_HC_TCD.format(i+2)
        TCD_HC_calibs.append(float(meta_frame['value'][cellname_HC_TCD.format(i+2)]))

        cellname_HC="Calibration factor FID_Inj{}"
        cellname_HC.format(i+2)
        FID_HC_calibs.append(float(meta_frame['value'][cellname_HC.format(i+2)]))
    # converting the 3 lists to a dataframe with injection number as index. HC = Hydrocarbon
    calib_factors = pd.DataFrame({'N2': N2_calibs, 'H2': H2_calibs,'HC_TCD':TCD_HC_calibs, 'HC_FID': FID_HC_calibs}, index=injection)
    calib_factors.loc[1] = [10,10,10,10] #placeholder
    calib_factors = calib_factors.sort_index()
    return calib_factors


def integrate_TCD_chromatogram(chrom_path,compound_frame_TCD):
    chromatogram = pd.read_csv(filepath_or_buffer=chrom_path, sep='\t',index_col=0,skiprows=42,thousands=r',')
    chromatogram = chromatogram.drop(columns=['Step (s)'])
    integrals = []
    for i in compound_frame_TCD.index.to_list():
        chrom_bl = chromatogram - chromatogram.iloc[chromatogram.index.get_loc(compound_frame_TCD['baseline_point'].loc[i])] #baseline_correction at selected point
        x = chrom_bl.index[(chrom_bl.index > compound_frame_TCD['lower_bound'].loc[i]) & (chrom_bl.index < compound_frame_TCD['upper_bound'].loc[i])]
        y = chrom_bl['Value (mV)'].loc[x]
        integrals.append(integrate.trapezoid(y,x))
    return integrals, chromatogram

def integrate_TCD_chromatogram_new(chrom_path,compound_frame_TCD):
    #1. Reads in a TCD chromatogram saved as .csv.
    #2. For every compound, the value at a point specified in in compound_frame_TCD is susbtracted as baseline.
    #3. The peak is integrated within the specified boundries.
    #4. The the integrals and the baseline corrected chromatogram are returned
    chromatogram = pd.read_csv(filepath_or_buffer=chrom_path, sep='\t',index_col=0,skiprows=42,thousands=r',')
    chromatogram = chromatogram.drop(columns=['Step (s)'])
    integrals = []
    for i in compound_frame_TCD.index.to_list():
        chrom_bl = chromatogram - chromatogram.iloc[chromatogram.index.get_loc(compound_frame_TCD['baseline_point'].loc[i])]
        x = chrom_bl.index[(chrom_bl.index > compound_frame_TCD['lower_bound'].loc[i]) & (chrom_bl.index < compound_frame_TCD['upper_bound'].loc[i])]
        y = chrom_bl['Value (mV)'].loc[x]
        bl_slope = (y.iloc[-1] - y.iloc[0])/(x[-1] - x[0])
        y = y - (bl_slope*(x-x[0])+y.iloc[0])

        integrals.append(integrate.trapezoid(y,x))
    return integrals, chrom_bl




def get_TCD_flows(filepath, compound_frame_TCD):
    #Reads in the TCD chromatograms and the metafile, integrates the peaks according to the sample list and calculates the molar flow.
    
    meta_frame = get_meta(filepath)
    calib_factors = get_calib_factors(meta_frame)
    calib_factors.loc[1] = [100,100,100,100] # blank row to not mess with changing indeces all the time
    calib_factors = calib_factors.sort_index()

    molar_flow = 2.018104822/50*float(meta_frame['value']['N2 flow rate: [ml/min]']) #mmol/min. 2.018 is the molar amount of N2 in 50ml at 25C and 10**5 Pa
    init = float(meta_frame['value']['GC start time:'])  # start time of filling the loops (min)
    inj_time = float(meta_frame['value']['Time sampling loop: [min]']) # time of recording of 1 injection in min

    filelist_TCD = sorted(glob.glob(filepath+'*TCD.txt*'))
    filelist_TCD.sort(key=len)
    integral_list = []
    #getting the integrals
    for i in enumerate(filelist_TCD):
        integrals, chromatogram = integrate_TCD_chromatogram(filelist_TCD[i[0]],compound_frame_TCD) 
        integral_list.append(np.array(integrals))
    TCD_frame = pd.DataFrame(integral_list, columns=compound_frame_TCD.index.to_list(),index=[i[0]+1 for i in enumerate(filelist_TCD)])
    #first row to nan
    TCD_frame.loc[1] = np.nan

    # Determining the molar flow of nitrogen which is used as internal standard.
    flow_frame_TCD_N2 = pd.DataFrame(data=molar_flow/(np.array(TCD_frame['Nitrogen']) / np.array(calib_factors['N2'])/100), columns=['Nitrogen flow'], index=TCD_frame.index)
    # Determining H2 flow
    H2_flow = np.array(TCD_frame['Hydrogen']) / np.array(calib_factors['H2'])/100 * flow_frame_TCD_N2['Nitrogen flow']
    TCD_frame.drop(columns=['Nitrogen','Hydrogen'], inplace=True)

    #multiply all the columns with flow_frame_TCD_N2_2 to get the flow corrected for potential volume expansion
    TCD_frame = TCD_frame.multiply(flow_frame_TCD_N2['Nitrogen flow'], axis=0)
    #molar flows are determined 
    for col in TCD_frame.columns:
        TCD_frame[col] = TCD_frame[col] / np.array(calib_factors['HC_TCD'])/100
    TCD_frame['Hydrogen'] = H2_flow
    TCD_frame['N2_corr'] = flow_frame_TCD_N2['Nitrogen flow']
    TCD_frame.index = (TCD_frame.index-1)*inj_time+init
    return TCD_frame 




def process_chromatogram(chrom_path,peaklist_path,time):
    #1. reads in a raw chromatogram exported by chromeleon to .csv
    #2. Takes the minimum value of the chromatogram as baseline
    #3. Integrates peaks according to the peaklist
    #4. Integrates the total chromatogram for the total carbon flow
    chromatogram = pd.read_csv(filepath_or_buffer=chrom_path, sep='\t',index_col=0,skiprows=42,thousands=r',')
    chromatogram = chromatogram.drop(columns=['Step (s)'])

    baseline = chromatogram['Value (pA)'].min() # left as is to allow different baseline
    chromatogram['bl_substracted'] = chromatogram['Value (pA)'] - baseline
    peaklist = pd.read_csv(filepath_or_buffer=peaklist_path, sep=',',index_col=0)
    integrals = []
    labels = peaklist.index.tolist()
    for i in range(len(peaklist)):
        integrals.append(integrate.trapezoid(chromatogram['bl_substracted'][(chromatogram.index>peaklist['lower_bound'].iloc[i]) & (chromatogram.index<peaklist['upper_bound'].iloc[i])],chromatogram.index[(chromatogram.index>peaklist['lower_bound'].iloc[i]) & (chromatogram.index<peaklist['upper_bound'].iloc[i])]))
    # append an integral of the total chromatogram
    labels.append('total')
    integrals.append(integrate.trapezoid(chromatogram['bl_substracted'],chromatogram.index))
    integral_frame_pA = pd.DataFrame({time: integrals}, index=labels)
    return chromatogram, integral_frame_pA


def get_integral_frame(filelist_FID, Peaklist_path, inj_time_timelist):
    integral_list = []
    for i in range(len(filelist_FID)):
        chromatogram, integral = process_chromatogram(filelist_FID[i],Peaklist_path,inj_time_timelist[i])
        integral_list.append(integral)

    integral_frame_FID = pd.concat(integral_list, axis=1)
    return integral_frame_FID





def calc_molarCflow(integral_frame_FID_pA, calib_factors, flow_frame_TCD):
    '''dividing by calibration factor and adjusting for N2 flow'''
    integral_frame_FID_mol = integral_frame_FID_pA / np.array(calib_factors['HC_FID'])[1:]/100 * flow_frame_TCD['N2_corr']
    return integral_frame_FID_mol.fillna(0)


def get_indiv_integrals(integral_frame_FID_mol, inj_time_timelist): 
    '''
    determines the masses of the individual compounds asusming a molar mass of 14 g/mol per carbon atom
    '''
    integral = integrate.trapezoid(integral_frame_FID_mol,inj_time_timelist)
    integral_frame = pd.DataFrame(integral, columns=['molar'], index=integral_frame_FID_mol.index)
    integral_frame['mass'] = integral_frame['molar']*14/1000 
    return integral_frame

def get_TCD_masses(integral_frame_TCD_mol, inj_time_timelist):
    #converts molar flows to masses
    
    dfT=integral_frame_TCD_mol.T.fillna(0)
    integral = integrate.trapezoid(dfT,inj_time_timelist)
    integral_frame = pd.DataFrame(integral, columns=['molar_TCD'], index=dfT.index)
    integral_frame['mass_TCD'] = integral_frame['molar_TCD']*[44/3,42/3,2]/1000
    return integral_frame


def get_liquid_yield(meta_frame):
    return float(meta_frame['value']['Mass of liquid sample: [g]'])


### TGA Section

def get_coke_amount(homedir, cat_mass = 0.625):
    '''determines the amount of coke formed on the catalyst from the weightloss of the spent catalyst determined by TGA and the total mass of cat used'''
    try:
        TGA_dir = homedir + '/TGA_burnoff/'
        TGA_list = glob.glob(TGA_dir + '*.txt')
        TGA_file = TGA_list[0]
    except:
        TGA_dir = homedir + '/TGA_Burnoff/'
        TGA_list = glob.glob(TGA_dir + '*.txt')
        TGA_file = TGA_list[0]
    names = ['Blank', 'Time', 'Unsubtracted weight', 'Baseline weight',
        'Program Temp.','Sample Temp.', 'Sample Purge Flow',
        'Balance purge flow']
    frame_full = pd.read_table(TGA_file,sep='\t', skiprows=98,header=None, names=names, engine='python', skipfooter=45)
    startline = frame_full[frame_full['Sample Temp.'] >= 150].index[0]
    endline = frame_full[frame_full['Sample Temp.'] >= 800].index[0]
    frame = frame_full.iloc[(startline-40):endline,:]
    min_mass = min(frame['Unsubtracted weight'])
    max_mass = max(frame['Unsubtracted weight'])
    percent_loss = (max_mass-min_mass)/max_mass*100

    total_coke = percent_loss/100*max_mass*cat_mass/min_mass
    return total_coke




def flatten_export(masses_export,cat): 
    new = masses_export.drop(['molar_TCD','molar'],axis = 1)
    new.drop(['Propane','Propylene','total'], inplace=True)
    new.sum(axis=1).to_frame(name=cat).T
    return new.sum(axis=1).to_frame(name=cat).T