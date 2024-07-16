import pandas as pd
import datetime
from scipy import integrate
import os


def get_water_content(TGA_file):
    names = ['Blank', 'Time', 'Unsubtracted weight', 'Baseline weight',
        'Program Temp.','Sample Temp.', 'Sample Purge Flow',
        'Balance purge flow']
    frame_full = pd.read_table(TGA_file,sep='\t', skiprows=38,header=None, names=names, engine='python', skipfooter=45)
    endline = frame_full[frame_full['Sample Temp.'] >= 800].index[0]
    frame = frame_full.iloc[0:endline,:]
    min_mass = min(frame['Unsubtracted weight'])
    max_mass = max(frame['Unsubtracted weight'])
    water_content = (max_mass-min_mass)/max_mass
    return water_content



def read_chromatogram(Path):
    #Reads the ascii file of the injection. Splits in into a metadate frame and a frame containing chromatograms of the 3 channels.
    #Sampling frequency must be equal for all channels.

    df_chromatogram = pd.read_csv(Path, sep='\t', header=None,skiprows=13)
    df_chromatogram_meta = pd.read_csv(Path, sep='\t', header=None,skiprows=0,nrows=13)

    #split the chromatogram into the different columns. This is achieved by splitting the frame at indexes matching 1/3 and 2/3 of the lenght
    df_Channel_1 = df_chromatogram.iloc[0:int(len(df_chromatogram)/3)].reset_index(drop=True)
    df_Channel_2 = df_chromatogram.iloc[int(len(df_chromatogram)/3):int(2*len(df_chromatogram)/3)].reset_index(drop=True)
    df_Channel_3 = df_chromatogram.iloc[int(2*len(df_chromatogram)/3):len(df_chromatogram)].reset_index(drop=True)

    #combine the 3 Channels into one frame, new index
    df_chromatogram = pd.concat([df_Channel_1,df_Channel_2,df_Channel_3],axis=1)
    sampling_freqs = [float(df_chromatogram_meta.iloc[7][0].split(',')[i]) for i in range(1,4)]
    if len(set(sampling_freqs)) != 1:
        raise ValueError('The sampling frequencies are not equal')

    sampling_freq = 1/sampling_freqs[0] #Hz

    df_chromatogram['Time[s]'] = df_chromatogram.index*sampling_freq
    #set time as index
    df_chromatogram = df_chromatogram.set_index('Time[s]')
    df_chromatogram.columns = ['FID_L','FID_M','TCD']
    return df_chromatogram, df_chromatogram_meta




def baseline_correct(frame):
    #for each channel, substract the minimum value
    for i in range(0,3):
        #substract the mean of the values between 288 and 294 seconds
        frame.iloc[:,0] = frame.iloc[:,0] - frame.iloc[288*50:294*50,0].mean()
        frame.iloc[:,1] = frame.iloc[:,1] - frame.iloc[288*50:294*50,1].mean()
        frame.iloc[:,2] = frame.iloc[:,2] - frame.iloc[10*50:15*50,2].mean()
        # frame.iloc[:,i] = frame.iloc[:,i] - frame.iloc[:,i].min()
    return frame




def integrate_peaks(frame, Peaks_channel2,Peaks_TCD,Peaks_channel1):
    df = frame['FID_M'] #only look into second channel
    integrals_list = []
    for i in range(0,len(Peaks_channel2)):
        #ideantify peaks and valleys in the range of the peak
        integral = integrate.trapezoid(df[(df.index > Peaks_channel2[i][1][0]) & (df.index < Peaks_channel2[i][1][1])])
        integrals_list.append(integral)
    df_TCD = frame['TCD']
    for i in range(0,len(Peaks_TCD)):
        integral = integrate.trapezoid(df_TCD[(df_TCD.index > Peaks_TCD[i][1][0]) & (df_TCD.index < Peaks_TCD[i][1][1])])
        integrals_list.append(integral)

    df_FIDL = frame['FID_L']
    for i in range(0,len(Peaks_channel1)):
        integral = integrate.trapezoid(df_FIDL[(df_FIDL.index > Peaks_channel1[i][1][0]) & (df_FIDL.index < Peaks_channel1[i][1][1])])
        integrals_list.append(integral)


    return integrals_list





def parse_log(logfile_path):
    Log = pd.read_csv(logfile_path, sep='\t', skiprows=1)
    Log = Log[['Date','Time','MFC 1 pv','MFC 2 pv','MFC 3 pv','MFC 4 pv', 'Oven Temperature','v11-reactor','v10-bubbler','v12-gc']]
    # for the v-11 reactor columns, replace all 0 with 'reactor' else 'bypass'
    Log['v11-reactor'] = Log['v11-reactor'].apply(lambda x: 'reactor' if x == 0 else 'bypass')
    Log.rename(columns={'MFC 1 pv':'N2_flow'}, inplace=True)
    Log.rename(columns={'MFC 4 pv':'He_Bubbler'}, inplace=True)
    Log.rename(columns={'MFC 3 pv':'He_Dilution'}, inplace=True)
    Log['timestamp'] = pd.to_datetime(Log['Date'] + ' ' + Log['Time'],format='%m/%d/%Y %I:%M:%S %p')
    Log['timestamp'] = Log.apply(lambda row: datetime.datetime.strptime(row['Date'] + ' ' + row['Time'], '%m/%d/%Y %I:%M:%S %p'), axis=1)
    return Log


def get_peaknames(Peaks_channel1,Peaks_channel2,Peaks_TCD):
    Peaknames=[Peaks_channel2[i][0] for i in range(0,len(Peaks_channel2))]
    Peaknames_TCD=[Peaks_TCD[i][0] for i in range(0,len(Peaks_TCD))]
    for i in range(0,len(Peaknames_TCD)):
        Peaknames.append(Peaknames_TCD[i])
    for i in range(0,len(Peaks_channel1)):
        Peaknames.append(Peaks_channel1[i][0])
    return Peaknames


def get_chrom_paths(data_dir,chromatogram_list):
    Paths = []
    for i in range(0,(len(chromatogram_list))):
        try:
            Paths.append(os.path.normpath(data_dir+'/chromatograms/' + chromatogram_list[i]))
        except:
            Paths.append(os.path.normpath(data_dir+'/Chromatograms/' + chromatogram_list[i]))
    return Paths

def process_chromatograms(Paths,Peaknames,Peaks_channel2,Peaks_TCD,Peaks_channel1):
    '''
    Reads in the chromatogram, baseline corrects it and integrates the peaks specified by the peak lists.
    Returns a frame with the peak integrals and injection timestamp
    '''
    Integral_Frame = pd.DataFrame(columns=Peaknames)

    for i in range(0,len(Paths)):
        df_chromatogram, df_chromatogram_meta = read_chromatogram(Paths[i])
        df_chromatogram = baseline_correct(df_chromatogram) # Baseline is the mean of a flat reagion in the chromatogram
        result_frame = integrate_peaks(df_chromatogram, Peaks_channel2,Peaks_TCD,Peaks_channel1)
        result_frame.append(pd.to_datetime(df_chromatogram_meta.iloc[6][0].split(',')[1]))
        Integral_Frame.loc[i] = result_frame
    return Integral_Frame


 
def get_temp_and_valves(Integral_Frame,Log):
    '''
    Gets the temperature and bubbler/bypass valves from the logfile.
    At the timestamp of the injection, look in the log file for the position of the valves and the temperature
    '''
    Integral_Frame['Temperature'] = Integral_Frame['Timestamp'].apply(lambda x: Log['Oven Temperature'][Log['timestamp'].sub(x).abs().idxmin()])
    Integral_Frame['v10-bubbler'] = Integral_Frame['Timestamp'].apply(lambda x: Log['v10-bubbler'][Log['timestamp'].sub(x).abs().idxmin()])
    Integral_Frame['v11-reactor'] = Integral_Frame['Timestamp'].apply(lambda x: Log['v11-reactor'][Log['timestamp'].sub(x).abs().idxmin()])
    return Integral_Frame





def get_mean_DMP(Arrhenius_Frame, Temperature):
    #list of indices where the temperature matches
    idx = Arrhenius_Frame[Arrhenius_Frame['Temperature'] == Temperature].index
    DMP_leftover = Arrhenius_Frame['DMP_norm'][idx[2:-2]].mean()
    return DMP_leftover

def get_mean_product(Arrhenius_Frame, Temperature):
    #list of indices where the temperature matches
    idx = Arrhenius_Frame[Arrhenius_Frame['Temperature'] == Temperature].index
    product_yield = Arrhenius_Frame['Products_norm'][idx[2:-2]].mean()
    return product_yield





