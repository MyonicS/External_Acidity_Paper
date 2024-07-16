import glob
import pandas as pd
import pytz
import spectrochempy as scp
import numpy as np
import xarray as xr
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


def parse_log(exp_path): #loading the log file and converting the time to a datetime object
    
    logpath = glob.glob(exp_path + 'log/*.txt')[0]
    log = pd.read_csv(logpath, sep='\t', header=1, names=['Date', 'Time', 'OvenSetpoint', 'OvenTemperature'],engine = 'python')
    log['DateTime'] = pd.to_datetime(log['Date'] + ' ' + log['Time'], format='%m/%d/%Y %I:%M:%S %p')
    #adjusting the timezone
    basename = os.path.basename(logpath)
    month, day = int(basename.split(' ')[0].split('-')[1]),int(basename.split(' ')[0].split('-')[2])
    if (month <=3 and day <= 31) or (month == 10 and day >= 29) or (month >10):
        offset = 60
    else:
        offset = 120

    log['DateTime'] = log['DateTime'].dt.tz_localize(pytz.FixedOffset(offset)) 
    return log




def get_timestamps(scp_ar): #extracting and converting timestamps from scp array
    time_array = scp_ar.y.values.magnitude # array of seconds
    timestamps = pd.to_timedelta(time_array, unit='s') + pd.Timestamp('1970-01-01')
    timestamps = timestamps.tz_localize('UTC').tz_convert('Europe/Amsterdam')
    return timestamps



def add_temp(scp_ar, log,timestamps):
    #adds the temperature as an extra coordinate in the scp xarray
    indices = np.searchsorted(log['DateTime'], timestamps, side='left')
    indices = np.clip(indices, 0, len(log['DateTime'])-1)
    temp_spec = np.array(log['OvenTemperature'][indices].reset_index(drop=True))
    temps = temp_spec
    temps = scp.Coord(temps, title="temperature", units='degree_Celsius')
    c_times=scp_ar.y.copy()
    scp_ar.y =[c_times, temps]
    return scp_ar


def xr_convert(scp_ar,background):
   spectrum = np.abs((-np.log10(scp_ar)+np.log10(background)).data)
   wavenumbers = scp_ar.x.values.magnitude
   timestamps = scp_ar.y['acquisition timestamp (GMT)'].values.magnitude
   df = pd.DataFrame(spectrum, index=timestamps, columns=wavenumbers)
   data_array = xr.DataArray(df, dims=["time", "wavenumber"], coords={"time": timestamps, "wavenumber": wavenumbers})
   data_array["time"] = pd.to_datetime(data_array.time.values, unit="s")
   data_array.coords["temperature"] = ("time", scp_ar.y.temperature.values.magnitude)
   return data_array   



def get_indices(exp_path,exp_name,save_indices=False,print_indices=True):
    if any('indices' in s for s in os.listdir(exp_path)):

        #load the indices file
        filen_name = glob.glob(exp_path + '*indices*')[0]

        indices = pd.read_csv(filen_name, skiprows=1, names=['start_bl', 'end_bl', 'start_dose', 'end_dose', 'start_desorb', 'end_desorb', 'start_dry', 'end_dry', '150_plateau'], engine='python')
        start_bl, end_bl, start_dose, end_dose, start_desorb, end_desorb, start_dry, end_dry, index_150_plateau = indices.iloc[0]
        if print_indices == True:
            print('indices file found')
            print(indices)
    else:
        print('no indices file found, please define the indices manually')
        #for Z11
        start_bl = 156
        end_bl = 280
        start_dose= 281
        end_dose = 370

        start_desorb = 371
        end_desorb = 542

        #extra: identify the drying region
        start_dry = 0
        end_dry = 155

        #extra:
        index_150_plateau = 420
    
        if save_indices == True:
            cutoff_idices = pd.DataFrame({'start_bl':start_bl,'end_bl':end_bl,'start_dose':start_dose,'end_dose':end_dose,'start_desorb':start_desorb,'end_desorb':end_desorb,'start_dry':start_dry,'end_dry':end_dry, '150_plateau':index_150_plateau}, index=[0])
            cutoff_idices.to_csv(exp_path+exp_name+'_cutoff_indices.csv',index=False)
    index_lib = {'start_bl':start_bl,'end_bl':end_bl,'start_dose':start_dose,'end_dose':end_dose,'start_desorb':start_desorb,'end_desorb':end_desorb,'start_dry':start_dry,'end_dry':end_dry, '150_plateau':index_150_plateau}
    return index_lib



def split_experiment(data_array, indices_lib):
    data_array_bl = data_array[indices_lib['start_bl']:indices_lib['end_bl']]
    data_array_dose = data_array[indices_lib['start_dose']:indices_lib['end_dose']]
    data_array_desorb = data_array[indices_lib['start_desorb']:indices_lib['end_desorb']]
    data_array_dry = data_array[indices_lib['start_dry']:indices_lib['end_dry']]
    return data_array_bl, data_array_dose, data_array_desorb, data_array_dry




def fit_integrate_peak(data, peak_loc,peak_window, fit_window, fit_select,plot =True,multi1=0,multi2=0,color='r'):
    #models
    def lorentzian(x, x0, gamma, A):
        return A * gamma**2 / ((x-x0)**2 + gamma**2)
    def gaussian(x, x0, sigma, A):
        return A * np.exp(-(x-x0)**2 / (2*sigma**2))
    
    if fit_select == 'lorentzian':
        fit_type = lorentzian
    elif fit_select == 'gaussian':
        fit_type = gaussian


    #isolate finding window
    data_window = data.sel(wavenumber=slice(peak_loc+peak_window, peak_loc-peak_window))
    peak_intensity = data_window.values.max()

    #fitting window
    fit_window = [data_window.idxmax().values-fit_window, data_window.idxmax().values+fit_window]
    data_window = data_window.sel(wavenumber=slice(fit_window[1], fit_window[0]))
    fit = curve_fit(fit_type, data_window['wavenumber'], data_window.values, p0=[data_window.idxmax().values, 30, 0.2])


    if fit_type == lorentzian:
        area_peak = np.abs(fit[0][1]*fit[0][2]*np.pi)
    elif fit_type == gaussian:
        area_peak = np.abs(np.sqrt(2*np.pi)*fit[0][1]*fit[0][2])

    if plot == True:
        ax = plt.gca()
        ax.axvline(fit_window[0], c='C1')
        ax.axvline(fit_window[1], c='C1')
        ax.plot(data.sel(wavenumber=slice(fit_window[1]+200, fit_window[0]-200))['wavenumber'], fit_type(data.sel(wavenumber=slice(fit_window[1]+200, fit_window[0]-200))['wavenumber'], *fit[0]), c=color,linestyle='--')
        return area_peak, peak_intensity
    elif plot == 'multi':
        fig1, ax = plt.subplots()
        line, = ax.plot(data.sel(wavenumber=slice(fit_window[1]+100, fit_window[0]-100))['wavenumber'], fit_type(data.sel(wavenumber=slice(fit_window[1]+100, fit_window[0]-100))['wavenumber'], *fit[0]), c=color,linestyle='--')
        #close the figure
        plt.close()
        return area_peak, peak_intensity, line
    elif plot == False:
        return area_peak, peak_intensity

def lorentzian(x, x0, gamma, A):
    return A * gamma**2 / ((x-x0)**2 + gamma**2)
def gaussian(x, x0, sigma, A):
    return A * np.exp(-(x-x0)**2 / (2*sigma**2))



def baseline_substract(bl_array,spectra_array): # for each spectrum, look up the specturm in the baseline array which has the closest temperature and substract it
    bl_temps = np.array(bl_array['temperature'])
    spectra_temps = np.array(spectra_array['temperature'])
    closest_indices = []
    for element in spectra_temps:
        closest_index = np.abs(bl_temps - element).argmin()
        closest_indices.append(closest_index)
    
    closest_temp_spectra = np.array(bl_array[closest_indices])
    spectra_array_corr = spectra_array - closest_temp_spectra
    return spectra_array_corr


def get_slice(data,high_wavenumber, low_wavenumber):
    return data.sel(wavenumber=slice(high_wavenumber, low_wavenumber))





# TPD analysis


#preparing TPD dataset
def cut_TPD_dataset(dataset):
    #selects only the spectra below 500C and substracts the minimum
    index_500C = np.where(dataset.temperature.values > 500.0)[0][0]
    dataset_500C = dataset[0:index_500C]
    return dataset_500C


#getting TPD profile
def get_tpd_BAS(dataset,pelettweight,peakloc=1545.0):
    dataset_p = cut_TPD_dataset(dataset)
    temps = dataset_p['temperature'].values
    integrals = []
    for i in range(len(dataset_p)):
        data_bl= linear_bl_corr(dataset_p[i])
        integral,intesity = fit_integrate_peak(data_bl, peakloc, 25, 15, 'lorentzian', plot=False)
        integrals.append(integral)
    df = pd.DataFrame({'temperature':temps,'integral':integrals,'integral_byweight':np.array(integrals)/pelettweight})
    return df
    


def linear_bl_corr(dataset, start = 1564, end =1508):
    return dataset - linear_baseline(dataset,start,end)


def linear_baseline(data, start, end):
    x = data.wavenumber
    x1 = start
    x2 = end
    y1 = data.sel(wavenumber=x1,method='nearest')
    y2 = data.sel(wavenumber=x2,method='nearest')
    m = (y2-y1)/(x2-x1)
    b = y1 - m*x1
    baseline = m*x+b
    return baseline