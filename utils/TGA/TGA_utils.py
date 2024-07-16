import pandas as pd
import numpy as np
import io
import re
from matplotlib import pyplot as plt
import os
import warnings

# 
class TGA_exp:
    #General class for TGA experiment with multiple stages
    def __init__(self, stage_files=None):
        self.stages = {}
        if stage_files is not None:
            for stage, file in stage_files.items():
                data = pd.read_csv(file)
                self.add_stage(stage, data)

    def add_stage(self, stage, data):
        if not isinstance(data, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame")
        self.stages[stage] = data

    def get_stage(self, stage):
        return self.stages.get(stage, None)
    def stage_names(self):
        return list(self.stages.keys())




class TGA_pyro(TGA_exp):
    #TGA class for a ramped plastic cracking experiments. Key stages are cracking(heating under N2) and burnoff(heating under O2).
    def __init__(self, stage_files=None):
        super().__init__(stage_files)
        self.Tmax = None
        self.T50 = None
    def cracking(self):
        return self.stages['stage4'] # If the TAa method is changed, this needs to be adjusted accordingly
    def burnoff(self):
        return self.stages['stage8']
    def m_cat(self):# returns the amount of catalyst
        return self.burnoff()['Unsubtracted weight'].min()
    def m_poly(self): # returns the amount of polymer
        return self.cracking()['Unsubtracted weight'].max() - self.m_cat()
    def m_coke(self):
        return self.cracking()['Unsubtracted weight'].min()-self.m_cat()
    def pct_loss(self):
        return self.m_poly()/(self.m_poly()+self.m_cat())
    def P_C_ratio(self):
        return self.m_poly()/self.m_cat()
    def coke_yield(self):
        return self.m_coke()/self.m_poly()


class TGA_pyro_iso(TGA_exp):
    #TGA for a isothermal plastic cracking experiments. 
    def __init__(self, stage_files=None):
        super().__init__(stage_files)
        self.Tmax = None
        self.T50 = None
    def cracking(self):
        return self.stages['stage5']
    def burnoff(self):
        return self.stages['stage7']
    def m_cat(self):# returns the amount of catalyst
        return self.burnoff()['Unsubtracted weight'].min()
    def m_poly(self): # returns the amount of polymer
        return self.cracking()['Unsubtracted weight'].max() - self.m_cat()
    def m_coke(self):
        return self.cracking()['Unsubtracted weight'].min()-self.m_cat()
    def pct_loss(self):
        return self.m_poly()/(self.m_poly()+self.m_cat())
    def P_C_ratio(self):
        return self.m_poly()/self.m_cat()
    def temp(self):
        return np.round(self.cracking()['Sample Temp.'].iloc[-1],0)



###

def parse_txt(filepath,type = 'general',calculate_DTGA = True): # type can be 'general' or 'pyro'
    #Parser for PerkinElmer TGA8000 ASCII output files. Splits the file into stages using an regex.
    if type == 'general':
        tga_exp_instance = TGA_exp()  # Create an instance of TGA_exp
    elif type == 'pyro':
        tga_exp_instance = TGA_pyro()
    elif type == 'pyro_iso':
        tga_exp_instance = TGA_pyro_iso()
    else:
        raise ValueError("type must be 'general' or 'pyro'")
    names = ['Blank', 'Time', 'Unsubtracted weight', 'Baseline weight',
            'Program Temp.', 'Sample Temp.', 'Sample Purge Flow',
            'Balance purge flow']

    def read_section(data):
        frame = pd.read_table(io.StringIO(data), sep='\t', header=None,skiprows=2,
                                names=names, engine='python')
        return frame

    with open(filepath) as full:
        text = full.read()
        split = re.split(r'(\d+\) TGA)', text) #important bit right here
        section_numbers = [int(element[0:-5]) for element in split[1::2]]
        sections = split[2::2]

    for i in range(len(section_numbers)): 
        tga_exp_instance.add_stage('stage'+str(section_numbers[i]), read_section(sections[i]))

    if calculate_DTGA == True:
        return calc_DTGA(tga_exp_instance)
    else:
        return tga_exp_instance

def calc_DTGA(tga_exp, averaging_window = 30): #smooths the derivative with a moving average over a window of 30 datapoints
    # calculates the derivative of the TGA curve for a plastic cracking experimen in units of total sample weight or just the plastic weight.
    # adds this derivative as a column to the dataframe of each stage.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for stage in [tga_exp.cracking(),tga_exp.burnoff()]:
            stage.drop(stage.index[stage['Sample Temp.'] == stage['Sample Temp.'].iloc[0]], inplace=True) # removing a couple datapoints to avoid infinities

            rel_weight_twl = (stage['Unsubtracted weight']/stage['Unsubtracted weight'].max()).to_numpy()
            rel_weight_pwl = ((stage['Unsubtracted weight']-tga_exp.m_cat())/tga_exp.m_poly()).to_numpy()

            stage['rel_weight_twl'] = rel_weight_twl
            stage['rel_weight_pwl'] = rel_weight_pwl

            temp = stage['Sample Temp.'].to_numpy()

            stage['DTGA_pwl'] = -np.gradient(rel_weight_pwl,temp)
            stage['DTGA_pwl']=stage['DTGA_pwl'].rolling(averaging_window,win_type='triang').mean()
            stage['DTGA_twl'] = -np.gradient(rel_weight_twl,temp)
            stage['DTGA_twl']=stage['DTGA_twl'].rolling(averaging_window,win_type='triang').mean()
    return tga_exp

def calc_Tmax(tga_exp,stage='cracking'):
    if stage == 'cracking':
        stage_select = tga_exp.cracking()
    elif stage == 'burnoff':
        stage_select = tga_exp.burnoff()
    Tmax = stage_select['Sample Temp.'].iloc[stage_select['DTGA_twl'].idxmax()]
    return Tmax

def calc_T50(tga_exp,stage='cracking'):
    if stage == 'cracking':
        stage_select = tga_exp.cracking()
    elif stage == 'burnoff':
        stage_select = tga_exp.burnoff()
    T50 = stage_select['Sample Temp.'].iloc[stage_select['rel_weight_pwl'].sub(0.5).abs().idxmin()]
    return T50




def get_color(min_rel_weight,cmap='viridis'):
    norm = plt.Normalize(0, 1.07)
    color = plt.get_cmap(cmap)(norm(min_rel_weight))
    return color