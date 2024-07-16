import matplotlib.pyplot as plt
import os
if os.name == 'posix':
    def set_plot_params():
        font = {'family' : 'Liberation Sans',
                'weight' : 'normal',
                'size'   : 7}
        plt.rc('font', **font)
        plt.rcParams['figure.figsize'] = [8.8/2.54,6.22/2.54]  # Set figure size


    def set_plot_params_2():
        font = {'family' : 'Liberation Sans',
                'weight' : 'normal',
                'size'   : 7}
        plt.rc('font', **font)
        plt.rcParams['figure.figsize'] = [18/2.54,12.73/2.54]
else:
    def set_plot_params():
        font = {'family' : 'Arial',
                'weight' : 'normal',
                'size'   : 7}
        plt.rc('font', **font)
        plt.rcParams['figure.figsize'] = [8.8/2.54,6.22/2.54]  # Set figure size


    def set_plot_params_2():
        font = {'family' : 'Arial',
                'weight' : 'normal',
                'size'   : 7}
        plt.rc('font', **font)
        plt.rcParams['figure.figsize'] = [18/2.54,12.73/2.54]


set_plot_params()