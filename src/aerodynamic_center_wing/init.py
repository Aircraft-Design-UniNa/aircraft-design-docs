# set the style of the notebook
from IPython.core.display import HTML
def css_styling():
    styles = open('./style/nbstyle.css', 'r').read()
    return HTML(styles)
css_styling()

# load libraries and set plot parameters
import math
import numpy as np
import tables as pt

from scipy.interpolate import interp1d, interp2d, Rbf

import h5py

from sympy import *

from IPython.display import display, Math, Latex, SVG

from cycler import cycler

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('pdf', 'png')
plt.rcParams['savefig.dpi'] = 96

plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.figsize'] = 10, 12
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 14

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"
plt.rcParams['font.serif'] = "cm"
# plt.rcParams['text.latex.preamble'] = "\usepackage{subdepth}, \usepackage{type1cm}"

plt.rc('axes', 
       prop_cycle=(
           cycler('color', ['black', 'darkblue', 'red', 'green', 'blue', 'brown', 'orange', 'darkgreen'])
           + cycler('linestyle',["-", "--", "-.", ":", ".", "h", "H","-."])
           + cycler('dashes',[[10000,1], [10,5], [10,5,2,5], [3,3], [20,5,5,5], [30,5], [10,10], [2,5,5,3]])
           )
    )

#----------------------------------------------------------
def plot_planform(c_r, c_t, b, Lambda_le, *args, **kwargs):

    # optional arguments
    c_mac = kwargs.get('mac', None)
    X_le_mac = kwargs.get('X_le_mac', None)
    Y_mac = kwargs.get('Y_mac', None)
    X_ac = kwargs.get('X_ac', None)
    
    xLineWing = [0,b/2,b/2,0]
    dy = (b/2)*math.tan(Lambda_le)
    yLineWing = [
        0, 0.25*c_r - 0.25*c_t + dy,
        0.25*c_r + 0.75*c_t + dy, c_r]

    # planform
    lineWing, = plt.plot(xLineWing, yLineWing, 'k-')
    # centerline
    centerLine, = plt.plot([0,0], [-1.1*c_r,2.1*c_r], 'b')
    centerLine.set_dashes([8, 4, 2, 4]) 
    # c/4 line
    pC4r = [0,0.25*c_r]
    pC4t = [b/2,dy + 0.25*c_r]
    quarterChordLine, = plt.plot([pC4r[0],pC4t[0]], [pC4r[1],pC4t[1]], 'k--')

    if ('mac' in kwargs) and ('X_le_mac' in kwargs) and ('Y_mac' in kwargs):
        c_mac = kwargs['mac']
        X_le_mac = kwargs['X_le_mac']
        Y_mac = kwargs['Y_mac']
        #print(mac)
        #print(X_le_mac)
        #print(Y_mac)
        lineMAC, = plt.plot([Y_mac, Y_mac], [X_le_mac, X_le_mac + c_mac], color="red", linewidth=2.5, linestyle="-")
        lineMAC.set_dashes([1000,1]) # HUUUUGE
        lineLEMAC, = plt.plot([0,b/2], [X_le_mac,X_le_mac], color="orange", linewidth=1.5, linestyle="-")
        lineLEMAC.set_dashes([10,2])
        lineTEMAC, = plt.plot([0,b/2], [X_le_mac + c_mac, X_le_mac + c_mac], color="orange", linewidth=1.5, linestyle="-")
        lineTEMAC.set_dashes([10,2])
        plt.scatter(Y_mac, X_le_mac, marker='o', s=40)
        ax = plt.gca()  # gca stands for 'get current axis'
        ax.annotate(
            r'$(Y_{\bar{c}},X_{\mathrm{le},\bar{c}}) = '
                +r'( {0:.3}'.format(Y_mac) + r'\,\mathrm{m}'+r',\,{0:.3}'.format(X_le_mac) + r'\,\mathrm{m} )$',
                         xy=(Y_mac, X_le_mac), xycoords='data',
                         xytext=(20, 30), textcoords='offset points', fontsize=12,
                         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2")) # 
        
    if 'X_ac' in kwargs:
        X_ac = kwargs['X_ac']
        #print(X_ac)
        plt.scatter(0, X_ac, marker='o', s=40)
        lineAC, = plt.plot([0,b/2], [X_ac,X_ac], color="brown", linewidth=3.5, linestyle="-")
        lineAC.set_dashes([10,2.5,3,2.5])
        ax = plt.gca()  # gca stands for 'get current axis'
        ax.annotate(r'$X_{\mathrm{ac,W}} = '+r'{0:.3}'.format(X_ac)+r'\,\mathrm{m} $',
                         xy=(b/2, X_ac), xycoords='data',
                         xytext=(20, 30), textcoords='offset points', fontsize=12,
                         arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2")) # 

    plt.axis('equal')
    
    #    xmajorLocator = MultipleLocator(2.0)
    #    xmajorFormatter = FormatStrFormatter('%.1f')
    #    xminorLocator = MultipleLocator(4)
    #    ax = plt.gca()  # gca stands for 'get current axis'
    #    ax.xaxis.set_major_locator(xmajorLocator)
    #    ax.xaxis.set_major_formatter(xmajorFormatter)
    #    # for the minor ticks, use no labels; default NullFormatter
    #    ax.xaxis.set_minor_locator(xminorLocator)
    
    plt.axis([-0.02*b/2, 1.1*b/2, -0.1*c_r, c_r + 1.1*dy])
    plt.gca().invert_yaxis()
    plt.title('Wing planform', fontsize=16)
    plt.xlabel('$y$ (m)', fontsize=16)
    plt.ylabel('$X$ (m)', fontsize=16)
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',c_r + 1.15*dy))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.1*b/2))

    plt.show()

#----------------------------------------------------------
def import_database_aerodynamic_center():
    fileName = "./resources/wing_aerodynamic_center.h5"
    f = h5py.File(fileName,'r',libver='latest')
    # K1
    data_K1 = f["(x_bar_ac_w)_k1_vs_lambda/data"]
    var0_K1 = f["(x_bar_ac_w)_k1_vs_lambda/var_0"]
    # K2
    data_K2 = f["(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/data"]
    var0_K2 = f["(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_0"]
    var1_K2 = f["(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_1"]
    var2_K2 = f["(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_2"]
    # xac/cr
    data_XacCr = f["(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/data"]
    var0_XacCr = f["(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_0"]
    var1_XacCr = f["(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_1"]
    var2_XacCr = f["(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_2"]
    return {
        'data_K1':data_K1, 
        'var0_K1':var0_K1, 
        'data_K2':data_K2, 
        'var0_K2':var0_K2, 
        'var1_K2':var1_K2, 
        'var2_K2':var2_K2,
        'data_XacCr':data_XacCr,
        'var0_XacCr':var0_XacCr,
        'var1_XacCr':var1_XacCr,
        'var2_XacCr':var2_XacCr
        }

#----------------------------------------------------------
def report_database_dimensions(database):
    shape_data_K1 = database['data_K1'].shape
    shape_var0_K1 = database['var0_K1'].shape
    
    shape_data_K2 = database['data_K2'].shape
    shape_var0_K2 = database['var0_K2'].shape
    shape_var1_K2 = database['var1_K2'].shape
    shape_var2_K2 = database['var2_K2'].shape
    
    shape_data_XacCr = database['data_XacCr'].shape
    shape_var0_XacCr = database['var0_XacCr'].shape
    shape_var1_XacCr = database['var1_XacCr'].shape
    shape_var2_XacCr = database['var2_XacCr'].shape
    
    print('(x_bar_ac_w)_k1_vs_lambda/var_0')
    print('shape of var0: {0}'.format(shape_var0_K1))
    print('(x_bar_ac_w)_k1_vs_lambda/data')
    print('shape of data: {0}'.format(shape_data_K1))
    """
    print('lambda --> K1')
    for i in range(shape_var0_K1[0]):
        print('{0}\t{1}'.format(dset_var0_K1[i],dset_data_K1[i]))
    """
    print('=====================================')
    print('(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_0')
    print('shape of var0: {0}'.format(shape_var0_K2))
    print('(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_1')
    print('shape of data: {0}'.format(shape_var1_K2))
    print('(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/var_2')
    print('shape of data: {0}'.format(shape_var2_K2))
    print('(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/data')
    print('shape of data: {0}'.format(shape_data_K2))
    print('=====================================')
    print("(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_0")
    print('shape of var0: {0}'.format(shape_var0_XacCr))
    print("(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_1")
    print('shape of data: {0}'.format(shape_var1_XacCr))
    print("(x_bar_ac_w)_x'_ac_over_root_chord_vs_tan_(L_LE)_over_beta_(AR_times_tan_(L_LE))_(lambda)/var_2")
    print('shape of data: {0}'.format(shape_var2_XacCr))
    print('(x_bar_ac_w)_k2_vs_L_LE_(AR)_(lambda)/data')
    print('shape of data: {0}'.format(shape_data_XacCr))

#----------------------------------------------------------
def plot_K1(var0_K1, data_K1):
    fig, ax = plt.subplots()
    plt.plot(var0_K1, data_K1, color="red", linewidth=2.5, linestyle="-")
    plt.title('Wing aerodynamic center --- effect of $\lambda$', fontsize=16)
    plt.xlabel('$\lambda$', fontsize=16)
    plt.ylabel('$K_1$', fontsize=16)
    plt.axis([0, 1, 0.8, 1.6])
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0.78))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.01))
    
    plt.show()

#----------------------------------------------------------
def plot_K2(var0_K2, var1_K2, var2_K2, data_K2, j_lambda=0):

    if j_lambda > 5 :
        print('Index j_lambda={0} out of range. Maximum allowed value 5'.format(j_lambda))
        return
    
    fig, ax = plt.subplots()
    
    #plt.gca().set_prop_cycle(
    #    cycler('color', ['c', 'm', 'y', 'k']) + cycler('lw', [1, 2, 3, 4]))
    
    idx_max_LambdaLE = 9
    
    for i_AR in range(0, 6):
        slice_ij = None
        slice_ij = data_K2[:,i_AR,j_lambda]
        line, = plt.plot(var2_K2, slice_ij, linewidth=2.5, linestyle="-")
        line.set_dashes([1000,1]) # HUUUUGE
        plt.annotate(r'$\mathrm{AR} = \,$'+r'{0}'.format(var1_K2[i_AR]),
                     xy=(var2_K2[idx_max_LambdaLE], slice_ij[idx_max_LambdaLE]), xycoords='data',
                     xytext=(40, 0), textcoords='offset points', fontsize=22,
                     arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"
    plt.title(
        'Wing aerodynamic center --- effect of $(\Lambda_{\mathrm{le}},\mathrm{AR})$, '
        +r'$\lambda = {0:.3}$'.format(var0_K2[j_lambda]),
        fontsize=22)
    
    plt.axis([0, 45, 0, 1.1*max(data_K2[:,5,j_lambda])])
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',-0.05))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.5))
    
    plt.xlabel('$\Lambda_{\mathrm{le}}$ (deg)', fontsize=22)
    plt.ylabel('$K_2$', fontsize=22)
    plt.show()

#----------------------------------------------------------
def multiplot_K2(var0_K2, var1_K2, var2_K2, data_K2):

    # j_lambda max = 5
    
    idx_max_LambdaLE = 9

    # fig, ax = plt.subplots()
    fig, axes = plt.subplots(3,2,figsize=(9, 11))
    
    j_lambda = 0
    for i, ax in enumerate(axes.flat, start=1):
        for i_AR in range(0, 6):
            slice_ij = None
            slice_ij = data_K2[:,i_AR,j_lambda]
            line, = ax.plot(var2_K2, slice_ij, linewidth=2.5, linestyle="-")
            line.set_dashes([1000,1]) # HUUUUGE
            ax.annotate(r'$\mathrm{AR} = \,$'+r'{0}'.format(var1_K2[i_AR]),
                         xy=(var2_K2[idx_max_LambdaLE], slice_ij[idx_max_LambdaLE]), xycoords='data',
                         xytext=(20, 0), textcoords='offset points', fontsize=12,
                         arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"
        ax.set_title(
            #'Wing aerodynamic center --- effect of $(\Lambda_{\mathrm{le}},\mathrm{AR})$, '
            r'$\lambda = {0:.3}$'.format(var0_K2[j_lambda]),
            fontsize=12)
        
        ymajorLocator = MultipleLocator(0.5)
        ymajorFormatter = FormatStrFormatter('%.1f')
        yminorLocator = MultipleLocator(5)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_major_formatter(ymajorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax.yaxis.set_minor_locator(yminorLocator)

        ax.axis([0, 45, 0, 1.1*max(data_K2[:,5,j_lambda])])
        
        # Moving spines
        #ax = plt.gca()  # gca stands for 'get current axis'
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['bottom'].set_position(('data',-0.05))
        ax.yaxis.set_ticks_position('left')
        ax.spines['left'].set_position(('data',-0.5))
        
        ax.set_xlabel('$\Lambda_{\mathrm{le}}$ (deg)', fontsize=12)
        ax.set_ylabel('$K_2$', fontsize=12)
        #ax[0,0].show()
        j_lambda += 1 

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.subplots_adjust(wspace=0.78)
    plt.show()
    

#----------------------------------------------------------
def plot_XacCr(var0_XacCr, var1_XacCr, var2_XacCr, data_XacCr, j_lambda=0):
    
    if j_lambda > 5 :
        print('Index j_lambda={0} out of range. Maximum allowed value 5'.format(j_lambda))
        return

    fig, ax = plt.subplots()
    
    #plt.gca().set_prop_cycle(
    #    cycler('color', ['c', 'm', 'y', 'k']) + cycler('lw', [1, 2, 3, 4]))
    
    idx_max_TanLambdaLE = 10
    
    for i_AR in range(0, 7):
        slice_ij = None
        slice_ij = data_XacCr[:,i_AR,j_lambda]
        line, = plt.plot(var2_XacCr, slice_ij, linewidth=2.5, linestyle="-")
        line.set_dashes([1000,1]) # HUUUUGE
        plt.annotate(r'$\mathrm{AR} \tan\Lambda_{\mathrm{le}} =\,$'+r'{0}'.format(var1_XacCr[i_AR,0]),
                     xy=(var2_XacCr[idx_max_TanLambdaLE], slice_ij[idx_max_TanLambdaLE]), xycoords='data',
                     xytext=(40, 0), textcoords='offset points', fontsize=22,
                     arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"
    plt.title(
        r'Wing aerodynamic center --- effect of $(\tan\Lambda_{\mathrm{le}}/\sqrt{1-M^2},\mathrm{AR}\tan\Lambda_{\mathrm{le}})$, '
        +'$\lambda = {0:.3}$'.format(var0_XacCr[j_lambda]),
        fontsize=22)
    
    plt.axis([0, 2.2, -0.05, 1.1*max(data_XacCr[:,6,j_lambda])])
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',-0.07))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.05))
    
    plt.xlabel(r'$\tan\Lambda_{\mathrm{le}}/\sqrt{1-M^2}$', fontsize=22)
    plt.ylabel('$X_{\mathrm{ac}}\'/c_{\mathrm{r}}$', fontsize=22)
    plt.show()

#----------------------------------------------------------
def multiplot_XacCr(var0_XacCr, var1_XacCr, var2_XacCr, data_XacCr):

    # j_lambda max = 5
    
    idx_max_TanLambdaLE = 10

    # fig, ax = plt.subplots()
    fig, axes = plt.subplots(3,2,figsize=(9, 11))
    
    j_lambda = 0
    for i, ax in enumerate(axes.flat, start=1):

        for i_AR in range(0, 7):
            slice_ij = None
            slice_ij = data_XacCr[:,i_AR,j_lambda]
            line, = ax.plot(var2_XacCr, slice_ij, linewidth=2.5, linestyle="-")
            line.set_dashes([1000,1]) # HUUUUGE
            ax.annotate(r'$\mathrm{AR}\tan\Lambda_{\mathrm{le}} = \,$'+r'{0}'.format(var1_XacCr[i_AR,0]),
                         xy=(var2_XacCr[idx_max_TanLambdaLE], slice_ij[idx_max_TanLambdaLE]), xycoords='data',
                         xytext=(20, 0), textcoords='offset points', fontsize=12,
                         arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"
        ax.set_title(
            #'Wing aerodynamic center --- effect of $(\Lambda_{\mathrm{le}},\mathrm{AR})$, '
            r'$\lambda = {0:.3}$'.format(var0_XacCr[j_lambda]),
            fontsize=12)

        ymajorLocator = MultipleLocator(0.25)
        ymajorFormatter = FormatStrFormatter('%.1f')
        yminorLocator = MultipleLocator(1)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_major_formatter(ymajorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax.yaxis.set_minor_locator(yminorLocator)

        ax.axis([0, 2, -0.05, 1.1*max(data_XacCr[:,6,j_lambda])])
        
        # Moving spines
        #ax = plt.gca()  # gca stands for 'get current axis'
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['bottom'].set_position(('data',-0.07))
        ax.yaxis.set_ticks_position('left')
        ax.spines['left'].set_position(('data',-0.05))
        
        ax.set_xlabel(r'$\tan\Lambda_{\mathrm{le}}/\sqrt{1-M^2}$', fontsize=12)
        ax.set_ylabel('$X_{\mathrm{ac}}\'/c_{\mathrm{r}}$', fontsize=12)

        #ax[0,0].show()
        j_lambda += 1 

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.subplots_adjust(wspace=0.85)
    plt.show()
    


#----------------------------------------------------------
def plot_interpolate_K1(var0_K1, data_K1, lam):
    g = interp1d(var0_K1, data_K1)
    K1 = g(lam)
    print('lambda = {0} --> K_1 = {1}'.format(lam,K1))
    fig, ax = plt.subplots()
    plt.plot(var0_K1, data_K1, color="brown", linewidth=2.5, linestyle="-")

    # interpolated data
    plt.scatter(lam, K1, marker='o', s=40)
    help_line, = plt.plot([lam,lam,0],[0,K1,K1], color="red", linewidth=1.5, linestyle="--")
    help_line.set_dashes([8, 4]) 
    
    plt.title('Wing aerodynamic center --- effect of $\lambda$', fontsize=22)
    plt.xlabel('$\lambda$', fontsize=16)
    plt.ylabel('$K_1$', fontsize=16)
    plt.axis([0, 1, 0.8, 1.6])
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0.78))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.01))
    
    plt.show()

#----------------------------------------------------------
def plot_interpolate_K2(var0_K2, var1_K2, var2_K2, data_K2, j_lambda, 
                        LamLE_deg, AR):

    if j_lambda > 5 :
        print('Index j_lambda={0} out of range. Maximum allowed value 5'.format(j_lambda))
        return

    # interpolation in 2D
    g = interp2d(var1_K2, var2_K2, data_K2[:,:,j_lambda], kind='linear')
    K2 = g(AR, LamLE_deg)
    print('Lambda_LE = {0} deg, AR = {1} --> K_2 = {2}'.format(LamLE_deg, AR, K2[0]))
    
    fig, ax = plt.subplots()
    
    #plt.gca().set_prop_cycle(
    #    cycler('color', ['c', 'm', 'y', 'k']) + cycler('lw', [1, 2, 3, 4]))
    
    idx_max_LambdaLE = 9
    
    for i_AR in range(0, 6):
        slice_ij = None
        slice_ij = data_K2[:,i_AR,j_lambda]
        line, = plt.plot(var2_K2, slice_ij, linewidth=2.5, linestyle="-")
        line.set_dashes([1000,1]) # HUUUUGE
        plt.annotate(r'$\mathrm{AR} = \,$'+r'{0}'.format(var1_K2[i_AR]),
                     xy=(var2_K2[idx_max_LambdaLE], slice_ij[idx_max_LambdaLE]), xycoords='data',
                     xytext=(40, 0), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"

    # interpolated data
    plt.scatter(LamLE_deg, K2, marker='o', s=40)
    help_line, = plt.plot([LamLE_deg,LamLE_deg,0],[0,K2,K2], color="red", linewidth=1.5, linestyle="--")
    help_line.set_dashes([8, 4]) 
    
    plt.title(
        'Wing aerodynamic center --- effect of $(\Lambda_{\mathrm{le}},\mathrm{AR})$, '
        +r'$\lambda = {0:.3}$'.format(var0_K2[j_lambda]),
        fontsize=16)
    
    plt.axis([0, 45, 0, 1.1*max(data_K2[:,5,j_lambda])])
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',-0.03))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.5))
    
    plt.xlabel('$\Lambda_{\mathrm{le}}$ (deg)', fontsize=16)
    plt.ylabel('$K_2$', fontsize=16)
    plt.show()

#----------------------------------------------------------
def plot_interpolate_XacCr(var0_XacCr, var1_XacCr, var2_XacCr, data_XacCr, j_lambda, 
                           LamLE_deg, AR, Mach):
    
    if j_lambda > 5 :
        print('Index j_lambda={0} out of range. Maximum allowed value 5'.format(j_lambda))
        return

    # interpolation in 2D
    g = interp2d(var1_XacCr, var2_XacCr, data_XacCr[:,:,j_lambda], kind='linear')
    y_ = AR*math.tan(LamLE_deg*np.pi/180)
    x_ = math.tan(LamLE_deg*np.pi/180)/math.sqrt(1 - math.pow(Mach,2))
    XacCr = g(y_,x_)
    print('x_1 = tan(Lambda_LE)/sqrt(1-M^2) = {0}\nx_2 = AR*tan(Lambda_LE) = {1}\n --> Xac\'/c_r = {2}'
        .format(x_, y_, XacCr[0]))

    fig, ax = plt.subplots()
    
    #plt.gca().set_prop_cycle(
    #    cycler('color', ['c', 'm', 'y', 'k']) + cycler('lw', [1, 2, 3, 4]))
    
    idx_max_TanLambdaLE = 10
    
    for i_AR in range(0, 7):
        slice_ij = None
        slice_ij = data_XacCr[:,i_AR,j_lambda]
        line, = plt.plot(var2_XacCr, slice_ij, linewidth=2.5, linestyle="-")
        line.set_dashes([1000,1]) # HUUUUGE
        plt.annotate(r'$\mathrm{AR} \tan\Lambda_{\mathrm{le}} =\,$'+r'{0}'.format(var1_XacCr[i_AR,0]),
                     xy=(var2_XacCr[idx_max_TanLambdaLE], slice_ij[idx_max_TanLambdaLE]), xycoords='data',
                     xytext=(40, 0), textcoords='offset points', fontsize=16,
                     arrowprops=dict(arrowstyle="->")) # , connectionstyle="arc3,rad=.5"

    # interpolated data
    plt.scatter(x_, XacCr, marker='o', s=40)
    help_line, = plt.plot([x_,x_,0],[0,XacCr,XacCr], color="red", linewidth=1.5, linestyle="-")
    help_line.set_dashes([8, 4]) 
    
    plt.title(
        r'Wing aerodynamic center --- effect of $(\tan\Lambda_{\mathrm{le}}/\sqrt{1-M^2},\mathrm{AR}\tan\Lambda_{\mathrm{le}})$, '
        +'$\lambda = {0:.3}$'.format(var0_XacCr[j_lambda]),
        fontsize=16)
    
    plt.axis([0, 2.2, -0.05, 1.1*max(data_XacCr[:,6,j_lambda])])
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',-0.07))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.05))
    
    plt.xlabel(r'$\tan\Lambda_{\mathrm{le}}/\sqrt{1-M^2}$', fontsize=16)
    plt.ylabel('$X_{\mathrm{ac}}\'/c_{\mathrm{r}}$', fontsize=16)
    plt.show()

#----------------------------------------------------------
def display_workflow_c_mac(S_ref, b, A_c, B_c, c_mac_law_integral_indefinite, mac):
    c_mac_law_integral_indefinite_latex = latex(c_mac_law_integral_indefinite)
    return Latex(
        r'\begin{align*}'
            + r'\bar{c} & {}=\frac{2}{S} \int_{0}^{b/2} c^2(y) \, \mathrm{d}y'
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'\,\mathrm{m}^2} \int_0^{' + '{0:.3}'.format(b/2) + '}' 
            + r'\Big(' 
            +   r'{0:.4}'.format(A_c) + r'\, y + ' + '{0:.4}'.format(B_c) + r'\,\text{m}'
            + r'\Big)^2\,\text{d}y' 
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'} \big(' 
            +   c_mac_law_integral_indefinite_latex 
            + r'\big)\Bigr|_0^{' + '{0:.4}'.format(b/2) + r'}'
            + '=' + '{0:.4}'.format(mac) + r'\,\text{m}'
        + r'\end{align*}'
    )

#----------------------------------------------------------
def display_workflow_X_le_mac(S_ref, b, A_c, B_c, A_xle, X_le_mac_law_integral_indefinite, X_le_mac):
    X_le_mac_law_integral_indefinite_latex = latex(X_le_mac_law_integral_indefinite)
    return Latex(
        r'\begin{align*}'
            + r'X_{\mathrm{le},\bar{c}} & {}=\frac{2}{S} \int_{0}^{b/2} X_{\mathrm{le}}(y) \, c(y) \, \mathrm{d}y'
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'\,\mathrm{m}^2} \int_0^{' + '{0:.3}'.format(b/2) + '}' 
            + r'{0:.4}'.format(A_xle) + r'\, y\,'            
            + r'\Big(' 
            +   r'{0:.4}'.format(A_c) + r'\, y + ' + '{0:.4}'.format(B_c) + r'\,\text{m}'
            + r'\Big)\,\text{d}y' 
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'} \big(' 
            +   X_le_mac_law_integral_indefinite_latex 
            + r'\big)\Bigr|_0^{' + '{0:.4}'.format(b/2) + r'}'
            + '=' + '{0:.4}'.format(X_le_mac) + r'\,\text{m}'
        + r'\end{align*}'
    )

#----------------------------------------------------------
def display_workflow_Y_mac(S_ref, b, A_c, B_c, Y_mac_law_integral_indefinite, Y_mac):
    Y_mac_law_integral_indefinite_latex = latex(Y_mac_law_integral_indefinite)
    return Latex(
        r'\begin{align*}'
            + r'Y_{\bar{c}} & {}=\frac{2}{S} \int_{0}^{b/2} y \, c(y) \, \mathrm{d}y'
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'\,\mathrm{m}^2} \int_0^{' + '{0:.3}'.format(b/2) + '}' 
            + r' y\,'            
            + r'\Big(' 
            +   r'{0:.4}'.format(A_c) + r'\, y + ' + '{0:.4}'.format(B_c) + r'\,\text{m}'
            + r'\Big)\,\text{d}y' 
            + r'\\[1em]' 
            + r'& {}= \frac{2}{'+ '{0:.4}'.format(S_ref) +r'} \big(' 
            +   Y_mac_law_integral_indefinite_latex 
            + r'\big)\Bigr|_0^{' + '{0:.4}'.format(b/2) + r'}'
            + '=' + '{0:.4}'.format(Y_mac) + r'\,\text{m}'
        + r'\end{align*}'
    )

#----------------------------------------------------------
def display_workflow_K2(lam_a, lam_b, K2_a, K2_b, taper_ratio, K2):
    return Latex(
        r'\begin{equation}'
        +  r'K_2 = K_2 \Big|_{\lambda=' + '{0:.3}'.format(lam_a) + r'}'
        +    r'+ \frac{'
        +      r'K_2 \Big|_{\lambda=' + '{0:.3}'.format(lam_b) + r'}'
        +      r'- K_2 \Big|_{\lambda=' + '{0:.3}'.format(lam_a) + r'}'
        +    r'}{'
        +      r'{0:.3}'.format(lam_b) + r'-{0:.3}'.format(lam_a)
        +    r'}'
        +    r'\,\big( {0:.3} - {1:.3} \big)'.format(taper_ratio,lam_a)
        +    r'= {0:.3}'.format(K2_a)
        +    r'+ \frac{'
        +      r'{0:.3}'.format(K2_b)
        +      r'- {0:.3}'.format(K2_a)
        +    r'}{'
        +      r'{0:.3}'.format(lam_b) + r'-{0:.3}'.format(lam_a)
        +    r'}'
        +    r'\,\big( {0:.3} - {1:.3} \big)'.format(taper_ratio,lam_a)
        +    r'=' + '{0:.3}'.format(K2)
        +r'\end{equation}'
    )

#----------------------------------------------------------
def display_workflow_XacCr(lam_a, lam_b, XacCr_a, XacCr_b, taper_ratio, XacCr):
    return Latex(
        r'\begin{equation}'
        +  r'\frac{X_{\mathrm{ac}}'+"'"+r'}{c_{\mathrm{r}}} = \frac{X_{\mathrm{ac}}'+"'"+r'}{c_{\mathrm{r}}} \Big|_{\lambda=' + '{0:.3}'.format(lam_a) + r'}'
        +    r'+ \frac{'
        +      r'\dfrac{X_{\mathrm{ac}}'+"'"+r'}{c_{\mathrm{r}}} \Big|_{\lambda=' + '{0:.3}'.format(lam_b) + r'}'
        +      r'- \dfrac{X_{\mathrm{ac}}'+"'"+r'}{c_{\mathrm{r}}} \Big|_{\lambda=' + '{0:.3}'.format(lam_a) + r'}'
        +    r'}{'
        +      r'{0:.3}'.format(lam_b) + r'-{0:.3}'.format(lam_a)
        +    r'}'
        +    r'\,\big( {0:.3} - {1:.3} \big)'.format(taper_ratio,lam_a)
        +    r'= {0:.3}'.format(XacCr_a)
        +    r'+ \frac{'
        +      r'{0:.3}'.format(XacCr_b)
        +      r'- {0:.3}'.format(XacCr_a)
        +    r'}{'
        +      r'{0:.3}'.format(lam_b) + r'-{0:.3}'.format(lam_a)
        +    r'}'
        +    r'\,\big( {0:.3} - {1:.3} \big)'.format(taper_ratio,lam_a)
        +    r'=' + '{0:.3}'.format(XacCr)
        +r'\end{equation}'
    )