# set the style of the notebook
from IPython.core.display import HTML
def css_styling():
    styles = open('./style/nbstyle.css', 'r').read()
    return HTML(styles)
css_styling()

# load libraries and set plot parameters
import math
import numpy as np
import pandas as pd
import tables as pt

import scipy
from scipy.interpolate import interp1d, interp2d

import h5py

# https://docs.python.org/3.5/library/shelve.html
import shelve

import sympy

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
def plot_planform(c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2, *args, **kwargs):

    fig = plt.subplots(figsize=(9, 9))
    
    # optional arguments
    c_mac = kwargs.get('mac', None)
    X_le_mac = kwargs.get('X_le_mac', None)
    Y_mac = kwargs.get('Y_mac', None)
    X_ac = kwargs.get('X_ac', None)
    
    xLineWing = [0, b_k/2, b/2, b/2, b_k/2, 0]
    dy_k = (b_k/2)*math.tan(Lambda_le_1)
    dy = dy_k + (b/2 - b_k/2)*math.tan(Lambda_le_2)
    yLineWing = [
        0, 
        dy_k, 
        dy,
        dy + c_t, 
        dy_k + c_k, 
        c_r]
        
    # planform
    lineWing, = plt.plot(xLineWing, yLineWing, 'k-')

    plt.scatter(xLineWing, yLineWing, marker='o', s=40)    
    
    # centerline
    centerLine, = plt.plot([0,0], [-0.2*c_r,2.1*c_r], 'b')
    centerLine.set_dashes([8, 4, 2, 4]) 
    # c/4 line
    pC4r = [0, 0.25*c_r]
    pC4k = [b_k/2, dy_k + 0.25*c_k]
    pC4t = [b/2, dy + 0.25*c_t]
    quarterChordLine, = plt.plot([pC4r[0],pC4k[0],pC4t[0]], [pC4r[1],pC4k[1],pC4t[1]], 'k--')
    plt.scatter([pC4r[0],pC4k[0],pC4t[0]], [pC4r[1],pC4k[1],pC4t[1]], marker='o', s=40)

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

    if ('X_le_r_eq' in kwargs) and ('c_r_eq' in kwargs):
        X_le_r_eq = kwargs['X_le_r_eq']
        c_r_eq = kwargs['c_r_eq']
        vertices = [(0, X_le_r_eq)] + [(b/2, dy)] + [(b/2, dy+c_t)] + [(0, X_le_r_eq + c_r_eq)]
        poly = Polygon(vertices, facecolor="yellow", alpha=0.5)
        poly.set_edgecolor("brown")
        poly.set_linewidth(2)
        ax0 = plt.gca()  # gca stands for 'get current axis'
        ax0.add_patch(poly)
        
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
    
    plt.axis([-0.02*b/2, 1.1*b/2, -0.05*c_r, 1.1*(dy + c_t)])
    plt.gca().invert_yaxis()
    plt.title('Wing planform', fontsize=16)
    plt.xlabel('$y$ (m)', fontsize=16)
    plt.ylabel('$X$ (m)', fontsize=16)
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('outward',10)) # outward by 10 points
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.07*b/2))

    plt.show()

#----------------------------------------------------------
def plot_wing_functions(c_r, c_k, c_t, 
                    eps_k, eps_t, alpha0l_r, alpha0l_k, alpha0l_t,
                    b_k, b, Lambda_le_1, Lambda_le_2, 
                    *args, **kwargs):
    """
    
    See: http://www.scipy-lectures.org/intro/matplotlib/matplotlib.html
    
    """

    # optional arguments
    f_chord = kwargs.get('f_chord', None)
    f_Xle = kwargs.get('f_Xle', None)
    f_twist = kwargs.get('f_twist', None)
    f_alpha0l = kwargs.get('f_alpha0l', None)
    f_S_integral = kwargs.get('f_S_integral', None)
    f_mac_integral = kwargs.get('f_mac_integral', None)
    f_Xle_mac_integral = kwargs.get('f_Xle_mac_integral', None)
    f_Y_mac_integral = kwargs.get('f_Y_mac_integral', None)
    f_alpha0L_integral_indefinite = kwargs.get('f_alpha0L_integral_indefinite', None)

    n_points = kwargs.get('f_chord', None)
    if ('n_points' in kwargs):
        n_points = kwargs['n_points']
    else:
        n_points = 20

    # define vectors
    vY0 = np.linspace(0, b/2, n_points, endpoint=True)
    vY1 = np.concatenate([vY0,[b_k/2]])
    vY = np.sort(np.unique(vY1))

    the_figsize = kwargs.get('figsize', None)
    if ('figsize' in kwargs):
        the_figsize = kwargs['figsize']
    else:
        the_figsize = (11,12)


    # Create a figure of size WxH inches, DPI dots per inch
    fig = plt.figure(figsize=the_figsize, dpi=300)
    
    # Create a new subplot from a grid of 1x1
    ax0 = plt.subplot(1, 1, 1)
    
    #fig, ax0 = plt.subplots()
    
    if ('f_chord' in kwargs):
        y = sympy.Symbol('y')
        f_chord = kwargs['f_chord']
        vChord = []
        for y in vY:
            vChord.append(f_chord(y,c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2))
        vChord = np.asarray(vChord)
        plt.plot(vY, vChord, color="red", linewidth=2.5, linestyle="-", label=r'local chord $c$ (m)')

    if ('f_Xle' in kwargs):
        y = sympy.Symbol('y')
        f_Xle = kwargs['f_Xle']
        vXle = []
        for y in vY:
            vXle.append(f_Xle(y, b_k, b, Lambda_le_1, Lambda_le_2))
        vXle = np.asarray(vXle)
        plt.plot(vY, vXle, color="green", linewidth=2.5, linestyle="-", label=r'local l.e. coordinate $X_{\mathrm{le}}$ (m)')
        
    if ('f_twist' in kwargs):
        y = sympy.Symbol('y')
        f_twist = kwargs['f_twist']
        vTwist = []
        for y in vY:
            vTwist.append(f_twist(y, eps_k, eps_t, b_k, b))
        vTwist = np.asarray(vTwist)
        plt.plot(vY, vTwist*180/np.pi, color="blue",  linewidth=2.5, linestyle="-", label=r"local $\epsilon_{\mathrm{g}}$ (deg)")

    if ('f_alpha0l' in kwargs):
        y = sympy.Symbol('y')
        f_alpha0l = kwargs['f_alpha0l']
        vAlpha0l = []
        for y in vY:
            vAlpha0l.append(f_alpha0l(y, alpha0l_r, alpha0l_k, alpha0l_t, b_k, b))
        vAlpha0l = np.asarray(vAlpha0l)
        plt.plot(vY, vAlpha0l*180/np.pi, color="brown",  linewidth=2.5, linestyle="-", label=r"local $\alpha_{0\ell}$ (deg)")

    if ('f_S_integral' in kwargs):
        f_S_integral = kwargs['f_S_integral']
        vS_integrand = []
        for y_ in vY:
            #print(y_)
            vS_integrand.append(f_S_integral(y_, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2))
        vS_integrand = np.asarray(vS_integrand)
        plt.plot(vY, vS_integrand,  linewidth=2.5, linestyle="-", dashes=[1000,1], marker="." ,
                 label=r"function $c(y)$ (m)")         
        # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
        vertices = [(0, 0)] + list(zip(vY, vS_integrand)) + [(b/2, 0)]
        poly = Polygon(vertices, facecolor="green", alpha=0.3, edgecolor="none")
        ax0.add_patch(poly)

    if ('f_mac_integral' in kwargs):
        f_mac_integral = kwargs['f_mac_integral']
        vMac_integrand = []
        for y_ in vY:
            vMac_integrand.append(f_mac_integral(y_, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2))
        vMac_integrand = np.asarray(vMac_integrand)
        plt.plot(vY, vMac_integrand, linewidth=2.5, linestyle="-", dashes=[1000,1], marker="." ,
                 label=r"function $c^2(y)$ (m${}^2$)")         
        # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
        vertices = [(0, 0)] + list(zip(vY, vMac_integrand)) + [(b/2, 0)]
        poly = Polygon(vertices, facecolor="orange", alpha=0.3, edgecolor="none")
        ax0.add_patch(poly)

    if ('f_Xle_mac_integral' in kwargs):
        f_Xle_mac_integral = kwargs['f_Xle_mac_integral']
        vXle_mac_integrand = []
        for y_ in vY:
            vXle_mac_integrand.append(f_Xle_mac_integral(y_, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2))
        vXle_mac_integrand = np.asarray(vXle_mac_integrand)
        plt.plot(vY, vXle_mac_integrand, linewidth=2.5, linestyle="-", dashes=[1000,1], marker="." ,
                 label=r"function $X_{\mathrm{le}}(y)\,c(y)$ (m${}^2$)")         
        # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
        vertices = [(0, 0)] + list(zip(vY, vXle_mac_integrand)) + [(b/2, 0)]
        poly = Polygon(vertices, facecolor="blue", alpha=0.3, edgecolor="none")
        ax0.add_patch(poly)

    if ('f_Y_mac_integral' in kwargs):
        f_Y_mac_integral = kwargs['f_Y_mac_integral']
        vY_mac_integrand = []
        for y_ in vY:
            vY_mac_integrand.append(f_Y_mac_integral(y_, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2))
        vY_mac_integrand = np.asarray(vY_mac_integrand)
        plt.plot(vY, vY_mac_integrand, linewidth=2.5, linestyle="-", dashes=[1000,1], marker="." ,
                 label=r"function $y\,c(y)$ (m${}^2$)")         
        # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
        vertices = [(0, 0)] + list(zip(vY, vY_mac_integrand)) + [(b/2, 0)]
        poly = Polygon(vertices, facecolor="brown", alpha=0.3, edgecolor="none")
        ax0.add_patch(poly)

    if ('f_alpha0L_integral' in kwargs):
        f_alpha0L_integral = kwargs['f_alpha0L_integral']
        vAlpha0L_integrand = []
        for y_ in vY:
            vAlpha0L_integrand.append(
                f_alpha0L_integral(y_, c_r, c_k, c_t, 
                                   eps_k, eps_t, alpha0l_r, alpha0l_k, alpha0l_t, 
                                   b_k, b, Lambda_le_1, Lambda_le_2)
                )
        vAlpha0L_integrand = np.asarray(vAlpha0L_integrand)
        plt.plot(vY, vAlpha0L_integrand*180/np.pi, linewidth=2.5, linestyle="-", dashes=[1000,1], marker="." ,
                 label=r"function $\big(\alpha_{0\ell} - \epsilon_{\mathrm{g}}\big) c$ (deg$\,$m)")         
        # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
        vertices = [(0, 0)] + list(zip(vY, vAlpha0L_integrand*180/np.pi)) + [(b/2, 0)]
        poly = Polygon(vertices, facecolor="red", alpha=0.3, edgecolor="none")
        ax0.add_patch(poly)

    # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
    #vertices = [(0, 0)] + list(zip(vY, vIntegrand*180/np.pi)) + [(b/2, 0)]
    #poly = Polygon(vertices, facecolor="orange", alpha=0.5, edgecolor="none")
    #ax0.add_patch(poly)
    
    plt.legend(loc='upper center', fontsize=18)
    
    if ('xmin' in kwargs):
        xmin = kwargs['xmin']
    else:
        xmin = 0

    if ('xmax' in kwargs):
        xmax = kwargs['xmax']
    else:
        xmax = 1.1*b/2

    if ('ymin' in kwargs):
        ymin = kwargs['ymin']
    else:
        ymin = -5

    if ('ymax' in kwargs):
        ymax = kwargs['ymax']
    else:
        ymax = 8
    
    plt.axis([xmin, xmax, ymin, ymax])
    
    # some annotations
    tipLine, = plt.plot([b/2,b/2], [0.9*ymin, 0.8*ymax], color="gray", linewidth=1.0, linestyle="-")
    tipLine.set_dashes([8, 4]) 
    plt.annotate(r'$y=\frac{1}{2}\,b$',
                 xy=(b/2, 0.65*ymin), xycoords='data',
                 xytext=(40, -40), textcoords='offset points', fontsize=22,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.5"))

    kinkLine, = plt.plot([b_k/2,b_k/2], [0.9*ymin, 0.8*ymax], color="gray", linewidth=1.0, linestyle="-")
    kinkLine.set_dashes([8, 4]) 
    plt.annotate(r'$y=\frac{1}{2}\,b_{\mathrm{k}}$',
                 xy=(b_k/2, 0.65*ymin), xycoords='data',
                 xytext=(40, -40), textcoords='offset points', fontsize=22,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.5"))
    
    zeroYLine, = plt.plot([xmin,xmax], [0,0],  color="gray", linewidth=1.0, linestyle="-")
    zeroYLine.set_dashes([8, 4]) 
    
    plt.title('Wing functions', fontsize=22)
    plt.xlabel('$y$ (m)', fontsize=22)
    #plt.ylabel('$X$ (m)', fontsize=22)
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    vshift_xaxis = kwargs.get('vshift_xaxis', None)
    if ('vshift_xaxis' in kwargs):
        vshift_xaxis = kwargs['vshift_xaxis']
    else:
        vshift_xaxis = 10

    ax.spines['bottom'].set_position(('outward',vshift_xaxis))
    
    hshift_yaxis = kwargs.get('hshift_yaxis', None)
    if ('hshift_yaxis' in kwargs):
        hshift_yaxis = kwargs['hshift_yaxis']
    else:
        hshift_yaxis = 10
        
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('outward',hshift_yaxis))    

    plt.show()

#----------------------------------------------------------
def retrieve_stored_data(file_name):
    store = shelve.open(file_name, flag='r')
    
    #print('-------- Database content, key => value --------')
    #for key in store:
    #    print(key, '=> {0}'.format(store[key]))
    
    #store.close()
    return store # not closed!

#----------------------------------------------------------
def data_summary(store):

    #print('-------- Database content, key => value --------')
    #for key in store:
    #    print(key, '=> {0}'.format(store[key]))

    c_r = store['c_r']
    c_k = store['c_k']
    c_t = store['c_t']
    eps_k = store['eps_k']
    eps_t = store['eps_t']
    alpha0l_r = store['alpha0l_r']
    alpha0l_k = store['alpha0l_k']
    alpha0l_t = store['alpha0l_t']
    b_k = store['b_k']
    b = store['b']
    Lambda_le_1 = store['Lambda_le_1']
    Lambda_le_2 = store['Lambda_le_2']
    # f_chord = store['f_chord']
    S_ref = store['S_ref']
    c_mac = store['c_mac']
    X_le_mac = store['X_le_mac']
    Y_mac = store['Y_mac']
    alpha0L = store['alpha0L']
    return Latex(
        r'\begin{array}{rl}'
        +  r'\text{root chord,}\, c_{\mathrm{r}}: & ' + r'{0}'.format(c_r) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{kink chord,}\, c_{\mathrm{k}}: & ' + r'{0}'.format(c_k) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{tip chord,}\, c_{\mathrm{t}}: & ' + r'{0}'.format(c_t) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{semispan, inner panel}\, \frac{1}{2}b_{\mathrm{k}}: & ' + r'{0}'.format(b_k/2) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{semispan,}\, \frac{1}{2}b: & ' + r'{0}'.format(b/2) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{leading edge sweep, inner panel,}\, \Lambda_{\mathrm{le},1}: &' 
        +    r'{0}'.format(Lambda_le_1*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{leading edge sweep, outer panel,}\, \Lambda_{\mathrm{le},2}: &' 
        +    r'{0}'.format(Lambda_le_2*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{kink section geometric twist,}\, \epsilon_{\mathrm{g,k}}: & ' + r'{0:.4}'.format(eps_k*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{tip section geometric twist,}\, \epsilon_{\mathrm{g,t}}: & ' + r'{0:.4}'.format(eps_t*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{root profile zero-lift angle of attack,}\, \alpha_{0\ell,\mathrm{r}}: & ' + r'{0:.4}'.format(alpha0l_r*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{kink profile zero-lift angle of attack,}\, \alpha_{0\ell,\mathrm{k}}: & ' + r'{0:.4}'.format(alpha0l_k*180/math.pi) + r'\,\text{deg}'
        +  r'\\'
        +  r'\text{tip profile zero-lift angle of attack,}\, \alpha_{0\ell,\mathrm{t}}: & ' + r'{0:.4}'.format(alpha0l_t*180/math.pi) + r'\,\text{deg}'
        +  r'\\[1em] \hline'
        +  r'\text{mean aerodynamic chord (MAC),}\, \bar{c}: & ' + r'{0:.4}'.format(c_mac) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{MAC leading edge longitudinal position,}\, X_{\mathrm{le},\bar{c}}: & ' + r'{0:.4}'.format(X_le_mac) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{MAC leading edge spanwise position,}\, Y_{\bar{c}}: & ' + r'{0:.4}'.format(Y_mac) + r'\,\text{m}'
        +  r'\\'
        +  r'\text{wing zero-lift angle of attack,}\, \alpha_{0L,\mathrm{W}}: & ' + r'{0:.3}'.format(alpha0L*180/math.pi) + r'\,\text{deg}'
        +r'\end{array}'
    )
