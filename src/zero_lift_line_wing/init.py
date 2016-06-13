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

from sympy import *

from IPython.display import display, Math, Latex, SVG

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('pdf', 'png')
plt.rcParams['savefig.dpi'] = 75

plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.figsize'] = 10, 6
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

#----------------------------------------------------------
def plot_planform(c_r, c_t, b):
    xLineWing = [0,b/2,b/2,0]
    yLineWing = [0,0.25*c_r-0.25*c_t,0.25*c_r+0.75*c_t,c_r]

    # planform
    lineWing, = plt.plot(xLineWing, yLineWing, 'k-')
    # centerline
    centerLine, = plt.plot([0,0], [-1.1*c_r,2.1*c_r], 'b')
    centerLine.set_dashes([8, 4, 2, 4]) 
    # c/4 line
    quarterChordLine, = plt.plot([0,1.05*b], [0.25*c_r,0.25*c_r], 'k--')

    plt.axis('equal')
    plt.axis([-0.1*b/2, 1.1*b/2, -0.1*c_r, 1.1*c_r])
    plt.gca().invert_yaxis()
    plt.title('Wing planform', fontsize=22)
    plt.xlabel('$y$ (m)', fontsize=22)
    plt.ylabel('$X$ (m)', fontsize=22)
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',2*c_r))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.1*b/2))

    plt.show()

#----------------------------------------------------------
def plot_wing_functions(c_r, c_t, b, A_c, B_c, A_alpha, B_alpha, A_eps):
    """
    
    See: http://www.scipy-lectures.org/intro/matplotlib/matplotlib.html
    
    """
    
    # define vectors
    vY = np.linspace(0, b/2, 10, endpoint=True)
    vChord = A_c*vY + B_c
    vAlpha0L = A_alpha*vY + B_alpha
    vEpsilon = A_eps*vY
    vIntegrand = (vAlpha0L - vEpsilon)*vChord
    
    # Create a figure of size WxH inches, DPI dots per inch
    fig = plt.figure(figsize=(9, 8), dpi=300)
    
    # Create a new subplot from a grid of 1x1
    ax0 = plt.subplot(1, 1, 1)
    
    #fig, ax0 = plt.subplots()
    
    
    plt.plot(vY, vChord, color="red", linewidth=2.5, linestyle="-", label="local chord (m)")
    plt.plot(vY, vAlpha0L*180/np.pi, color="green",  linewidth=2.5, linestyle="-", label=r"local $\alpha_{0\ell}$ (deg)")
    plt.plot(vY, vEpsilon*180/np.pi, color="blue",  linewidth=2.5, linestyle="-", label=r"local $\epsilon_{\mathrm{g}}$ (deg)")
    plt.plot(vY, vIntegrand*180/np.pi, color="brown",  linewidth=2.5, linestyle="-", 
             label=r"Integrand function $\big(\alpha_{0\ell}$ - \epsilon_{\mathrm{g}}\big) c$ (deg m)")
    
    # shaded region --> http://matplotlib.org/examples/showcase/integral_demo.html
    vertices = [(0, 0)] + list(zip(vY, vIntegrand*180/np.pi)) + [(b/2, 0)]
    poly = Polygon(vertices, facecolor="orange", alpha=0.5, edgecolor="none")
    ax0.add_patch(poly)
    
    plt.legend(loc='upper center', fontsize=18)
    
    tipLine, = plt.plot([b/2,b/2], [-4, 2], color="gray", linewidth=1.0)
    tipLine.set_dashes([8, 4]) 
    plt.annotate(r'$y=\frac{b}{2}$',
                 xy=(b/2, -4), xycoords='data',
                 xytext=(40, -40), textcoords='offset points', fontsize=22,
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.5"))
    
    zeroYLine, = plt.plot([-1,b/2], [0,0],  color="gray", linewidth=1.0)
    zeroYLine.set_dashes([8, 4]) 
    
    plt.axis([0, 1.1*b/2, -5, 6])
    plt.title('Wing functions', fontsize=22)
    plt.xlabel('$y$ (m)', fontsize=22)
    #plt.ylabel('$X$ (m)', fontsize=22)
    
    # Moving spines
    ax = plt.gca()  # gca stands for 'get current axis'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',-6))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',-0.1))    
    plt.show()
    
#----------------------------------------------------------
def display_workflow_S(b, A_c, B_c, c_law_integral_indefinite, c_law_integral_definite):
    c_law_integral_indefinite_latex = latex(c_law_integral_indefinite)
    return Latex(
        r'\begin{align*}'
            r'S & {}= 2 \int_0^{' + '{0:.3}'.format(b/2) + '}' 
            + r'\Big(' 
            +   r'{0:.4}'.format(A_c) + r'\, y + ' + '{0:.4}'.format(B_c) + r'\,\text{m}'
            + r'\Big)\,\text{d}y' 
            + r'\\[1em]' 
            + r'& {}= 2 \big(' 
            +   c_law_integral_indefinite_latex 
            + r'\big)\Bigr|_0^{' + '{0:.4}'.format(b/2) + r'} \,\text{m}^2'
            + '=' + '{0:.4}'.format(c_law_integral_definite) + r'\,\text{m}^2'
        + r'\end{align*}'
    )
    
#----------------------------------------------------------
def display_workflow_alpha0L(b, S_ref, A_c, B_c, A_alpha, B_alpha, A_eps,
                             alpha0L_law_integral_indefinite,
                             alpha0L):
    alpha0L_law_integral_indefinite_latex = latex(alpha0L_law_integral_indefinite)
    return Latex(
        r'\begin{multline*}'
            + r'\alpha_{0L,\mathrm{W}} = \frac{2}{' + '{0:.4}'.format(S_ref) + r'\,\text{m}^2}'
            + r'\int_0^{' + '{0:.3}'.format(b/2) + r'}' 
            + r'\Big(' 
            +   r'{0:.4}'.format(A_alpha) + r'\, \frac{\text{rad}}{\text{m}}\,y ' 
            +     '{0:.4}'.format(B_alpha) + r'\,\text{rad}'
            +     '{0:.4}'.format(A_eps) + r'\,\frac{\text{rad}}{\text{m}}y'            
            + r'\Big)\,'
            + r'\Big(' 
            +   r'{0:.4}'.format(A_c) + r'\, y + ' + '{0:.4}'.format(B_c) + r'\,\text{m}'
            + r'\Big)\,'
            + r'\text{d}y' 
            + r'\\[1em]' 
            + r'= \frac{2}{' + '{0:.4}'.format(S_ref) + r'} \big(' 
            +   alpha0L_law_integral_indefinite_latex 
            + r'\big)\Bigr|_0^{' + '{0:.4}'.format(b/2) + r'} \,\text{rad}'
            + r'\\[1em]' 
            + r'= {0:.4}'.format(alpha0L) + r'\,\text{rad}'
            + r'= {0:.4}'.format(alpha0L*180/math.pi) + r'\,\text{deg}'
        + r'\end{multline*}'
        )