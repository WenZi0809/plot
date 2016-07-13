# Created on Mon Jan 12 20:09:28 2015

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

## Necessary Modules
#--------------------------------------------------------------------------
from __future__ import division
import numpy as np
from tadlib.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from tadlib import analyze
from palettable.colorbrewer.qualitative import Dark2_8
hexcolors = Dark2_8.hex_colors
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
#--------------------------------------------------------------------------
## Utilities
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0.9)
    ax.spines['left'].set_linestyle('dotted')
    
    ax.spines['bottom'].set_lw(0.5)
    ax.spines['bottom'].set_color('#b3b3b3')
    ax.spines['bottom'].set_alpha(0.9)

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])

def convert(sig, Res, fold):
    """
    Signal conversion for lower-resolution plotting.
    
    """
    x = np.arange(sig.size) // (fold // Res)
    w = sig
    
    return np.bincount(x, weights = w)
#--------------------------------------------------------------------------
# Hi-C Related ...
HiCFolder = '/data3/3dgenome/data/HiC/Corrected-hg19'
TADFolder = '/public/home/3dgenome/wenzi/CTCF/gm12878'
loopFolder = '/public/home/3dgenome/wenzi/CTCF/gm12878'
Cell = 'GM12878'
enzyme = 'MboI'
ResHiC = 5000
ResLabel = str(ResHiC//1000) + 'K_c'
Pre = '-'.join([Cell, enzyme, 'allReps-filtered', ResLabel])
HiCFil = Pre + '-sparse.npz'
#TADFil = '.'.join([Pre, 'hmmdomain', 'txt'])
HiCSource = os.path.join(HiCFolder, HiCFil)
TADSource = os.path.join(TADFolder, 'GM12878-MboI-allReps-filtered-5K_c.arrowhead.txt')
HiCData = np.load(HiCSource)
TADData = analyze.TAD(TADSource).data
chrom = ['22']
interval = 400
# Output Figure Name
oPrefix = 'GM12878-HeatmapsSignals-NEW_CTCF_C-Fos'
#--------------------------------------------------------------------------
# One-dimensional Tracks
sigFolder = '/public/home/3dgenome/wenzi/CTCF/gm12878/22'
Sigs =  ['Ctcf','C-fos']
#Sigs = ['RNAPII', 'H3K4me3', 'H3K27ac', 'H3K4me1']
ResSig = 5000
step = ResHiC * interval / ResSig
fold = 5000
database = {}
for s in Sigs:
    desFil = Cell + '_' + s + '.npy'
    database[s] = np.load(os.path.join(sigFolder, desFil))
#--------------------------------------------------------------------------
# Read Loop Data
loopType = np.dtype({'names':['chr','S1','E1','S2','E2'],
                     'formats':['S2', np.int, np.int, np.int, np.int]})
loopFil = 'GM12878-MboI-allReps-filtered-5K_c.loops.txt'
loopSource = os.path.join(loopFolder, loopFil)
loops = np.loadtxt(loopSource, dtype = loopType)
#--------------------------------------------------------------------------
# About the figure
size = (12, 16)
width = 0.618; Left = (1 - width) / 2
HB = 0.1; HH = width * size[0] / size[1] # HeatMap Bottom / HeatMap Height
SB = HB + HH # Bottom of tracks
ST = 0.9 # The Top of tracks
SH = (ST - SB) / (len(Sigs)) # Height of each track
#--------------------------------------------------------------------------
# Draw
for i in chrom:
    ## Data by Chromosome
    inter = HiCData[i]
    tads = TADData[TADData['chr'] == i]
    loopbychrom = loops[loops['chr'] == i]
    current = {}
    for name in Sigs:
        current[name] = database[name]['value'][database[name]['chr'] == i]
    # Use PDF Backend
    pp = PdfPages(oPrefix + '-' + str(i) + '.pdf')
    ## One figure for each genomic interval
    startHiC = 0
    startSig = 0
    for idx in range(max(inter['bin1'].max(), inter['bin2'].max()) // interval):
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left, HB, width, HH])
        EndHiC = startHiC + interval
        ## Get Contact Matrix
        Matrix = getmatrix(inter, startHiC, EndHiC)
        nonzero = Matrix[np.nonzero(Matrix)]
        if nonzero.size <= 100:
            plt.close(fig)
            startHiC = EndHiC
            startSig += step
            continue
        vmax = np.percentile(nonzero, 95)
        ## Mask original matrix using loop data
        mask = (loopbychrom['S1'] >= startHiC * ResHiC) & (loopbychrom['E2'] < EndHiC * ResHiC)
        extract = loopbychrom[mask]
        Bool = np.zeros_like(Matrix, dtype = bool)
        for point in extract:
            xi = point['S1'] // ResHiC - startHiC
            xj = point['E1'] // ResHiC - startHiC
            yi = point['S2'] // ResHiC - startHiC
            yj = point['E2'] // ResHiC - startHiC
            Bool[xi:xj, yi:yj] = 1
        Matrix = np.ma.array(Matrix, mask = Bool)
        ## Heatmap Plotting
        sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0, interval, 6).astype(int))
        pos = [((startHiC + t) * ResHiC) for t in ticks]
        labels = [properU(p) for p in pos]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ## Domain Boundaries
        mask = (tads['end'] > startHiC * ResHiC) & (tads['start'] < EndHiC * ResHiC)
        extract = tads[mask]
        pairs = [(bi['start']//ResHiC - startHiC, bi['end']//ResHiC - startHiC)
                 for bi in extract]
        for corner in pairs:
            if (corner[0] <= 0):
                ax.plot([0, corner[1]], [corner[1], corner[1]], color = '#b3b3b3',
                        linewidth = 1)
                ax.plot([corner[1], corner[1]], [0, corner[1]], color = '#b3b3b3',
                        linewidth = 1)
            elif (corner[1] >= interval):
                ax.plot([corner[0], corner[0]], [corner[0], interval],
                        color = '#b3b3b3', linewidth = 1)
                ax.plot([corner[0], interval], [corner[0], corner[0]],
                        color = '#b3b3b3', linewidth = 1)
            else:
                ax.plot([corner[0], corner[0]], [corner[0], corner[1]],
                        color = '#b3b3b3', linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[0], corner[0]],
                        color = '#b3b3b3', linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[1], corner[1]],
                        color = '#b3b3b3', linewidth = 1)
                ax.plot([corner[1], corner[1]], [corner[0], corner[1]],
                        color = '#b3b3b3', linewidth = 1)
        
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        caxis_H(ax)
        
        ## Colorbar
        ax = fig.add_axes([Left + width + 0.08, HB, 0.03, HH])
        fig.colorbar(sc, cax = ax, format = '%d')
        
        startHiC = EndHiC
        
        ## Signal Tracks
        EndSig = startSig + step
        cb = SB
        count = 0 # Used for color index
        bounds = [(e['start']//fold - startSig//(fold//ResSig), e['end']//fold - startSig//(fold//ResSig)) for e in extract]
        for name in Sigs:
            sig = current[name][startSig:EndSig]
            newS = convert(sig, ResSig, fold)
            ax = fig.add_axes([Left, cb, width, SH])
            ax.fill_between(np.arange(newS.size), newS, color = hexcolors[count],
                            edgecolor = 'none')
            ax.set_ylabel(name, labelpad = 50, rotation = 'horizontal',
                          style = 'italic', size = 10)
                              
            ylim = ax.get_ylim()
            
            # Axis Control
            caxis_S(ax, hexcolors[count])
                        
            ## Reset arguments
            cb += SH; count += 1
        
        startSig = EndSig
        
        pp.savefig(fig)
        
        plt.close(fig)
        
    pp.close()
