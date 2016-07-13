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
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
colors = ['#097054','#FFDE00','#6599FF','#FF9900','#993300','#99CC99',
          '#003366','#92CD00','#FFCF79','#2C6700','#996699','#666666']
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
    ax.spines['left'].set_alpha(0.75)
    ax.spines['left'].set_linestyle('dotted')
    
    ax.spines['bottom'].set_lw(0.5)
    ax.spines['bottom'].set_color('#666666')
    ax.spines['bottom'].set_alpha(0.75)

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
# External Data
HiCFolder = '/data1/xtwang/workspace/new_hic_pipeline/Corrected-hg19'
sigFolder = '/data1/xtwang/workspace/ENCODE'
TADFolder = '/data1/xtwang/workspace/calTADs'
Cell = 'GM12878'
enzyme = 'MboI'
ResHiC = 5000
ResLabel = str(ResHiC//1000) + 'K_c'
Pre = '-'.join([Cell, enzyme, 'allReps-filtered', ResLabel])
HiCFil = Pre + '-sparse.npz'
TADFil = '.'.join([Pre, 'hmmdomain', 'txt'])
HiCSource = os.path.join(HiCFolder, HiCFil)
TADSource = os.path.join(TADFolder, TADFil)
HiCData = np.load(HiCSource)
TADData = analyze.TAD(TADSource).data
Sigs = ['Rnaseq','CTCF']
chrom = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
interval = 1000
ResSig = 5000 # Signal Data Resolution
step = ResHiC * interval / ResSig
fold = 5000 # Figure resolution

# Output Figure Nameb
oPrefix = 'GM12878-MboI-allReps-filtered-5K_c-sparse_10M'

chroms = 'hg19.txt'
 
Type = np.dtype({'names':['chr','start1','start2'],'formats':['S2',np.int,np.int]})
scnvs = np.loadtxt('GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt',dtype=Type,usecols=(0,1,4),skiprows=1)
scnv_color = 'm'
#--------------------------------------------------------------------------
# Signals
sigDest = os.path.join(sigFolder, Cell)
sigFils = os.listdir(sigDest)
database = {}
for f in sigFils:
    if f.endswith('.stat'):
        parse = f.split('_')[1]
        if parse in Sigs:
            database[parse] = np.loadtxt(os.path.join(sigDest, f))
chromlens = np.loadtxt(chroms, dtype = np.int)
cumlens = list(np.cumsum(chromlens//ResSig))
cumlens.insert(0, 0)
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
    # Hi-C
    inter = HiCData[i]
    # Signals
    if i == 'X':
        itoN = '23'
    else:
        itoN = i
    cidx = int(itoN) - 1
    current = {}
    for name in database:
        current[name] = database[name][cumlens[cidx]:cumlens[cidx+1]]
    # TADs
    tads = TADData[TADData['chr'] == i]
    chr_scnv = scnvs[scnvs['chr'] == i]
    # Use PDF Backend
    pp = PdfPages(oPrefix + '-' + str(i) + '.pdf')
    ## One figure for each genomic interval
    startHiC = 0
    startSig = 0
    for idx in range(max(inter['bin1'].max(), inter['bin2'].max()) // interval):
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left, HB, width, HH])
        EndHiC = startHiC + interval
        ## Heatmap
        Matrix = getmatrix(inter, startHiC, EndHiC)
        nonzero = Matrix[np.nonzero(Matrix)]
        if nonzero.size <= 10:
            plt.close(fig)
            startHiC = EndHiC
            startSig += step
            continue
        vmax = np.percentile(nonzero, 95)
        sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        mask_scnvs = (chr_scnv['start1']//ResHiC >= startHiC) & (chr_scnv['start2']//ResHiC < EndHiC)
        extract_scnvs = chr_scnv[mask_scnvs]
        scnv_plotx = extract_scnvs['start1']//ResHiC - startHiC + 0.5
        scnv_ploty = extract_scnvs['start2']//ResHiC - startHiC + 0.5
##        ax.scatter(scnv_plotx, scnv_ploty, marker = '.', edgecolor = 'none', c = scnv_color, s = 18)
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
        pairs = [(i['start']//ResHiC - startHiC, i['end']//ResHiC - startHiC)
                 for i in extract]
        for corner in pairs:
            if (corner[0] <= 0):
                ax.plot([0, corner[1]], [corner[1], corner[1]], color = '#666666',
                        linewidth = 1)
                ax.plot([corner[1], corner[1]], [0, corner[1]], color = '#666666',
                        linewidth = 1)
            elif (corner[1] >= interval):
                ax.plot([corner[0], corner[0]], [corner[0], interval],
                        color = '#666666', linewidth = 1)
                ax.plot([corner[0], interval], [corner[0], corner[0]],
                        color = '#666666', linewidth = 1)
            else:
                ax.plot([corner[0], corner[0]], [corner[0], corner[1]],
                        color = '#666666', linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[0], corner[0]],
                        color = '#666666', linewidth = 1)
                ax.plot([corner[0], corner[1]], [corner[1], corner[1]],
                        color = '#666666', linewidth = 1)
                ax.plot([corner[1], corner[1]], [corner[0], corner[1]],
                        color = '#666666', linewidth = 1)
        
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        caxis_H(ax)
        
        ## Colorbar
        ax = fig.add_axes([Left + width + 0.08, HB, 0.03, HH])
        fig.colorbar(sc, cax = ax, format = '%d')
        
        startHiC = EndHiC # Reset the start site
        
        # Other One-dimensional information
        EndSig = startSig + step
        cb = SB # Current Bottom
        count = 0 # Used for color index
        bounds = [(i['start']//fold - startSig//(fold//ResSig), i['end']//fold - startSig//(fold//ResSig)) for i in extract]
        for name in current:
            sig = current[name][startSig:EndSig]
            newS = convert(sig, ResSig, fold)
            ax = fig.add_axes([Left, cb, width, SH])
            ax.plot(newS, color = colors[count])
            ax.set_ylabel(name, labelpad = 70, rotation = 'horizontal',
                          style = 'italic', size = 10)
                              
            ylim = ax.get_ylim()
                
            ## Domain Boundarys
#            for bi in range(len(bounds)):
#                if bi == 0:
#                    if bounds[bi][0] > 0:
#                        ax.vlines(bounds[bi][0], ylim[0], ylim[1], colors = '#666666')
#                else:
#                    if bounds[bi][0] != bounds[bi-1][1]:
#                        ax.vlines(bounds[bi][0], ylim[0], ylim[1], colors = '#666666')
#                if bounds[bi][1] < step//(fold//ResSig):
#                    ax.vlines(bounds[bi][1], ylim[0], ylim[1], colors = '#666666')
                
            # Axis Control
            caxis_S(ax, colors[count])
                        
            ## Reset arguments
            cb += SH; count += 1
        
        startSig = EndSig
        
        
        pp.savefig(fig)
        
        plt.close(fig)
        
    pp.close()
