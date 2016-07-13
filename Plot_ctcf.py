# Created on Mon Jan 12 20:09:28 2015

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

## Necessary Modules
#--------------------------------------------------------------------------
from __future__ import division
import numpy as np
from tadlib.analyze import getmatrix
import matplotlib, collections
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
# Color Map for CTCFs
from palettable.colorbrewer.sequential import YlGnBu_3, YlOrRd_3
CTCF_colors = YlGnBu_3.get_mpl_colormap()
peak_colors = YlOrRd_3.get_mpl_colormap()
#--------------------------------------------------------------------------
## Utilities
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color('k')
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
#--------------------------------------------------------------------------
# Hi-C Related ...
HiCFolder = '/data3/3dgenome/data/HiC/Corrected-hg19'
sigFolder = '/public/home/3dgenome/wenzi/Plot_ctcf'
TADFolder = '/public/home/3dgenome/wenzi/Plot_ctcf'
Cell = 'GM12878'
enzyme = 'MboI'
ResHiC = 10000
ResLabel = str(ResHiC//1000) + 'K_c'
Pre = '-'.join([Cell, enzyme, 'allReps-filtered', ResLabel])
HiCFil = Pre + '-sparse.npz'
TADFil = '.'.join([Pre, 'hmmdomain', 'txt'])
HiCSource = os.path.join(HiCFolder, HiCFil)
TADSource = os.path.join(TADFolder, 'GM12878-MboI-hierDomains-chr1.txt')
HiCData = np.load(HiCSource)
TADData = analyze.TAD(TADSource).data
chrom = ['1']
interval = 400
# Output Figure Name
oPrefix = 'GM12878-MboI-allReps-filtered-10K_c-sparse'
chroms = 'hg19.txt'
#--------------------------------------------------------------------------
# CTCF Motif Matches
CTCF_dtype = np.dtype({'names':['chr','start','end','strength','strand'],
                       'formats':['S5', np.int, np.int, np.float, 'S1']})
CTCF_data = np.loadtxt('gm12878_100.gff', dtype = CTCF_dtype,
                       usecols = [0, 3, 4, 5, 6])
for i in CTCF_data:
    i['strand'] = '>' if i['strand'] == '-' else '<'

parse_CTCF = [line.rstrip().split()[-1] for line in open('gm12878_100.gff') if not line.startswith('#')]
parse_q = np.array([float(mem[(mem.find('qvalue')+7):].split(';')[0]) for mem in parse_CTCF])
CTCF_data = CTCF_data[parse_q <= 0.184]
# CTCF / Rad21 Peaks
peak_dtype = np.dtype({'names':['chr','pos','sig','submit'],
                       'formats':['S5', np.int, np.float, np.int]})
CTCF_peak = np.loadtxt('wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak',
                       dtype = peak_dtype, usecols = [0, 1, 6, -1])
Rad21_peak = np.loadtxt('wgEncodeAwgTfbsHaibGm12878Rad21V0416101UniPk.narrowPeak',
                        dtype = peak_dtype, usecols = [0, 1, 6, -1])
CTCF_peak['pos'] += CTCF_peak['submit']
Rad21_peak['pos'] += Rad21_peak['submit']
#--------------------------------------------------------------------------
# About the figure
size = (12, 16)
width = 0.618; Left = (1 - width) / 2
HB = 0.1; HH = width * size[0] / size[1] # HeatMap Bottom / HeatMap Height
SB = HB + HH # Bottom of tracks
ST = 0.9 # The Top of tracks
SSH = (ST - SB) / 9 # Standard height of each track
#--------------------------------------------------------------------------
# Draw
for i in chrom:
    ## Data by Chromosome
    inter = HiCData[i]
    Current_CTCF = CTCF_data[CTCF_data['chr'] == ('chr' + i)]
    Current_cpeak = CTCF_peak[CTCF_peak['chr'] == ('chr' + i)]
    Current_rpeak = Rad21_peak[Rad21_peak['chr'] == ('chr' + i)]
    tads = TADData[TADData['chr'] == i]
    Color_norm1 = matplotlib.colors.Normalize(vmin = Current_CTCF['strength'].min(),
                                              vmax = Current_CTCF['strength'].max())
    Color_norm2 = matplotlib.colors.Normalize(vmin = Current_cpeak['sig'].min(),
                                              vmax = Current_cpeak['sig'].max())
    Color_norm3 = matplotlib.colors.Normalize(vmin = Current_rpeak['sig'].min(),
                                              vmax = Current_rpeak['sig'].max())
    startBins1 = Current_CTCF['start'] // ResHiC
    startBins2 = Current_cpeak['pos'] // ResHiC
    startBins3 = Current_rpeak['pos'] // ResHiC
    # Use PDF Backend
    pp = PdfPages(oPrefix + '-' + str(i) + '.pdf')
    ## One figure for each genomic interval
    startHiC = 0
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
            continue
        vmax = np.percentile(nonzero, 95)
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
        
        CTCF_Mask = (startBins1 >= startHiC) & (startBins1 < EndHiC)
        cpeak_Mask = (startBins2 >= startHiC) & (startBins2 < EndHiC)
        rpeak_Mask = (startBins3 >= startHiC) & (startBins3 < EndHiC)
        startBP = startHiC * ResHiC
        endBP = EndHiC * ResHiC
        startHiC = EndHiC # Reset the start site
        
        cpeak_Extract = Current_cpeak[cpeak_Mask]
        rpeak_Extract = Current_rpeak[rpeak_Mask]
        
        # Plot directional CTCF matches ...
        CTCF_Extract = Current_CTCF[CTCF_Mask]
        startBins_Extract = startBins1[CTCF_Mask]
        CTCF_collections = collections.Counter(startBins_Extract)
        newStruct = {}
        for ci in CTCF_collections:
            newStruct[ci] = CTCF_Extract[startBins_Extract==ci]
            newStruct[ci].sort(order = ['start'])
        
        most_common = CTCF_collections.most_common(1)
        
        bounds = [(bi['start'] - startBP, bi['end'] - startBP) for bi in extract]
        if len(most_common) > 0:
            MSH = (ST - SB) / most_common[0][1]
            SH = min(SSH, MSH)
            cb = SB
            if len(cpeak_Extract) > 0:
                ax_temp1 = fig.add_axes([Left, cb, width, SH])
                for ci in cpeak_Extract:
                    local_pos = ci['pos'] - startBP
                    local_color = peak_colors(Color_norm2(ci['sig']))
                    ax_temp1.vlines(local_pos, 0.4, 0.6, colors = local_color, linewidth = 3)
                cb += SH
                ax_temp1.set_xlim(0, endBP - startBP)
                ax_temp1.set_ylim(0, 1)
                ax_temp1.set_ylabel('CTCFPeak', labelpad = 45, rotation = 'horizontal',
                                    style = 'italic', size = 10)
                
                for bi in range(len(bounds)):
                    if bi == 0:
                        if bounds[bi][0] > 0:
                            ax_temp1.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    else:
                        if bounds[bi][0] != bounds[bi-1][1]:
                            ax_temp1.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    if bounds[bi][1] < endBP:
                        ax_temp1.vlines(bounds[bi][1], 0, 1, colors = '#666666', alpha = 0.5)
                        
                caxis_S(ax_temp1)
                
            if len(rpeak_Extract) > 0:
                ax_temp2 = fig.add_axes([Left, cb, width, SH])
                for ri in rpeak_Extract:
                    local_pos = ri['pos'] - startBP
                    local_color = peak_colors(Color_norm3(ri['sig']))
                    ax_temp2.vlines(local_pos, 0.4, 0.6, colors = local_color, linewidth = 3)
                cb += SH
                ax_temp2.set_xlim(0, endBP - startBP)
                ax_temp2.set_ylim(0, 1)
                ax_temp2.set_ylabel('Rad21Peak', labelpad = 45, rotation = 'horizontal',
                                    style = 'italic', size = 10)
                
                for bi in range(len(bounds)):
                    if bi == 0:
                        if bounds[bi][0] > 0:
                            ax_temp2.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    else:
                        if bounds[bi][0] != bounds[bi-1][1]:
                            ax_temp2.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    if bounds[bi][1] < endBP:
                        ax_temp2.vlines(bounds[bi][1], 0, 1, colors = '#666666', alpha = 0.5)
                        
                caxis_S(ax_temp2)
                
            track_num = most_common[0][1]
            
            ax_list = []
            for ti in range(track_num):
                ax_list.append(fig.add_axes([Left, cb, width, SH]))
                cb += SH
            
            for ci in newStruct:
                CTCF_count = len(newStruct[ci])
                for mi in range(CTCF_count):
                    local_ax = ax_list[mi]
                    local_start = newStruct[ci][mi]['start'] - startBP
                    local_color = CTCF_colors(Color_norm1(newStruct[ci][mi]['strength']))
                    local_ax.plot(local_start, 0.5, newStruct[ci][mi]['strand'],
                                  markersize = 4, markeredgecolor = 'none',
                                  markerfacecolor = local_color)
                    local_ax.set_xlim(0, endBP - startBP)
            
            for ti in range(track_num):
                local_ax = ax_list[ti]
                ylim = local_ax.get_ylim()
                for bi in range(len(bounds)):
                    if bi == 0:
                        if bounds[bi][0] > 0:
                            local_ax.vlines(bounds[bi][0], ylim[0], ylim[1], colors = '#666666', alpha = 0.5)
                    else:
                        if bounds[bi][0] != bounds[bi-1][1]:
                            local_ax.vlines(bounds[bi][0], ylim[0], ylim[1], colors = '#666666', alpha = 0.5)
                    if bounds[bi][1] < endBP:
                        local_ax.vlines(bounds[bi][1], ylim[0], ylim[1], colors = '#666666', alpha = 0.5)
                        
                caxis_S(local_ax)
        else:
            SH = SSH
            cb = SB
            if len(cpeak_Extract) > 0:
                ax_temp1 = fig.add_axes([Left, cb, width, SH])
                for ci in cpeak_Extract:
                    local_pos = ci['pos'] - startBP
                    local_color = peak_colors(Color_norm2(ci['sig']))
                    ax_temp1.vlines(local_pos, 0.4, 0.6, colors = local_color, linewidth = 3)
                cb += SH
                ax_temp1.set_xlim(0, endBP - startBP)
                ax_temp1.set_ylim(0, 1)
                ax_temp1.set_ylabel('CTCFPeak', labelpad = 45, rotation = 'horizontal',
                                    style = 'italic', size = 10)
                for bi in range(len(bounds)):
                    if bi == 0:
                        if bounds[bi][0] > 0:
                            ax_temp1.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    else:
                        if bounds[bi][0] != bounds[bi-1][1]:
                            ax_temp1.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    if bounds[bi][1] < endBP:
                        ax_temp1.vlines(bounds[bi][1], 0, 1, colors = '#666666', alpha = 0.5)
                caxis_S(ax_temp1)
                
            if len(rpeak_Extract) > 0:
                ax_temp2 = fig.add_axes([Left, cb, width, SH])
                for ri in rpeak_Extract:
                    local_pos = ri['pos'] - startBP
                    local_color = peak_colors(Color_norm3(ri['sig']))
                    ax_temp2.vlines(local_pos, 0.4, 0.6, colors = local_color, linewidth = 3)
                cb += SH
                ax_temp2.set_xlim(0, endBP - startBP)
                ax_temp2.set_ylim(0, 1)
                ax_temp2.set_ylabel('Rad21Peak', labelpad = 45, rotation = 'horizontal',
                                    style = 'italic', size = 10)
                for bi in range(len(bounds)):
                    if bi == 0:
                        if bounds[bi][0] > 0:
                            ax_temp2.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    else:
                        if bounds[bi][0] != bounds[bi-1][1]:
                            ax_temp2.vlines(bounds[bi][0], 0, 1, colors = '#666666', alpha = 0.5)
                    if bounds[bi][1] < endBP:
                        ax_temp2.vlines(bounds[bi][1], 0, 1, colors = '#666666', alpha = 0.5)
                caxis_S(ax_temp2)
        
        pp.savefig(fig)
        
        plt.close(fig)
        
    pp.close()
