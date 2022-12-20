# quick visualize functions
from _plotly_utils.colors.sequential import Viridis
#from . import compartment
from . import structure
import pandas as pd
import numpy as np

def plotPsCurve(mcoolsPath:list,celltypeNames:list,chroms:list,resolution=100000,title="P(s) curve",plotType="interaction",base=1.1,log_x=True,log_y=True):
    """
    plotPsCurve function take bin and value
    """
    import plotly.express as px
    from IPython.display import Image

    #Calculate P(s) data, get a 3 column pd.DataFrame with (bin,resolution,celltype)
    psDataAll = []
    for i in range(len(mcoolsPath)):
        psDataAll.append(compartment.getPsData(mcoolsPath[i],["chr"+str(i+1) for i in range(len(chroms))],resolution=resolution,celltype=celltypeNames[i],base=base)) 
    merged = pd.concat(psDataAll)

    data =  pd.merge(merged,merged.groupby("celltype").sum(),how="left",on="celltype").assign(prob= lambda df: df.aveCount_x/df.aveCount_y)

    fig = px.line(x=data["bin_x"]*resolution,y=data["prob"],color=data["celltype"],title=title,log_x=log_x,log_y=log_y).update_layout(template='simple_white')
    fig.update_layout(width=800,height=600)
    fig.update_layout(xaxis_title="Genomic Distance(bp)",
                    yaxis_title="Contact Probability")
    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def plotMatrix(matrix:np.ndarray,if_log=False,title="Matrix",plotType="static",range_color=None,
                axis2genome=False,genome_coord1=None,genome_coord2=None,resolution = None):
    """
    plotMatrix function for plot hic contact matrix
    """
    import plotly.express as px 
    from IPython.display import Image

    if(if_log == True): 
        matrix = np.log10(matrix)

    fig = px.imshow(matrix,color_continuous_scale=px.colors.sequential.Viridis,range_color=range_color)
    fig = fig.update_layout(template='simple_white').update_layout(width=650,height=600)

    if (axis2genome):
        import re
        #manually change axis
        posx = re.split("[:-]",genome_coord2)
        xvals = np.percentile([np.round(i) for i in range(0,matrix.shape[1])],(0,25,50,75,100),interpolation='midpoint')
        xtexts = xvals*resolution + int(posx[1].replace(",","")) + resolution/2
        xtexts = [genomecoord2human(i) for i in xtexts]

        posy = re.split("[:-]",genome_coord1)
        yvals = np.percentile([np.round(i) for i in range(0,matrix.shape[0])],(0,25,50,75,100),interpolation='midpoint')
        ytexts = yvals*resolution + int(posy[1].replace(",","")) + resolution/2
        ytexts = [genomecoord2human(i) for i in ytexts]

        fig = fig.update_xaxes(ticktext = xtexts,tickvals = xvals).update_yaxes(ticktext = ytexts,tickvals = yvals)

    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def genomecoord2human(n):
    symbols = ('bp','Kb','Mb','Gb')
    exp = int(np.log10(n)/3)
    return str(round(n/(1000**exp),2))+symbols[exp]

def plotRegionFromMCOOLS(filepath:str,resolution:int,genome_coord1:str,genome_coord2=None,
    if_log=False,balance=False,title="Matrix",plotType="static",range_color=None,color_continuous_scale = None):
    """
    plotMatrix function for plot hic contact matrix
    """
    import plotly.express as px 
    from IPython.display import Image
    import re
    import cooler
    import numpy as np

    cool = filepath+"::/resolutions/"+str(resolution)

    if (genome_coord2 == None):
        genome_coord2 = genome_coord1
    c = cooler.Cooler(cool)
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    if(if_log == True): 
        matrix = np.log10(matrix+1)
    if(color_continuous_scale == None):
        color_continuous_scale = px.colors.sequential.Viridi    
    fig = px.imshow(matrix,color_continuous_scale=color_continuous_scale,range_color=range_color)
    fig = fig.update_layout(title=title)
    fig = fig.update_layout(template='simple_white').update_layout(width=650,height=600)
    #fig = fig.update_layout(xaxis_title=genome_coord2,yaxis_title=genome_coord1)

    #manually change axis
    posx = re.split("[:-]",genome_coord2)
    xvals = np.percentile([np.round(i) for i in range(0,matrix.shape[1])],(0,25,50,75,100),interpolation='midpoint')
    xtexts = xvals*resolution + int(posx[1].replace(",","")) + resolution/2
    xtexts = [genomecoord2human(i) for i in xtexts]

    posy = re.split("[:-]",genome_coord1)
    yvals = np.percentile([np.round(i) for i in range(0,matrix.shape[0])],(0,25,50,75,100),interpolation='midpoint')
    ytexts = yvals*resolution + int(posy[1].replace(",","")) + resolution/2
    ytexts = [genomecoord2human(i) for i in ytexts]

    fig = fig.update_xaxes(ticktext = xtexts,tickvals = xvals).update_yaxes(ticktext = ytexts,tickvals = yvals)

    # static plot have better performance in jupyter
    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def plotMatrixWithGeneTrack(genes:pd.DataFrame,mat:np.ndarray, extent:list, site_bedpe=None,clow=0, chigh=1,if_colorbar=True,cmap="bwr"):
    """
    input: genes - genes to plot, should be in pd.DataFrame format contain "chrom-start-end-id-symbol-strand", sensitive to order
              mat - matrix to plot, should be in numpy.ndarray format
              extent - resize the matrix axis number to the extent size. List contain 4 elements, [x_start,x_end,y_end,y_start].
              site_bedpe - takes a standard bedpe format file and draw circles in these positions

              output: matplotlib figure
    """
    
    start1,end1,start2,end2 = extent

    from matplotlib.ticker import EngFormatter
    bp_formatter = EngFormatter('b')
    import matplotlib.pyplot as plt
    from dna_features_viewer import GraphicFeature, GraphicRecord

    def format_ticks(ax, x=True, y=True, rotate=True):
        if y:
            ax.yaxis.set_major_formatter(bp_formatter)
        if x:
            ax.xaxis.set_major_formatter(bp_formatter)
            ax.xaxis.tick_bottom()
        if rotate:
            ax.tick_params(axis='x',rotation=45)
    def symbol2num(x):
        if x=="+":
            return +1
        else: return -1

    fig, axs = plt.subplots(figsize = (9,12), gridspec_kw={'height_ratios': [1, 4]},nrows=2,
                            sharex=True,sharey=False,constrained_layout=True)

    # plot gene track
    features=[ GraphicFeature(start=gene[1], end=gene[2], strand=symbol2num(gene[5]), color="#ffffff",box_linewidth=0,
            fontdict={"size":16,"family":"sans-serif"},
            label=gene[4]) for gene in genes.values.tolist() if gene[1] >= start1 and gene[2] <= end1
    ]
    record = GraphicRecord(sequence_length=end1-start1, features=features)
    ax = axs[0]
    _ = record.plot(ax=ax)
    plt.xlim(start1,end1)
    format_ticks(ax)

    # plot hic track
    ax = axs[1]
    im2 = ax.matshow(mat ,extent=extent,vmin = clow, vmax = chigh , cmap=cmap)

    resolution = (end1-start1)/mat.shape[0]

    if site_bedpe is not None :
        for site in site_bedpe.values.tolist():
            ax.add_artist(plt.Circle(((site[1]+site[2])/2-resolution,(site[4]+site[5])/2-resolution),radius=2*resolution,color="black",alpha=0.5,fill=False))

    if(if_colorbar==True):
        plt.colorbar(im2,ax=ax,label=None)
    plt.locator_params(nbins=8)
    format_ticks(ax)
    
    
    plt.rc('font', size=15)
    plt.rc('axes', titlesize=15)     # fontsize of the axes title
    plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
    plt.rc('legend', fontsize=15)    # legend fontsize
    plt.rc('figure', titlesize=15)  # fontsize of the figure title
    
    return fig


def plotRegionWIthGeneTrack(matrix1coolPath:list,matrix2coolPath:list, genome_coord:str, gene_bed:str,site_bedpe=None,
                                                clow=-0.2,chigh=0.2,if_colorbar=False,if_diff=False,cmap="bwr"):
    import re
    import matplotlib.pyplot as plt

    coordList = re.split("[:-]",genome_coord)

    mat1 = structure.get3dProximityStackMatrix(matrix1coolPath,genome_coord)
    if(matrix1coolPath == matrix2coolPath):
        mat2 = mat1
    else:
        mat2 = structure.get3dProximityStackMatrix(matrix2coolPath,genome_coord)

    if(if_diff):
        mat = structure.get3dProximityStackMatrix(matrix1coolPath,genome_coord) - structure.get3dProximityStackMatrix(matrix2coolPath,genome_coord)
    else:
        mat = np.triu(mat1) + np.tril(mat2)
    
    genes  = pd.read_csv(gene_bed,sep="\t",header=None)
    genes  = genes[genes[0] == coordList[0]]

    fig = plotMatrixWithGeneTrack(genes = genes,mat = mat,site_bedpe=site_bedpe, 
                                            extent = [int(coordList[1]),int(coordList[2]),int(coordList[2]),int(coordList[1])],clow=clow,chigh=chigh,if_colorbar=if_colorbar,cmap=cmap
                        )
    plt.locator_params(nbins=4)
    return fig