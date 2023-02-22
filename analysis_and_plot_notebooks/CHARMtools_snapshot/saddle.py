# code for compartmnet strength analysis, mainly adapted from 
import cooler 
import cooltools
import warnings
from cytoolz import merge
import numpy as np
import time
import logging
log = logging.getLogger(__name__)


def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

#     warnings.warn(
#         "Generating a saddleplot will be deprecated in future versions, "
#         + "please see https://github.com/open2c_examples for examples on how to plot saddles.",
#         DeprecationWarning,
#     )

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    # Old version
    # hist = np.bincount(x, minlength=len(binedges) + 1)

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout
    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure
    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins
    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist
    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist
    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    # Colorbar
    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["cbar"] = cb = plt.colorbar(img, **cbar_kws)
        # cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
        # # do linspace between vmin and vmax of 5 segments and trunc to 1 decimal:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        grid["cbar"] = cb = plt.colorbar(img, format=MinOneMaxFormatter(), **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings
    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    plt.grid(False)
    plt.axis("off")
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid

def saddle_strength(S, C):
    """
    modified by skelviper
    also returen A and B score, which is A/B region/inter
    Parameters
    ----------
    S, C : 2D arrays, square, same shape
        Saddle sums and counts, respectively

    Returns
    -------
    1D array
    Ratios of cumulative corner interaction scores, where the saddle data is
    grouped over the AA+BB corners and AB+BA corners with increasing extent.

    """
    m, n = S.shape
    if m != n:
        raise ValueError("`saddledata` should be square.")

    ratios = np.zeros(n)
    A = np.zeros(n)
    B = np.zeros(n)
    for k in range(1, n):
        intra_sum = np.nansum(S[0:k, 0:k]) + np.nansum(S[n - k : n, n - k : n])
        intra_count = np.nansum(C[0:k, 0:k]) + np.nansum(C[n - k : n, n - k : n])
        intra = intra_sum / intra_count
        A_sum = np.nansum(S[n - k : n, n - k : n])
        A_count = np.nansum(C[n - k : n, n - k : n])
        A_strength = A_sum / A_count
        B_sum = np.nansum(S[0:k, 0:k]) 
        B_count = np.nansum(C[0:k, 0:k])
        B_strength = B_sum / B_count

        inter_sum = np.nansum(S[0:k, n - k : n]) + np.nansum(S[n - k : n, 0:k])
        inter_count = np.nansum(C[0:k, n - k : n]) + np.nansum(C[n - k : n, 0:k])
        inter = inter_sum / inter_count

        ratios[k] = intra / inter
        A[k] = A_strength / inter
        B[k] = B_strength / inter
    return np.array([ratios,A,B]).T

def iterativeCorrection(matrix, v=None, M=50, tolerance=1e-5, verbose=False):
    """
    adapted from cytonised version in mirnylab
    original code from: ultracorrectSymmetricWithVector
    https://bitbucket.org/mirnylab/mirnylib/src/924bfdf5ed344df32743f4c03157b0ce49c675e6/mirnylib/numutils_new.pyx?at=default
    Main method for correcting DS and SS read data.
    Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction
    :param matrix: a scipy sparse matrix
    :param tolerance: Tolerance is the maximum allowed relative
                      deviation of the marginals.
    """
    if verbose:
        log.setLevel(logging.INFO)

    total_bias = np.ones(matrix.shape[0], 'float64')

    if np.isnan(matrix.sum()):
        log.warn("[iterative correction] the matrix contains nans, they will be replaced by zeros.")
        matrix.data[np.isnan(matrix.data)] = 0

    matrix = matrix.astype(float)
    W = matrix.tocoo()

    if np.abs(matrix - matrix.T).mean() / (1. * np.abs(matrix.mean())) > 1e-10:
        raise ValueError("Please provide symmetric matrix!")

    start_time = time.time()
    log.info("starting iterative correction")
    for iternum in range(M):
        iternum += 1
        s = np.array(W.sum(axis=1)).flatten()
        mask = (s == 0)
        s = s / np.mean(s[~mask])

        total_bias *= s
        deviation = np.abs(s - 1).max()

        s = 1.0 / s

        # The following code  is an optimization of this
        # for i in range(N):
        #     for j in range(N):
        #         W[i,j] = W[i,j] / (s[i] * s[j])

        W.data *= np.take(s, W.row)
        W.data *= np.take(s, W.col)
        if np.any(W.data > 1e100):
            log.error("*Error* matrix correction is producing extremely large values. "
                      "This is often caused by bins of low counts. Use a more stringent "
                      "filtering of bins.")
            exit(1)
        if verbose:
            if iternum % 5 == 0:
                end_time = time.time()
                estimated = (float(M - iternum) * (end_time - start_time)) / iternum
                m, sec = divmod(estimated, 60)
                h, m = divmod(m, 60)
                log.info("pass {} Estimated time {:.0f}:{:.0f}:{:.0f}".format(iternum, h, m, sec))
                log.info("max delta - 1 = {} ".format(deviation))

        if deviation < tolerance:
            log.info("[iterative correction] {} iterations used\n".format(iternum + 1))
            break

    # scale the total bias such that the sum is 1.0
    corr = total_bias[total_bias != 0].mean()
    total_bias /= corr
    W.data = W.data * corr * corr
    if np.any(W.data > 1e10):
        log.error("*Error* matrix correction produced extremely large values. "
                  "This is often caused by bins of low counts. Use a more stringent "
                  "filtering of bins.")
        exit(1)

    return W.tocsr(), total_bias


def calc_compartment_strength(coolerPath,view_df,compartment,Q_LO=0,Q_HI=1,N_GROUPS=5,I_th=1):
    try:
        clr = cooler.Cooler(coolerPath)
        expected = cooltools.expected_cis(clr,view_df)
        interaction_sum, interaction_count = cooltools.api.saddle.saddle(clr,expected,compartment,'cis', n_bins=N_GROUPS,qrange=(Q_LO,Q_HI),view_df=view_df)
        strength = saddle_strength(interaction_sum[1:-1,1:-1], interaction_count[1:-1,1:-1])[I_th]
        #strength = saddle_strength(iterativeCorrection(csr_matrix(interaction_sum[1:-1,1:-1]))[0].toarray(),
        #                           iterativeCorrection(csr_matrix(interaction_count[1:-1,1:-1]))[0].toarray())[I_th]
    except:
        strength = np.array([np.nan,np.nan,np.nan])
    return np.append(coolerPath,strength)