# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 19, 2018

"""
Plotting functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import os
import socket
from copy import deepcopy

import numpy as np
import numpy.ma as ma
from scipy.interpolate import splev, splprep
from scipy.ndimage.filters import gaussian_filter

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
mpl.use('Agg')

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.patches import Arrow

tRed = list(np.array([226,101,95]) / 255.)
tBlue = list(np.array([96,149,201]) / 255.)
tGreen = list(np.array([170,196,109]) / 255.)

import getdist
from getdist import plots, mcsamples

import ternary
from ternary.heatmapping import polygon_generator

import shapely.geometry as geometry

from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Delaunay

from utils.enums import DataType, str_enum
from utils.enums import Likelihood, ParamTag, StatCateg, Texture
from utils.misc import get_units, make_dir, solve_ratio, interval
from utils.fr import angles_to_u, angles_to_fr, SCALE_BOUNDARIES


BAYES_K = 1.   # Strong degree of belief.
# BAYES_K = 3/2. # Very strong degree of belief.
# BAYES_K = 2.   # Decisive degree of belief


LV_ATMO_90PC_LIMITS = {
    3: (2E-24, 1E-1),
    4: (2.7E-28, 3.16E-25),
    5: (1.5E-32, 1.12E-27),
    6: (9.1E-37, 2.82E-30),
    7: (3.6E-41, 1.77E-32),
    8: (1.4E-45, 1.00E-34)
}


PS = 8.203E-20 # GeV^{-1}
PLANCK_SCALE = {
    5: PS,
    6: PS**2,
    7: PS**3,
    8: PS**4
}


if os.path.isfile('./plot_llh/paper.mplstyle'):
    plt.style.use('./plot_llh/paper.mplstyle')
elif os.environ.get('GOLEMSOURCEPATH') is not None:
    plt.style.use(os.environ['GOLEMSOURCEPATH']+'/GolemFit/scripts/paper/paper.mplstyle')
if 'submitter' in socket.gethostname():
    rc('text', usetex=False)

mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{xcolor}',
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']
mpl.rcParams['text.latex.unicode'] = True


def gen_figtext(args):
    """Generate the figure text."""
    t = r'$'
    if args.data is DataType.REAL:
        t += r'\textbf{IceCube\:Preliminary}' + '$\n$'
    elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
        t += r'{\rm\bf IceCube\:Simulation}' + '$\n$'
        t += r'\rm{Injected\:composition}'+r'\:=\:({0})_\oplus'.format(
            solve_ratio(args.injected_ratio).replace('_', ':')
        ) + '$\n$'
    t += r'{\rm Source\:composition}'+r'\:=\:({0})'.format(
        solve_ratio(args.source_ratio).replace('_', ':')
    ) + r'_\text{S}'
    t += '$\n$' + r'{\rm Dimension}'+r' = {0}$'.format(args.dimension)
    return t


def texture_label(x, dim):
    cpt = r'c' if dim % 2 == 0 else r'a'
    if x == Texture.OEU:
        # return r'$\mathcal{O}_{e\mu}$'
        return r'$\mathring{'+cpt+r'}_{e\mu}^{('+str(int(dim))+r')}$'
    elif x == Texture.OET:
        # return r'$\mathcal{O}_{e\tau}$'
        return r'$\mathring{'+cpt+r'}_{\tau e}^{('+str(int(dim))+r')}$'
    elif x == Texture.OUT:
        # return r'$\mathcal{O}_{\mu\tau}$'
        return r'$\mathring{'+cpt+r'}_{\mu\tau}^{('+str(int(dim))+r')}$'
    else:
        raise AssertionError


def cmap_discretize(cmap, N):
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def get_limit(scales, statistic, args, mask_initial=False, return_interp=False):
    max_st = np.max(statistic)
    print 'scales, stat', zip(scales, statistic)
    if args.stat_method is StatCateg.BAYESIAN:
        if (statistic[0] - max_st) > np.log(10**(BAYES_K)):
            raise AssertionError('Discovered LV!')

    try:
        tck, u = splprep([scales, statistic], s=0)
    except:
        print 'Failed to spline'
        # return None
        raise
    sc, st = splev(np.linspace(0, 1, 1000), tck)

    if mask_initial:
        scales_rm = sc[sc >= scales[1]]
        statistic_rm = st[sc >= scales[1]]
    else:
        scales_rm = sc
        statistic_rm = st

    min_idx = np.argmin(scales)
    null = statistic[min_idx]
    # if np.abs(statistic_rm[0] - null) > 0.8:
    #     print 'Warning, null incompatible with smallest scanned scale! For ' \
    #         'DIM {0} [{1}, {2}, {3}]!'.format(
    #             args.dimension, *args.source_ratio
    #         )
    #     null = statistic_rm[0]
    if args.stat_method is StatCateg.BAYESIAN:
        reduced_ev = -(statistic_rm - null)
        print '[reduced_ev > np.log(10**(BAYES_K))]', np.sum([reduced_ev > np.log(10**(BAYES_K))])
        al = scales_rm[reduced_ev > np.log(10**(BAYES_K))]
    else:
        assert 0
    if len(al) == 0:
        print 'No points for DIM {0} [{1}, {2}, {3}]!'.format(
            args.dimension, *args.source_ratio
        )
        return None
    re = -(statistic-null)[scales > al[0]]
    if np.sum(re < np.log(10**(BAYES_K)) - 0.1) >= 2:
        print 'Warning, peaked contour does not exclude large scales! For ' \
            'DIM {0} [{1}, {2}, {3}]!'.format(
                args.dimension, *args.source_ratio
            )
        return None
    if np.sum(re >= np.log(10**(BAYES_K)) + 0.0) < 2:
        print 'Warning, only single point above threshold! For ' \
            'DIM {0} [{1}, {2}, {3}]!'.format(
                args.dimension, *args.source_ratio
            )
        return None

    if return_interp:
        return (scales_rm, reduced_ev)

    # Divide by 2 to convert to standard SME coefficient
    lim = al[0] - np.log10(2.)
    # lim = al[0]
    print 'limit = {0}'.format(lim)
    return lim


def heatmap(data, scale, vmin=None, vmax=None, style='triangular'):
    for k, v in data.items():
        data[k] = np.array(v)
    style = style.lower()[0]
    if style not in ["t", "h", 'd']:
        raise ValueError("Heatmap style must be 'triangular', 'dual-triangular', or 'hexagonal'")

    vertices_values = polygon_generator(data, scale, style)

    vertices = []
    for v, value in vertices_values:
        vertices.append(v)
    return vertices


def get_tax(ax, scale, ax_labels, rot_ax_labels=False, fontsize=23):
    ax.set_aspect('equal')

    # Boundary and Gridlines
    fig, tax = ternary.figure(ax=ax, scale=scale)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color='grey', multiple=scale/5., linewidth=0.5, alpha=0.7, ls='--')
    # tax.gridlines(color='grey', multiple=scale/10., linewidth=0.2, alpha=1, ls=':')

    # Set Axis labels and Title
    if rot_ax_labels: roty, rotz = (-60, 60)
    else: roty, rotz = (0, 0)
    tax.bottom_axis_label(
        ax_labels[0], fontsize=fontsize+4,
        position=(0.55, -0.10/2, 0.5), rotation=0
    )
    tax.right_axis_label(
        ax_labels[1], fontsize=fontsize+4,
        position=(2./5+0.1, 3./5+0.06, 0), rotation=roty
    )
    tax.left_axis_label(
        ax_labels[2], fontsize=fontsize+4,
        position=(-0.15, 3./5+0.1, 2./5), rotation=rotz
    )

    # Remove default Matplotlib axis
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()

    # Set ticks
    ticks = np.linspace(0, 1, 6)
    tax.ticks(ticks=ticks, locations=ticks*scale, axis='lr', linewidth=1,
              offset=0.03, fontsize=fontsize, tick_formats='%.1f')
    tax.ticks(ticks=ticks, locations=ticks*scale, axis='b', linewidth=1,
              offset=0.02, fontsize=fontsize, tick_formats='%.1f')
    # tax.ticks()

    tax._redraw_labels()

    return tax


def project(p):
    """Convert from flavour to cartesian."""
    a, b, c = p
    x = a + b/2.
    y = b * np.sqrt(3)/2.
    return [x, y]


def project_toflavour(p, nbins):
    """Convert from cartesian to flavour space."""
    x, y = p
    b = y / (np.sqrt(3)/2.)
    a = x - b/2.
    return [a, b, nbins-a-b]


def tax_fill(ax, points, nbins, **kwargs):
    pol = np.array(map(project, points))
    ax.fill(pol.T[0]*nbins, pol.T[1]*nbins, **kwargs)


def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull
    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])
    coords = np.array([point.coords[0]
                       for point in points])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = np.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = np.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = np.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = np.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points


def flavour_contour(frs, nbins, coverage, ax=None, smoothing=0.4,
                    hist_smooth=0.05, plot=True, fill=False, oversample=1.,
                    delaunay=False, d_alpha=1.5, d_gauss=0.08, debug=False,
                    **kwargs):
    """Plot the flavour contour for a specified coverage."""
    # Histogram in flavour space
    os_nbins = nbins * oversample
    H, b = np.histogramdd(
        (frs[:,0], frs[:,1], frs[:,2]),
        bins=(os_nbins+1, os_nbins+1, os_nbins+1),
        range=((0, 1), (0, 1), (0, 1))
    )
    H = H / np.sum(H)

    # 3D smoothing
    H_s = gaussian_filter(H, sigma=hist_smooth)

    # Finding coverage
    H_r = np.ravel(H_s)
    H_rs = np.argsort(H_r)[::-1]
    H_crs = np.cumsum(H_r[H_rs])
    thres = np.searchsorted(H_crs, coverage/100.)
    mask_r = np.zeros(H_r.shape)
    mask_r[H_rs[:thres]] = 1
    mask = mask_r.reshape(H_s.shape)

    # Get vertices inside covered region
    binx = np.linspace(0, 1, os_nbins+1)
    interp_dict = {}
    for i in xrange(len(binx)):
        for j in xrange(len(binx)):
            for k in xrange(len(binx)):
                if mask[i][j][k] == 1:
                    interp_dict[(i, j, k)] = H_s[i, j, k]
    vertices = np.array(heatmap(interp_dict, os_nbins))
    points = vertices.reshape((len(vertices)*3, 2))
    if debug:
        ax.scatter(*(points/float(oversample)).T, marker='o', s=3, alpha=1.0, color=kwargs['color'], zorder=9)

    pc = geometry.MultiPoint(points)
    if not delaunay:
        # Convex full to find points forming exterior bound
        polygon = pc.convex_hull
        ex_cor = np.array(list(polygon.exterior.coords))
    else:
        # Delaunay
        concave_hull, edge_points = alpha_shape(pc, alpha=d_alpha)
        polygon = geometry.Polygon(concave_hull.buffer(1))
        if d_gauss == 0.:
            ex_cor = np.array(list(polygon.exterior.coords))
        else:
            ex_cor = gaussian_filter(
                np.array(list(polygon.exterior.coords)), sigma=d_gauss
            )

    # Join points with a spline
    tck, u = splprep([ex_cor.T[0], ex_cor.T[1]], s=0., per=1, k=1)
    xi, yi = map(np.array, splev(np.linspace(0, 1, 300), tck))

    # Spline again to smooth
    if smoothing != 0:
        tck, u = splprep([xi, yi], s=smoothing, per=1, k=3)
        xi, yi = map(np.array, splev(np.linspace(0, 1, 600), tck))

    xi /= float(oversample)
    yi /= float(oversample)
    ev_polygon = np.dstack((xi, yi))[0]

    # Remove points interpolated outside flavour triangle
    f_ev_polygon = np.array(map(lambda x: project_toflavour(x, nbins), ev_polygon))

    xf, yf, zf = f_ev_polygon.T
    mask = np.array((xf < 0) | (yf < 0) | (zf < 0) | (xf > nbins) |
                    (yf > nbins) | (zf > nbins))
    ev_polygon = np.dstack((xi[~mask], yi[~mask]))[0]

    # Plot
    if plot:
        if fill:
            ax.fill(
                ev_polygon.T[0], ev_polygon.T[1],
                label=r'{0}\%'.format(int(coverage)), **kwargs
            )
        else:
            ax.plot(
                ev_polygon.T[0], ev_polygon.T[1],
                label=r'{0}\%'.format(int(coverage)), **kwargs
            )
    else:
        return ev_polygon

def plot_Tchain(Tchain, axes_labels, ranges):
    """Plot the Tchain using getdist."""
    Tsample = mcsamples.MCSamples(
        samples=Tchain, labels=axes_labels, ranges=ranges
    )

    Tsample.updateSettings({'contours': [0.90, 0.99]})
    Tsample.num_bins_2D=10
    Tsample.fine_bins_2D=50
    Tsample.smooth_scale_2D=0.05

    g = plots.getSubplotPlotter()
    g.settings.num_plot_contours = 2
    g.settings.axes_fontsize = 10
    g.settings.figure_legend_frame = False
    g.settings.lab_fontsize = 20
    g.triangle_plot(
        [Tsample], filled=True,# contour_colors=['green', 'lightgreen']
    )
    return g


def chainer_plot(infile, outfile, outformat, args, llh_paramset, fig_text=None,
                 labels=None, ranges=None):
    """Make the triangle plot."""
    if hasattr(args, 'plot_elements'):
        if not args.plot_angles and not args.plot_elements:
            return
    elif not args.plot_angles:
        return

    if not isinstance(infile, np.ndarray):
        raw = np.load(infile)
    else:
        raw = infile
    print 'raw.shape', raw.shape
    print 'raw', raw

    make_dir(outfile), make_dir
    if fig_text is None:
        fig_text = gen_figtext(args)

    if labels is None: axes_labels = llh_paramset.labels
    else: axes_labels = labels
    if ranges is None: ranges = llh_paramset.ranges

    if args.plot_angles:
        print "Making triangle plots"
        Tchain = raw
        g = plot_Tchain(Tchain, axes_labels, ranges)

        mpl.pyplot.figtext(0.6, 0.7, fig_text, fontsize=20)

        # for i_ax_1, ax_1 in enumerate(g.subplots):
        #     for i_ax_2, ax_2 in enumerate(ax_1):
        #         if i_ax_1 == i_ax_2:
        #             itv = interval(Tchain[:,i_ax_1], percentile=90.)
        #             for l in itv:
        #                 ax_2.axvline(l, color='gray', ls='--')
        #                 ax_2.set_title(r'${0:.2f}_{{{1:.2f}}}^{{+{2:.2f}}}$'.format(
        #                     itv[1], itv[0]-itv[1], itv[2]-itv[1]
        #                 ), fontsize=10)

        # if not args.fix_mixing:
        #     sc_index = llh_paramset.from_tag(ParamTag.SCALE, index=True)
        #     itv = interval(Tchain[:,sc_index], percentile=90.)
        #     mpl.pyplot.figtext(
        #         0.5, 0.3, 'Scale 90% Interval = [1E{0}, 1E{1}], Center = '
        #         '1E{2}'.format(itv[0], itv[2], itv[1])
        #     )

        for of in outformat:
            print 'Saving', outfile+'_angles.'+of
            g.export(outfile+'_angles.'+of)

    if not hasattr(args, 'plot_elements'):
        return

    if args.plot_elements:
        print "Making triangle plots"
        if args.fix_mixing_almost:
            raise NotImplementedError
        nu_index = llh_paramset.from_tag(ParamTag.NUISANCE, index=True)
        fr_index = llh_paramset.from_tag(ParamTag.MMANGLES, index=True)
        sc_index = llh_paramset.from_tag(ParamTag.SCALE, index=True)
        if not args.fix_source_ratio:
            sr_index = llh_paramset.from_tag(ParamTag.SRCANGLES, index=True)

        nu_elements = raw[:,nu_index]
        fr_elements = np.array(map(flat_angles_to_u, raw[:,fr_index]))
        sc_elements = raw[:,sc_index]
        if not args.fix_source_ratio:
            sr_elements = np.array(map(angles_to_fr, raw[:,sr_index]))
        if args.fix_source_ratio:
            Tchain = np.column_stack(
                [nu_elements, fr_elements, sc_elements]
            )
        else:
            Tchain = np.column_stack(
                [nu_elements, fr_elements, sc_elements, sr_elements]
            )

        trns_ranges = np.array(ranges)[nu_index,].tolist()
        trns_axes_labels = np.array(axes_labels)[nu_index,].tolist()
        if args.fix_mixing is not MixingScenario.NONE:
            trns_axes_labels += \
                [r'\mid \tilde{U}_{e1} \mid'    , r'\mid \tilde{U}_{e2} \mid'    , r'\mid \tilde{U}_{e3} \mid'     , \
                 r'\mid \tilde{U}_{\mu1} \mid'  , r'\mid \tilde{U}_{\mu2} \mid'  , r'\mid \tilde{U}_{\mu3} \mid'   , \
                 r'\mid \tilde{U}_{\tau1} \mid' , r'\mid \tilde{U}_{\tau2} \mid' , r'\mid \tilde{U}_{\tau3} \mid']
            trns_ranges += [(0, 1)] * 9
        if not args.fix_scale:
            trns_axes_labels += [np.array(axes_labels)[sc_index].tolist()]
            trns_ranges += [np.array(ranges)[sc_index].tolist()]
        if not args.fix_source_ratio:
            trns_axes_labels += [r'\phi_e', r'\phi_\mu', r'\phi_\tau']
            trns_ranges += [(0, 1)] * 3

        g = plot_Tchain(Tchain, trns_axes_labels, trns_ranges)

        if args.data is DataType.REAL:
            plt.text(0.8, 0.7, 'IceCube Preliminary', color='red', fontsize=15,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            plt.text(0.8, 0.7, 'IceCube Simulation', color='red', fontsize=15,
                     ha='center', va='center')

        mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            print 'Saving', outfile+'_elements'+of
            g.export(outfile+'_elements.'+of)


def plot_statistic(data, outfile, outformat, args, scale_param, label=None):
    """Make MultiNest factor or LLH value plot."""
    print 'Making Statistic plot'
    fig_text = gen_figtext(args)
    if label is not None: fig_text += '\n' + label

    print 'data', data
    print 'data.shape', data.shape
    print 'outfile', outfile
    try:
        scales, statistic = ma.compress_rows(data).T
        lim = get_limit(deepcopy(scales), deepcopy(statistic), args, mask_initial=True)
        tck, u = splprep([scales, statistic], s=0)
    except:
        return
    sc, st = splev(np.linspace(0, 1, 1000), tck)
    scales_rm = sc[sc >= scales[1]]
    statistic_rm = st[sc >= scales[1]]

    min_idx = np.argmin(scales)
    null = statistic[min_idx]
    # fig_text += '\nnull lnZ = {0:.2f}'.format(null)

    if args.stat_method is StatCateg.BAYESIAN:
        reduced_ev = -(statistic_rm - null)
    elif args.stat_method is StatCateg.FREQUENTIST:
        reduced_ev = -2*(statistic_rm - null)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    xlims = SCALE_BOUNDARIES[args.dimension]
    ax.set_xlim(xlims)
    ax.set_xlabel(r'${\rm log}_{10}\left[\Lambda^{-1}_{'+ \
                  r'{0}'.format(args.dimension)+r'}'+ \
                  get_units(args.dimension)+r'\right]$', fontsize=16)

    if args.stat_method is StatCateg.BAYESIAN:
        ax.set_ylabel(r'$\text{Bayes\:Factor}\:\left[\text{ln}\left(B_{0/1}\right)\right]$')
    elif args.stat_method is StatCateg.FREQUENTIST:
        ax.set_ylabel(r'$-2\Delta {\rm LLH}$')

    # ymin = np.round(np.min(reduced_ev) - 1.5)
    # ymax = np.round(np.max(reduced_ev) + 1.5)
    # ax.set_ylim((ymin, ymax))

    ax.scatter(scales[1:], -(statistic[1:]-null), color='r')
    ax.plot(scales_rm, reduced_ev, color='k', linewidth=1, alpha=1, ls='-')

    ax.axhline(y=np.log(10**(BAYES_K)), color='red', alpha=1., linewidth=1.2, ls='--')
    ax.axvline(x=lim, color='red', alpha=1., linewidth=1.2, ls='--')

    at = AnchoredText(
        fig_text, prop=dict(size=10), frameon=True, loc=4
    )
    at.patch.set_boxstyle("round,pad=0.1,rounding_size=0.5")
    ax.add_artist(at)
    make_dir(outfile)
    for of in outformat:
        print 'Saving as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def plot_table_sens(data, outfile, outformat, args, show_lvatmo=True):
    print 'Making TABLE sensitivity plot'
    argsc = deepcopy(args)

    dims = args.dimensions
    srcs = args.source_ratios
    if args.texture is Texture.NONE:
        textures = [Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]

    if len(srcs) > 3:
        raise NotImplementedError

    xlims = (-60, -20)
    ylims = (0.5, 1.5)

    colour = {2:'red', 1:'blue', 0:'green'}
    rgb_co = {2:[1,0,0], 1:[0,0,1], 0:[0.0, 0.5019607843137255, 0.0]}

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(len(dims), 1)
    gs.update(hspace=0.15)

    first_ax = None
    legend_log = []
    legend_elements = []

    for idim, dim in enumerate(dims):
        print '|||| DIM = {0}'.format(dim)
        argsc.dimension = dim
        gs0 = gridspec.GridSpecFromSubplotSpec(
            len(textures), 1, subplot_spec=gs[idim], hspace=0
        )

        for itex, tex in enumerate(textures):
            argsc.texture = tex
            ylabel = texture_label(tex, dim)
            # if angles == 2 and ian == 0: continue
            ax = fig.add_subplot(gs0[itex])

            if first_ax is None:
                first_ax = ax

            ax.set_xlim(xlims)
            ax.set_ylim(ylims)
            ax.set_yticks([], minor=True)
            ax.set_yticks([1.], minor=False)
            ax.set_yticklabels([ylabel], fontsize=13)
            ax.yaxis.tick_right()
            for xmaj in ax.xaxis.get_majorticklocs():
                ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.2, linewidth=1)
            ax.get_xaxis().set_visible(False)
            if itex == len(textures) - 2:
                ax.spines['bottom'].set_alpha(0.6)
            elif itex == len(textures) - 1:
                ax.text(
                    -0.04, 1.0, '$d = {0}$'.format(dim), fontsize=16,
                    rotation='90', verticalalignment='center',
                    transform=ax.transAxes
                )
                # dim_label_flag = False
                ax.spines['top'].set_alpha(0.6)
                ax.spines['bottom'].set_alpha(0.6)

            for isrc, src in enumerate(srcs):
                print '== src', src
                argsc.source_ratio = src

                if dim in PLANCK_SCALE.iterkeys():
                    ps = np.log10(PLANCK_SCALE[dim])
                    if ps < xlims[0]:
                        ax.annotate(
                            s='', xy=(xlims[0], 1), xytext=(xlims[0]+1, 1),
                            arrowprops={'arrowstyle': '->, head_length=0.2',
                                        'lw': 1, 'color':'purple'}
                        )
                    elif ps > xlims[1]:
                        ax.annotate(
                            s='', xy=(xlims[1]-1, 1), xytext=(xlims[1], 1),
                            arrowprops={'arrowstyle': '<-, head_length=0.2',
                                        'lw': 1, 'color':'purple'}
                        )
                    else:
                        ax.axvline(x=ps, color='purple', alpha=1., linewidth=1.5)

                try:
                    scales, statistic = ma.compress_rows(data[idim][isrc][itex]).T
                except: continue
                lim = get_limit(deepcopy(scales), deepcopy(statistic), argsc, mask_initial=True)
                if lim is None: continue

                ax.axvline(x=lim, color=colour[isrc], alpha=1., linewidth=1.5)
                ax.add_patch(patches.Rectangle(
                    (lim, ylims[0]), 100, np.diff(ylims), fill=True,
                    facecolor=colour[isrc], alpha=0.3, linewidth=0
                ))

                if isrc not in legend_log:
                    legend_log.append(isrc)
                    label = r'$\left('+r'{0}'.format(solve_ratio(src)).replace('_',':')+ \
                            r'\right)_{\text{S}}\:\text{at\:source}$'
                    legend_elements.append(
                        Patch(facecolor=rgb_co[isrc]+[0.3],
                              edgecolor=rgb_co[isrc]+[1], label=label)
                    )

            if itex == len(textures)-1 and show_lvatmo:
                LV_lim = np.log10(LV_ATMO_90PC_LIMITS[dim])
                ax.add_patch(patches.Rectangle(
                    (LV_lim[1], ylims[0]), LV_lim[0]-LV_lim[1], np.diff(ylims),
                    fill=False, hatch='\\\\'
                ))

    ax.get_xaxis().set_visible(True)
    ax.set_xlabel(r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} (\Lambda^{-1}_{d}\:/\:{\rm GeV}^{-d+4})\: ]$',
                 labelpad=5, fontsize=19)
    ax.tick_params(axis='x', labelsize=16)

    purple = [0.5019607843137255, 0.0, 0.5019607843137255]
    if show_lvatmo:
        legend_elements.append(
            Patch(facecolor='none', hatch='\\\\', edgecolor='k', label='IceCube, Nature.Phy.14,961(2018)')
        )
    legend_elements.append(
        Patch(facecolor=purple+[0.7], edgecolor=purple+[1], label='Planck Scale Expectation')
    )
    legend = first_ax.legend(
        handles=legend_elements, prop=dict(size=11), loc='upper left',
        title='Excluded regions', framealpha=1., edgecolor='black',
        frameon=True
    )
    first_ax.set_zorder(10)
    plt.setp(legend.get_title(), fontsize='11')
    legend.get_frame().set_linestyle('-')

    ybound = 0.595
    if args.data is DataType.REAL:
        # fig.text(0.295, 0.684, 'IceCube Preliminary', color='red', fontsize=13,
        fig.text(0.278, ybound, r'\bf IceCube Preliminary', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)
    elif args.data is DataType.REALISATION:
        fig.text(0.278, ybound-0.05, r'\bf IceCube Simulation', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)
    else:
        fig.text(0.278, ybound, r'\bf IceCube Simulation', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)

    make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def plot_x(data, outfile, outformat, args, normalise=False):
    """Limit plot as a function of the source flavour ratio for each operator
    texture."""
    print 'Making X sensitivity plot'
    dim = args.dimension
    if dim < 5: normalise = False
    srcs = args.source_ratios
    x_arr = np.array([i[0] for i in srcs])
    if args.texture is Texture.NONE:
        textures = [Texture.OEU, Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]

    # Rearrange data structure
    r_data = np.full((
        data.shape[1], data.shape[0], data.shape[2], data.shape[3]
    ), np.nan)

    for isrc in xrange(data.shape[0]):
        for itex in xrange(data.shape[1]):
            r_data[itex][isrc] = data[isrc][itex]
    r_data = ma.masked_invalid(r_data)
    print r_data.shape, 'r_data.shape'

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)

    ylims = SCALE_BOUNDARIES[dim]
    if normalise:
        if dim == 5: ylims = (-24, -8)
        if dim == 6: ylims = (-12, 8)
        if dim == 7: ylims = (0, 20)
        if dim == 8: ylims = (12, 36)
    else:
        if dim == 3: ylims = (-28, -22)
        if dim == 4: ylims = (-35, -26)
        if dim == 5: SCALE_BOUNDARIES[5]
    xlims = (0, 1)

    colour = {0:'red', 2:'blue', 1:'green'}
    rgb_co = {0:[1,0,0], 2:[0,0,1], 1:[0.0, 0.5019607843137255, 0.0]}

    legend_log = []
    legend_elements = []
    labelsize = 13
    largesize = 17

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    xticks = [0, 1/3., 0.5, 2/3., 1]
    # xlabels = [r'$0$', r'$\frac{1}{3}$', r'$\frac{1}{2}$', r'$\frac{2}{3}$', r'$1$']
    xlabels = [r'$0$', r'$1 / 3$', r'$1/2$', r'$2/3$', r'$1$']
    ax.set_xticks([], minor=True)
    ax.set_xticks(xticks, minor=False)
    ax.set_xticklabels(xlabels, fontsize=largesize)
    if dim != 4 or dim != 3:
        yticks = range(ylims[0], ylims[1], 2) + [ylims[1]]
        ax.set_yticks(yticks, minor=False)
    if dim == 3 or dim == 4:
        yticks = range(ylims[0], ylims[1], 1) + [ylims[1]]
        ax.set_yticks(yticks, minor=False)
    # for ymaj in ax.yaxis.get_majorticklocs():
    #     ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.2, linewidth=1)
    for xmaj in xticks:
        if xmaj == 1/3.:
            ax.axvline(x=xmaj, ls='--', color='gray', alpha=0.5, linewidth=0.3)
        # else:
        #     ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.2, linewidth=1)

    ax.text(
        (1/3.)+0.01, 0.01,  r'$(0.33:0.66:0)_\text{S}$', fontsize=labelsize,
        transform=ax.transAxes, rotation='vertical', va='bottom'
    )
    ax.text(
        0.96, 0.01, r'$(1:0:0)_\text{S}$', fontsize=labelsize,
        transform=ax.transAxes, rotation='vertical', va='bottom', ha='left'
    )
    ax.text(
        0.01, 0.01, r'$(0:1:0)_\text{S}$', fontsize=labelsize,
        transform=ax.transAxes, rotation='vertical', va='bottom'
    )
    yl = 0.55
    if dim == 3: yl = 0.65
    ax.text(
        0.03, yl, r'${\rm \bf Excluded}$', fontsize=largesize,
        transform=ax.transAxes, color = 'g', rotation='vertical', zorder=10
    )
    ax.text(
        0.95, 0.55, r'${\rm \bf Excluded}$', fontsize=largesize,
        transform=ax.transAxes, color = 'b', rotation='vertical', zorder=10
    )

    for itex, tex in enumerate(textures):
        print '|||| TEX = {0}'.format(tex)
        lims = np.full(len(srcs), np.nan)

        for isrc, src in enumerate(srcs):
            x = src[0]
            print '|||| X = {0}'.format(x)
            args.source_ratio = src
            d = r_data[itex][isrc]
            if np.sum(d.mask) > 2: continue
            scales, statistic = ma.compress_rows(d).T
            lim = get_limit(deepcopy(scales), deepcopy(statistic), args, mask_initial=True)
            if lim is None: continue
            if normalise:
                lim -= np.log10(PLANCK_SCALE[dim])
            lims[isrc] = lim

        lims = ma.masked_invalid(lims)
        size = np.sum(~lims.mask)
        if size == 0: continue

        print 'x_arr, lims', zip(x_arr, lims)
        if normalise:
            zeropoint = 100
        else:
            zeropoint = 0
        lims[lims.mask] = zeropoint

        l0 = np.argwhere(lims == zeropoint)[0]
        h0 = len(lims) - np.argwhere(np.flip(lims, 0) == zeropoint)[0]
        lims[int(l0):int(h0)] = zeropoint

        x_arr_a = [x_arr[0]-0.1] + list(x_arr)
        x_arr_a = list(x_arr_a) + [x_arr_a[-1]+0.1]
        lims = [lims[0]] + list(lims)
        lims = list(lims) + [lims[-1]]

        s = 0.2
        g = 2
        if dim == 3 and tex == Texture.OUT:
            s = 0.4
            g = 4
        if dim in (4,5) and tex == Texture.OUT:
            s = 0.5
            g = 5
        if dim == 7 and tex == Texture.OET:
            s = 1.6
            g = 2
        if dim == 7 and tex == Texture.OUT:
            s = 2.0
            g = 20
        if dim == 8 and tex == Texture.OET:
            s = 0.8
            g = 6
        if dim == 8 and tex == Texture.OUT:
            s = 1.7
            g = 8

        # ax.scatter(x_arr_a, lims, color='black', s=1)
        tck, u = splprep([x_arr_a, lims], s=0, k=1)
        x, y = splev(np.linspace(0, 1, 200), tck)
        tck, u = splprep([x, y], s=s)
        x, y = splev(np.linspace(0, 1, 400), tck)
        y = gaussian_filter(y, sigma=g)
        ax.fill_between(x, y, zeropoint, color=rgb_co[itex]+[0.3])
        # ax.scatter(x, y, color='black', s=1)
        # ax.scatter(x_arr_a, lims, color=rgb_co[itex], s=8)

        if itex not in legend_log:
            legend_log.append(itex)
            # label = texture_label(tex, dim)[:-1] + r'\:{\rm\:texture}$'
            label = texture_label(tex, dim)[:-1] + r'\:({\rm this\:work})$'
            legend_elements.append(
                Patch(facecolor=rgb_co[itex]+[0.3],
                      edgecolor=rgb_co[itex]+[1], label=label)
            )

    LV_lim = np.log10(LV_ATMO_90PC_LIMITS[dim])
    if normalise:
        LV_lim -= np.log10(PLANCK_SCALE[dim])
    ax.add_patch(patches.Rectangle(
        (xlims[0], LV_lim[1]), np.diff(xlims), LV_lim[0]-LV_lim[1],
        fill=False, hatch='\\\\'
    ))

    if dim in PLANCK_SCALE:
        ps = np.log10(PLANCK_SCALE[dim])
        if normalise and dim == 6:
            ps -= np.log10(PLANCK_SCALE[dim])
            ax.add_patch(Arrow(
                0.24, -0.009, 0, -5, width=0.12, capstyle='butt',
                facecolor='purple', fill=True, alpha=0.8,
                edgecolor='darkmagenta'
            ))
            ax.add_patch(Arrow(
                0.78, -0.009, 0, -5, width=0.12, capstyle='butt',
                facecolor='purple', fill=True, alpha=0.8,
                edgecolor='darkmagenta'
            ))

            ax.text(
                0.26, 0.5, r'${\rm \bf Quantum\:Gravity\:Frontier}$',
                fontsize=largesize-2, transform=ax.transAxes, va='top',
                ha='left', color='purple'
            )
        if dim > 5:
            ax.axhline(y=ps, color='purple', alpha=1., linewidth=1.5)

    cpt = r'c' if dim % 2 == 0 else r'a'
    if normalise:
        ft = r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} \left (\mathring{'+cpt+r'}^{(' + \
                r'{0}'.format(args.dimension)+r')}\cdot{\rm E}_{\:\rm P}'
        if dim > 5: ft += r'^{\:'+ r'{0}'.format(args.dimension-4)+ r'}'
        ft += r'\right )\: ]$'
        fig.text(
            0.01, 0.5, ft, ha='left',
            va='center', rotation='vertical', fontsize=largesize
        )
    else:
        fig.text(
            0.01, 0.5,
            r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} \left (\mathring{'+cpt+r'}^{(' +
            r'{0}'.format(args.dimension)+r')}\:' + get_units(args.dimension) +
            r'\right )\: ]$', ha='left',
            va='center', rotation='vertical', fontsize=largesize
        )

    ax.set_xlabel(
        r'${\rm Source\:Composition}\:[\:\left (\:x:1-x:0\:\right )_\text{S}\:]$',
        labelpad=10, fontsize=largesize
    )
    ax.tick_params(axis='x', labelsize=largesize-1)

    purple = [0.5019607843137255, 0.0, 0.5019607843137255]
    # legend_elements.append(
    #     Patch(facecolor=purple+[0.7], edgecolor=purple+[1], label='Planck Scale Expectation')
    # )
    legend_elements.append(
        Patch(facecolor='none', hatch='\\\\', edgecolor='k', label='IceCube [TODO]')
    )
    legend = ax.legend(
        handles=legend_elements, prop=dict(size=labelsize-2),
        loc='upper center', title='Excluded regions', framealpha=1.,
        edgecolor='black', frameon=True, bbox_to_anchor=(0.5, 1)
    )
    plt.setp(legend.get_title(), fontsize=labelsize)
    legend.get_frame().set_linestyle('-')

    # ybound = 0.14
    # if args.data is DataType.REAL:
    #     fig.text(0.7, ybound, r'\bf IceCube Preliminary', color='red', fontsize=13,
    #              ha='center', va='center', zorder=11)
    # elif args.data is DataType.REALISATION:
    #     fig.text(0.7, ybound-0.05, r'\bf IceCube Simulation', color='red', fontsize=13,
    #              ha='center', va='center', zorder=11)
    # else:
    #     fig.text(0.7, ybound, r'\bf IceCube Simulation', color='red', fontsize=13,
    #              ha='center', va='center', zorder=11)

    make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile + '.' + of)
        fig.savefig(outfile + '.' + of, bbox_inches='tight', dpi=150)
