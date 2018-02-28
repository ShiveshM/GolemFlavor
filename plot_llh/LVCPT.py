
# coding: utf-8

## The Theory

import numpy
import MinimalTools as MT
import PhysConst as PC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mco
import matplotlib as mpl
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy as sp
from numpy import linalg as LA

use_cython = False

if use_cython:
    import cython.cLVCPT as clv

mpl.rc('font', family='serif', size=20)

pc = PC.PhysicsConstants()

degree = np.pi/180.0
pc.th12 = 33.36*degree#33.48*degree
pc.th23 = 45.*degree#42.3*degree
pc.th13 = 8.66*degree#8.5*degree
pc.delta1 = 0.0#300.0*degree#306.*degree # perhaps better just set to 0.
pc.dm21sq = 7.5e-5
pc.dm31sq = 2.47e-3#2.457e-3
pc.Refresh()

MT.calcU(pc)
DELTAM2 = MT.flavorM2(pc)

def Hamiltonian(Enu, LVATERM = np.zeros((3,3), dtype=numpy.complex),
                LVCTERM = np.zeros((3,3), dtype=numpy.complex)):
    return DELTAM2/(2.0*Enu) + LVATERM + Enu*LVCTERM 

def OscProbFromMixingMatrix(alpha, beta, MixMatrix):
    return sum([(np.absolute(MixMatrix[i][alpha])*np.absolute(MixMatrix[i][beta]))**2 for i in range(pc.numneu)] )
    #return sum([(np.absolute(MixMatrix[i][alpha]))**2*(np.absolute(MixMatrix[i][beta]))**2 for i in range(pc.numneu)] )
    #prob = 0.0;
    #for i in range(pc.numneu) :
    #    prob += (np.absolute(MixMatrix[i][alpha]))**2*(np.absolute(MixMatrix[i][beta]))**2
    #return prob

def OscProb(alpha, Enu, LVATERM = np.zeros((3,3), dtype=numpy.complex),
                LVCTERM = np.zeros((3,3), dtype=numpy.complex)):
    eigvals, eigvec = MT.eigenvectors(Hamiltonian(Enu, LVATERM=LVATERM, LVCTERM=LVCTERM))
    #print eigvec.dtype
    if use_cython:
        return [ clv.OscProbFromMixingMatrix(alpha,beta,eigvec) for beta in range(pc.numneu)]
    else:
        return [ OscProbFromMixingMatrix(alpha,beta,eigvec) for beta in range(pc.numneu)]

def FlavorRatio(initial_flavor_ratio, Enu, LVATERM = np.zeros((3,3), dtype=numpy.complex),
                LVCTERM = np.zeros((3,3), dtype=numpy.complex)):
    final_flavor_ratio = [0.0]*pc.numneu
    osc_prob_array = [OscProb(beta,Enu,LVATERM=LVATERM,LVCTERM=LVCTERM) for beta in range(pc.numneu)]

    for alpha in range(pc.numneu):
        for beta,phi in enumerate(initial_flavor_ratio):
            final_flavor_ratio[alpha] += osc_prob_array[beta][alpha]*phi
    return final_flavor_ratio

def RRR(initial_flavor_ratio, Enu, LVATERM = np.zeros((3,3), dtype=numpy.complex),
                LVCTERM = np.zeros((3,3), dtype=numpy.complex)):
    ffr = FlavorRatio(initial_flavor_ratio,Enu,LVATERM=LVATERM,LVCTERM=LVCTERM)
    return ffr[1]/ffr[0]
def SSS(initial_flavor_ratio, Enu, LVATERM = np.zeros((3,3), dtype=numpy.complex),
                LVCTERM = np.zeros((3,3), dtype=numpy.complex)):
    ffr = FlavorRatio(initial_flavor_ratio,Enu,LVATERM=LVATERM,LVCTERM=LVCTERM)
    return ffr[2]/ffr[1]

def PointToList(p1,p2):
    return [[p1[0],p2[0]],[p1[1],p2[1]]]

def PointFromFlavor(origin,scale,flavor_ratio_list):
    nu_e_vec = np.array([1.,0.])*scale
    nu_mu_vec = np.array([1./2.,np.sqrt(3.)/2.])*scale
    nu_tau_vec = np.array([-1./2.,np.sqrt(3.)/2.])*scale
    fpos = origin + flavor_ratio_list[0]*nu_e_vec + flavor_ratio_list[1]*nu_mu_vec
    return [fpos[0],fpos[1]]

def MakeFlavorTriangle(list_of_flavor_ratios, scale = 8,
                       p = np.array([0.,0.]), save_file = False, PlotPoints = False, PlotTrayectories = False, figure = None, alpha = 1.0,
                       filename = "triangle",icolor = "green", icolormap = "Greens", divisions = 5, initial_flavor_ratio = [1,0,0],
                       term = "a", subdivisions = False, triangle_collection = None, color_scale = "lin", return_fig = True, addtext = "",
                       add_default_text = True, ilw = 1., slw = 0.75, output_format = "eps", inner_line_color = "k", plot_color_bar = False):
    # i will be nice ...
    list_of_flavor_ratios = np.array(list_of_flavor_ratios)

    if figure == None:
        fig = plt.figure(figsize=(scale,scale), frameon = False)
    else:
        fig = figure

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # delete extra lines
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)

    s0 = np.array([1.,0.])*scale
    s1 = np.array([1./2.,np.sqrt(3.)/2.])*scale
    s2 = np.array([1./2.,-np.sqrt(3.)/2.])*scale

    # make triangle outer frame

    plt.plot(*PointToList(p, p+s0), color = "k", lw = 4)
    plt.plot(*PointToList(p, p+s1), color = "k", lw = 2)
    plt.plot(*PointToList(p+s0, p+s1), color = "k", lw = 2)

    # put outer triangle labels

    # ax.text((p+s0*0.5+s0*0.025)[0], (p+s0*0.5-np.array([0,0.15*scale]))[1],r"$\alpha^{\oplus}_{e}$",
    #             horizontalalignment="center",fontsize = 50, zorder = 10)
    # ax.text((p+s1*0.5-0.2*s0)[0], (p+s1*0.5+0.1*s0)[1]+scale*0.1,r"$\alpha^{\oplus}_{\tau}$",
    #             horizontalalignment="center",fontsize = 50, zorder = 10, rotation = 60.)
    # ax.text((p+s1*0.5 + 0.7*s0)[0], (p+s1*0.5 + 0.6*s0)[1]+0.05*scale,r"$\alpha^{\oplus}_{\mu}$",
    #             horizontalalignment="center",fontsize = 50, zorder = 10, rotation = -60)

    ax.text((p+s0*0.5+s0*0.025)[0], (p+s0*0.5-np.array([0,0.15*scale]))[1],r"$f_{e,\oplus}$",
                horizontalalignment="center",fontsize = 50, zorder = 10)
    ax.text((p+s1*0.5-0.2*s0)[0], (p+s1*0.5+0.1*s0)[1]+scale*0.1,r"$f_{\tau, \oplus}$",
                horizontalalignment="center",fontsize = 50, zorder = 10, rotation = 60.)
    ax.text((p+s1*0.5 + 0.7*s0)[0], (p+s1*0.5 + 0.6*s0)[1]+0.05*scale,r"$f_{\mu, \oplus}$",
                horizontalalignment="center",fontsize = 50, zorder = 10, rotation = -60)

    # construct triangle grid
    fsl = 35
    for i in range(divisions+1):
        subsize = 1./float(divisions)

        ax.text((p+s0*subsize*float(i))[0], (p+s0*subsize*float(i)-np.array([0,0.05*scale]))[1],str(i*subsize),
                horizontalalignment="center",fontsize = fsl)
        plt.plot(*PointToList(p+s0*subsize*float(i), p+s1+s2*subsize*float(i)), color = inner_line_color, lw = ilw, ls = "dashed", zorder = -1)
        ax.text((p+s1-s1*subsize*float(i)-np.array([0.06*scale,0.0]))[0], (p+s1-s1*subsize*float(i))[1],str(i*subsize),
                horizontalalignment="center",fontsize = fsl)
        plt.plot(*PointToList(p+s0*subsize*float(divisions-i), p+s1-s1*subsize*float(i)), color = inner_line_color, lw = ilw, ls = "dashed", zorder = -1)

        ax.text((p+s1+s2*subsize*float(i)+np.array([0.05*scale,0.0]))[0], (p+s1+s2*subsize*float(i))[1],str((divisions-i)*subsize),
                horizontalalignment="center",fontsize = fsl)
        plt.plot(*PointToList(p+s1*subsize*float(divisions-i), p+s1+s2*subsize*float(i)), color = inner_line_color, lw = ilw, ls = "dashed", zorder = -1)

        if subdivisions and i < divisions:
            plt.plot(*PointToList(p+s0*subsize*float(i+0.5), p+s1+s2*subsize*float(i+0.5)), color = inner_line_color, lw = slw, ls = "dotted", zorder = -1)
        if subdivisions and i > 0:
            plt.plot(*PointToList(p+s0*subsize*float(divisions-(i-0.5)), p+s1-s1*subsize*float(i-0.5)), color = inner_line_color, lw = slw, ls = "dotted", zorder = -1)
            plt.plot(*PointToList(p+s1*subsize*float(divisions-(i-0.5)), p+s1+s2*subsize*float(i-0.5)), color = inner_line_color, lw = slw, ls = "dotted", zorder = -1)


    # plot triangle collection
    if (triangle_collection != None):
        # get total number of points
        total_points = float(sum([ triangle.number_of_points for triangle in triangle_collection]))
        max_points = float(max([ triangle.number_of_points for triangle in triangle_collection]))
        color_map = plt.get_cmap(icolormap)
        for triangle in triangle_collection:
            if triangle.number_of_points > 0:
                xx,yy = zip(*triangle.coordinates)
                if color_scale == "lin":
                    plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(0.75), alpha = np.sqrt(float(triangle.number_of_points)/max_points))
                    #plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(float(triangle.number_of_points)/max_points), alpha = alpha)
                elif color_scale == "log":
                    plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(0.75), alpha = (np.log10(float(triangle.number_of_points))/np.log10(max_points)))
                    #plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(0.7), alpha = (np.log10(float(triangle.number_of_points))/np.log10(max_points))**(2./3.))
                    #plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(np.log10(float(triangle.number_of_points))/np.log10(max_points)), alpha = alpha)
                    #plt.fill(xx,yy,lw = 0., zorder = -0.8, color = color_map(np.log10(float(triangle.number_of_points))))
                else : 
                    raise NameError('Error. Love CA.')

        if(plot_color_bar):
            # the color bar magic
            # location set on 0 to 1 scales.
            left = 0.1
            bottom = -0.25
            width = 0.8
            height = 0.025
            cbaxes = fig.add_axes([left,bottom,width,height])
            if color_scale == "lin":
                norm = mpl.colors.Normalize(vmin = 0., vmax = max_points)
            elif color_scale == "log":
                norm = mpl.colors.Normalize(vmin = 0., vmax = 1.0)
            else :
                raise NameError('Error. Love CA.')
            mpl.rcParams.update({'font.size': 10})
            triangle_colorbar = mpl.colorbar.ColorbarBase(cbaxes, cmap = color_map, norm = norm,
                                                          orientation = "horizontal", spacing = "proportional",
							 # )
                                                          format ='%.0e')
            cbaxes.set_xlabel("Likelihood", fontsize = 12)

    # plot flavor ratio points
    if PlotTrayectories :
        if len(list_of_flavor_ratios.shape) == 3 :
            for flavor_ratio_l in list_of_flavor_ratios:
                flv_ratio_coords = map(lambda f:PointFromFlavor(p,scale,np.array(f)),flavor_ratio_l)
                xc, yc = zip(*flv_ratio_coords)
                plt.plot(xc,yc, color = icolor,
                        ms = 10, linewidth = 4, zorder = 0)
        elif len(list_of_flavor_ratios.shape) == 2 :
            flv_ratio_coords = map(lambda f:PointFromFlavor(p,scale,np.array(f)),list_of_flavor_ratios)
            xc, yc = zip(*flv_ratio_coords)

            plt.plot(xc,yc, color = icolor,
                    ms = 10, linewidth = 4, zorder = 0)
        else:
            raise NameError('Check your input flavor list array and the joined flag. Love CA.')
    elif PlotPoints:
        if len(list_of_flavor_ratios.shape) !=2 :
            print len(list_of_flavor_ratios.shape)
            raise NameError('Check your input flavor list array and the joined flag. Love CA.')
        for i,flavor_ratio in enumerate(list_of_flavor_ratios):
            if len(icolor) != len(list_of_flavor_ratios):
                icolor_ = icolor
            else:
                icolor_ = icolor[i]
            plt.plot(*PointFromFlavor(p,scale,np.array(flavor_ratio)), color = icolor_,
                     marker = 'o', ms = 10, linewidth = 0,
                     markeredgecolor=None,markeredgewidth=0.1, zorder = 1000)

    # put back color scale axis
    if add_default_text:
        ax.text((s0/5.+0.9*s1)[0],(s0/5.+0.9*s1)[1],
                "LV "+term+"-term scan with\n $\ \phi_e:\ \phi_\\mu:\ \phi_\\tau = "+str(initial_flavor_ratio[0])+":\ "+str(initial_flavor_ratio[1])+":\ "+str(initial_flavor_ratio[2])+"$"+" \n "+addtext,
                fontsize = 20)

    if(save_file):
        # save figure
        plt.savefig("./plots/"+filename+"."+output_format, dpi = 600, bbox_inches='tight')
    else:
        # show figure
        plt.show()
    if return_fig:
        return fig


def s_bario(p,p0,p1,p2):
    return (p0[1]*p2[0] - p0[0]*p2[1] + (p2[1] - p0[1])*p[0] + (p0[0] - p2[0])*p[1])

def t_bario(p,p0,p1,p2):
    return (p0[0]*p1[1] - p0[1]*p1[0] + (p0[1] - p1[1])*p[0] + (p1[0] - p0[0])*p[1])

def IsInTriangle(p,p0,p1,p2,area):
    s = s_bario(p,p0,p1,p2)
    t = t_bario(p,p0,p1,p2)
    #print s,t,2.0*area - s - t
    return  s >= -1.e-15 and  t >= -1.0e-15 and s+t <= 2.0*area


class Triangle:
    coordinates = []
    area = 0.0
    number_of_points = 0.0
    n_t = 0
    i = 0 
    j = 0
    orientation = ""
    
    def IsPointIn(self,point):
        p0 = self.coordinates[0]
        p1 = self.coordinates[1]
        p2 = self.coordinates[2]
        return IsInTriangle(point,p0,p1,p2,self.area)


def GenerateTriangles(scale, divisions, p = np.array([0.,0.])):
    s0 = np.array([1.,0.])*scale/float(divisions)
    s1 = np.array([1./2.,np.sqrt(3.)/2.])*scale/float(divisions)
    s2 = np.array([1./2.,-np.sqrt(3.)/2.])*scale/float(divisions)
    
    area = np.sqrt(3)*(LA.norm(s0)/2.0)**2
    
    n_t = 0
    
    triangle_collection = []
    for i in range(divisions):
        for j in range(divisions-i):
            lower_triangle = Triangle()

            p0_l = p + i*s0 + j*s1 
            p1_l = p0_l + s0
            p2_l = p0_l + s1
            
            lower_triangle.coordinates = [p0_l,p1_l,p2_l]
            lower_triangle.n_t = n_t
            lower_triangle.i = i
            lower_triangle.j = j
            lower_triangle.orientation = "L"
            lower_triangle.area = area
            
            n_t += 1
            # append to triangle collection
            triangle_collection.append(lower_triangle)
            
            upper_triangle = Triangle()

            p0_u = p2_l
            p1_u = p1_l
            p2_u = p1_l + s1
            
            upper_triangle.coordinates = [p0_u,p1_u,p2_u]
            upper_triangle.n_t = n_t
            upper_triangle.i = i
            upper_triangle.j = j
            upper_triangle.orientation = "U"
            upper_triangle.area = area
            
            n_t += 1
            # append to triangle collection
            triangle_collection.append(upper_triangle)
    return triangle_collection

def AddPointToTriangleCollectionLegacy(flavor_ratio, triangle_collection,
                  p = np.array([0.,0.]), scale = 8, divisions = 10):
    point = PointFromFlavor(p,scale,np.array(flavor_ratio))
    electron = 0; tau = 2;
    # the silly way
    for triangle in triangle_collection:
        if(triangle.IsPointIn(point)):
            triangle.number_of_points += 1.

def AddPointToTriangleCollection(flavor_ratio, triangle_collection,
                  p = np.array([0.,0.]), scale = 8, divisions = 10):
    point = PointFromFlavor(p,scale,np.array(flavor_ratio))
    electron = 0; muon = 1; tau = 2;
    # the smart way
    u_i = int(flavor_ratio[electron]*float(divisions))
    u_j = int(flavor_ratio[muon]*float(divisions))
    index = u_i*(2*divisions-u_i+1) + 2*u_j
    if triangle_collection[index].IsPointIn(point):
        triangle_collection[index].number_of_points += 1.
    else:
        triangle_collection[index+1].number_of_points += 1.
# legacy
    #elif triangle_collection[index+1].IsPointIn(point):
    #    triangle_collection[index+1].number_of_points += 1.
    #else:
    #    print "Math error."
    #    print point, "\n",u_i, u_j, "\n", triangle_collection[index].coordinates, "\n", triangle_collection[index+1].coordinates
    #    raise NameError("Error triangle location math")

class AnarchySampling:
    def __init__(self, n_sample, LV_scale_1, LV_scale_2, term):
        self.n_sample = n_sample
        self.th12_sample = np.arcsin(np.sqrt(np.random.uniform(0.,1., size=n_sample)))
        self.th13_sample = np.arccos(np.sqrt(np.sqrt(np.random.uniform(0.,1., size=n_sample))))
        self.th23_sample = np.arcsin(np.sqrt(np.random.uniform(0.,1., size=n_sample)))
        self.delta_sample = np.random.uniform(0.,2.*np.pi, size=n_sample)

        self.LV_scale_1 = LV_scale_1
        self.LV_scale_2 = LV_scale_2

        self.term = term

def GenerateFlavorRatioPoints(Initial_Flavor_Ratio, SamplingObject, gamma = 2.0,
                              Log10Emax = 7., Log10Emin = 4.0, Epoints = 30,
                              save_list = False, save_avg = True):
    flavor_tray_list = []
    flavor_avg_list = []

    # energy things

    Erange = np.logspace(Log10Emin,Log10Emax,Epoints) # in GeV
    Emin = Erange[0]
    Emax = Erange[-1]

    if gamma == 1 or gamma == 1.0:
        spectral_normalization = np.log(Emax)-np.log(Emin)
    else:
        spectral_normalization = (Emax**(1.-gamma) - Emin**(1.-gamma))/(1.-gamma)

    spectral_function = lambda Enu: Enu**(-gamma)/spectral_normalization

    # loop over random parameters
    for i in range(SamplingObject.n_sample):
        lv_term = MT.LVP()

        lv_term.th12 = SamplingObject.th12_sample[i]
        lv_term.th13 = SamplingObject.th13_sample[i]
        lv_term.th23 = SamplingObject.th23_sample[i]
        lv_term.delta1 = SamplingObject.delta_sample[i]

        lv_term.LVS21 = SamplingObject.LV_scale_1
        lv_term.LVS31 = SamplingObject.LV_scale_2

        lv_term.Refresh()

        LVTERM = MT.LVTerm(lv_term);

        if SamplingObject.term == "a":
            flavor_ratio_list = np.array(map(lambda Enu : FlavorRatio(Initial_Flavor_Ratio, Enu*pc.GeV, LVATERM = LVTERM), Erange))
        elif SamplingObject.term == "c":
            flavor_ratio_list = np.array(map(lambda Enu : FlavorRatio(Initial_Flavor_Ratio, Enu*pc.GeV, LVCTERM = LVTERM), Erange))
        else :
            raise NameError('Only a or c term.'+ str(term))

        if save_avg: 
            if Epoints != 1:
                flavor_avg = [0.]*lv_term.numneu
                for alpha in range(lv_term.numneu):
                    #inter = interpolate.interp1d(Erange,flavor_ratio_list[:,alpha])
                    inter = interpolate.UnivariateSpline(Erange,flavor_ratio_list[:,alpha])
                    flavor_avg[alpha] = integrate.quad(lambda Enu : inter(Enu)*spectral_function(Enu),
                                                       Emin,Emax, limit = 75, epsrel = 1e-2, epsabs = 1.0e-2)[0]
                    #flavor_avg[alpha] = integrate.quadrature(lambda Enu : inter(Enu)*spectral_function(Enu),
                    #                                         Emin,Emax, maxiter = 75, rtol = 1e-3, tol = 1.e-3)[0]
                flavor_avg_list.append(flavor_avg)
            else:
                flavor_avg = flavor_ratio_list[0]
                flavor_avg_list.append(flavor_avg)

        if save_list:
            flavor_tray_list.append(flavor_ratio_list)

    if save_list and save_avg:
        return flavor_tray_list, flavor_avg_list
    elif save_list:
        return flavor_tray_list
    elif save_avg:
        return flavor_avg_list
    else :
        print "Math is broken."
        return None
