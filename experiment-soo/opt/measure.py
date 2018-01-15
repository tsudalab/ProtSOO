#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI lolicon.jimmy@gmail.com
# REMEMBER: to update the pso answer file location
# REMEMBER: the protein angle calculation is specified for the LYZ protein, please update the two C atom id for another case.

from __future__ import print_function       # for python 2.7
import sys, os, argparse, math, re

# These are the "Color Blind 10" colors as RGB.  
#blind10=[(255,128,14),(171,171,171),(95,158,209),(89,89,89),(0,107,164),(255,188,121),(207,207,207),(200,82,0),(162,200,236),(137,137,137)]
blind10=[(0,107,164),(200,82,0),(255,128,14),(137,137,137),(171,171,171),(162,200,236),(89,89,89),(255,188,121),(95,158,209),(207,207,207)]
tableau5=[(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189)]
tableau10=[(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)] 
tableau10l=[(174, 199, 232), (255, 187, 120), (152, 223, 138), (255, 152, 150), (197, 176, 213), (196, 156, 148), (247, 182, 210), (199, 199, 199), (219, 219, 141), (158, 218, 229)] 
tableau10m=[(114, 158, 206), (255, 158, 74), (103, 191, 92), (237, 102, 93), (173, 139, 201), (168, 120, 110), (237, 151, 202), (162, 162, 162), (205, 204, 93), (109, 204, 218)] 
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
color_scheme = {"blind10": blind10, "tableau20": tableau20, "tableau10": tableau10, "tableau5": tableau5, "tableau10light": tableau10l, "tableau10medium": tableau10m}

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
# Set the default color cycle
#mpl.rcParams['axes.color_cycle'] = blind10
# Alternately, we could use rc:
# mpl.rc('axes', color_cycle=['r','k','c'])
import numpy as np

import pdb


class measure:
    def __init__(self):
        self.loacts = ['helping', 'measure', 'rmwtr', 'energy', 'resh', 'sumup', 'avgex']
        self.parser = argparse.ArgumentParser(description="This is measure.py")
        self.parser.add_argument("act", help="run something main/etc...", choices=self.loacts)
        self.parser.add_argument("--pdb", help="set the input protein pdb filename.", required=False, default="system.pdb")
        self.parser.add_argument("--smooth", help="smooth the curve using WINDOW size average.", required=False, default=0, type=int, metavar="WINDOW")
        self.parser.add_argument("--xvg", help="load the xvg files.", nargs="*", default=[], type=str)
        self.parser.add_argument("--bvg", help="load the xvg files but plot as bold line. (it is the first to be processing)", nargs="*", default=[], type=str)
        self.parser.add_argument("--txt", help="load the txt files. (txt is plotting before xvg)", nargs="*", default=[], type=str)
        self.parser.add_argument("--bxt", help="load the txt files but plot as bold line. (it is the second to be processing)", nargs="*", default=[], type=str)
        self.parser.add_argument("--lin-width", help="set the width of normal lines.", default=1.0, type=float)
        self.parser.add_argument("--bxt-width", help="set the width of the bold lines. it is multiples of normal width.", default=2.0, type=float)
        self.parser.add_argument("--xten", help="x asix divided by 10^x(ten base) where x is the input, higher priority then xuntil", type=int, default=0)
        self.parser.add_argument("--xuntil", help="display just until this number (showing on the plotting).", type=float, default=0)
        self.parser.add_argument("--yoffset", help="offset all the y values by adding this number.", type=float, default=0)
        self.parser.add_argument("--ylim", help="manually choose the y asix limit.", nargs=2, type=float)
        self.parser.add_argument("--savetxt", help="save the xvg files into numpy file format. Only available in action: energy", default=False, action='store_true')
        self.parser.add_argument("--savepdf", help="save plot as pdf file.", default=False, action='store_true')
        self.parser.add_argument("--txtsize", help="set the text size of axes label.", default=24, type=int)
        self.parser.add_argument("--nbins", help="set the nbins for matplotlib.", default=7, type=int)
        self.parser.add_argument("--xlabel", help="set the x axis label.", default="", type=str)
        self.parser.add_argument("--ylabel", help="set the y axis label.", default="", type=str)
        legend_group = self.parser.add_mutually_exclusive_group()
        legend_group.add_argument("--no-legend", help="do NOT display the legend, useful when plotting a set of graphics.", default=False, action='store_true')
        legend_group.add_argument("--right-legend", help="show the legend on right hand side.", default=False, action='store_true')
        legend_group.add_argument("--top-legend", help="show the legend on top.", default=False, action='store_true')
        self.parser.add_argument("--legends", help="given the legend name manually.", type=str, nargs="*")
        self.parser.add_argument("--name", help="put a text at left top corner.", type=str, default="")
        self.parser.add_argument("--ltxsize", help="set the legend text font size.", default=12, type=int)
        self.parser.add_argument("--save", help="equivalent to --savetxt --savepdf", default=False, action='store_true')
        self.parser.add_argument("--vs", help="concatenate the arrays vertically, for joining different types of energy. only for xvg file, and default is **horizontally**.", default=False, action='store_true')
        self.parser.add_argument("-z", "--compress", help="compress the txt files.", action="store_const", const=".gz", default="", dest="z")
        self.parser.add_argument("-b", "--batch", help="does not show the plot graph, useful for batch processing.", default=False, action="store_true", dest="b")
        self.parser.add_argument("-o", "--output", help="set the output pdf file name.", default="energy.pdf", dest="o")
        #self.parser.add_argument("--color-scheme", help="set the color scheme", default="blind10", choices=["blind10", "tableau20"])
        self.parser.add_argument("--color-scheme", help="set the color scheme", default="tableau10", choices=color_scheme.keys())
        self.args = self.parser.parse_args()

        # backup the system stdout because pymol will override it.
        self.pout = sys.stdout

        # pymol environment TODO: only for MacOS, pymol was installed from MacPort
        self.moddir='/usr/lib64/python2.7/site-packages'
        self.moddir='/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/'
        sys.path.insert(0, self.moddir)
        os.environ['PYMOL_PATH'] = os.path.join(self.moddir, 'pymol/pymol_path')

        # pymol launching
        import pymol
        #pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
        pymol.pymol_argv = ['pymol','-qc'] 
        pymol.finish_launching()
        self.cmd = pymol.cmd
        self.cmd.set("retain_order", 1)
        self.cmd.set("pdb_retain_ids", 1)

        #self.cmd.do('print("Hajimaru yo ~♪")')
        #self.hello("はじまるよ～♪ \n")     

        self.values = []
        self.legend = []
        self.crank = -1
        self.cthem = color_scheme[self.args.color_scheme]
        self.clenf = len(self.cthem)
        self.curc = ""  #current color code 

        self.ylimp = -7777777
        self.ylimn = +7777777
        


    def nxc(self):
        """get the index of next color"""
        self.crank += 1
        self.crank %= self.clenf
        #convert into matplotlib accepted format
        r = [ x/255.0 for x in self.cthem[self.crank] ]
        return r

    def hello(self, msg="はじまるよ～♪ "):
        #TODO convert hello print-out to logger based
        self.pout.write(msg)
        self.pout.write("\n")
        #logger.debug("[hello] {}".format(msg))
        #logger.debug("hello function is deprecated.")

    def load_file(self, filenm):
        self.filenm = ""
        if not os.path.isfile(filenm):
            return 0

        self.cmd.delete("*")    #Clean everything first
        try:
            self.cmd.load(filenm, "system")
            self.filenm = filenm
            self.hello("Loaded file: {}".format(self.filenm))
            self.cmd.select("surface", "resn PTS or resn SUB")
            self.cmd.select("protein", "not resn SOL and not resn PTS and not resn SUB")
            self.hello("Selected surface and protein objects.")
        except:
            self.hello("Loading file ERROR!!!")
            return 0

        return 1

    def cutSuf(self):
        self.xyzedges("resn PTS and name C1")
        x_off_pts = (self.xedgep-self.xedgen)*0.2
        y_off_pts = (self.yedgep-self.yedgen)*0.2

        self.xyzedges("resn SUB")
        x_off_sub = (self.xedgep-self.xedgen)*0.2
        y_off_sub = (self.yedgep-self.yedgen)*0.2

        remove_count = 0

        try:
            atoms = self.cmd.get_model("resn PTS or resn SUB")
        except:
            self.hello("Can not find model({})".format(s))
            return 0

        for atom in atoms.atom:
            try:
                x = float(atom.coord[0])
                y = float(atom.coord[1])
                z = float(atom.coord[2])
                resn = atom.resn
                resi = atom.resi
                atid = atom.id
                name = atom.name
            except:
                x=y=z=0
                name = rean = atid = resi = ""
            if resn == "PTS" and name == "C1" and \
                    (x < self.xedgen+x_off_pts or x > self.xedgep-x_off_pts or \
                    y < self.yedgen+y_off_pts or y > self.yedgep-y_off_pts):
                        self.cmd.remove("resn PTS and resi {}".format(resi))
                        remove_count+=1
            elif resn == "SUB" and \
                    (x < self.xedgen+x_off_sub or x > self.xedgep-x_off_sub or \
                    y < self.yedgen+y_off_sub or y > self.yedgep-y_off_sub):
                        self.cmd.remove("id {}".format(atid))
                        remove_count+=1
        self.hello("Removed {} lines of PTS/SUB".format(remove_count))

        


    def cmxlength(self, p="protein"):
        """measure the maximum length but without print out"""
        pass

    def mxlength(self, p="protein"):
        """measure the maximum length of protein"""
        try:
            atoms = self.cmd.get_model(p)
        except:
            #self.hello("Can not find model({}).".format(p))
            return 1

        maxplen = 0
        maxatom = maxbtom = 0  #in case of no atoms in the model
        for atoma in atoms.atom:
            try:
                a = [float(i) for i in atoma.coord[:3]]
            except ValueError:
                a = [0.0] * 3
            for atomb in atoms.atom:
                try:
                    b = [float(i) for i in atomb.coord[:3]]
                except ValueError:
                    b = [0.0] * 3
                deltasq = [(i-j)**2 for i,j in zip(a, b)]
                distanc = math.sqrt(sum(deltasq))
                if distanc > maxplen:
                    maxplen = distanc
                    maxatom = "{} {} - {}({})".format(atoma.resn, atoma.resi, atoma.name, atoma.index)
                    maxatid = atoma.index
                    maxatrd = "{} {}".format(atoma.resn, atoma.resi)
                    maxbtom = "{} {} - {}({})".format(atomb.resn, atomb.resi, atomb.name, atomb.index)
                    maxbtid = atomb.index
                    maxbtrd = "{} {}".format(atomb.resn, atomb.resi)
        #self.hello("Found the maximum length {} between {} and {}.".format(maxplen, maxatom, maxbtom))
        return maxplen

    def xyzedges(self, s="surface"):
        """measure the X/Y/Z edges and length."""
        #PyMol already implemented get_extent function to get the X Y Z min/max values. see: http://www.pymolwiki.org/index.php/Get_extent
        try:
            atoms = self.cmd.get_model(s)
        except:
            #self.hello("Can not find model({})".format(s))
            return 0

        self.xedgep = -999999999
        self.xedgen = +999999999
        self.yedgep = -999999999
        self.yedgen = +999999999
        self.zedgep = -999999999
        self.zedgen = +999999999

        for atom in atoms.atom:
            try:
                x = float(atom.coord[0])
                y = float(atom.coord[1])
                z = float(atom.coord[2])
            except ValueError:
                x=y=z=0
            if x > self.xedgep:
                self.xedgep = x
            if x < self.xedgen:
                self.xedgen = x
            if y > self.yedgep:
                self.yedgep = y
            if y < self.yedgen:
                self.yedgen = y
            if z > self.zedgep:
                self.zedgep = z
            if z < self.zedgen:
                self.zedgen = z

        #self.hello("Found the {} x axis edge: [{}, {}]; y axis edge: [{}, {}]; z axis edge: [{}, {}].".format(s, self.xedgen, self.xedgep, self.yedgen, self.yedgep, self.zedgen, self.zedgep))
        #self.hello("x edge length: {}; y edge length: {}; z edge length: {};".format(self.xedgep-self.xedgen,self.yedgep-self.yedgen, self.zedgep-self.zedgen))
        self.centerp = [(j-i)/2 +i for i,j in [(self.xedgen, self.xedgep), (self.yedgen, self.yedgep), (self.zedgen, self.zedgep)]]
        #self.hello("The center of {} is: [{}, {}, {}]".format(s, self.centerp[0], self.centerp[1], self.centerp[2]))


    def dElipsoid(self):
        #Now, draw the Elipsoid.
        # The semi-axes of Elipsoid
        self.epsax = [(j-i)/2 for i,j in [(self.xedgen, self.xedgep), (self.yedgen, self.yedgep), (self.zedgen, self.zedgep)]]
        """
        pointa = self.centerp[:]
        pointa[0] = self.xedgen
        pointa = [self.xedgep, self.yedgep, self.zedgep]
        pointb = [i-j*1.5 for i,j in zip(pointa, self.centerp)]
        self.cmd.do('print("point B is: {}")'.format(pointb))
        one = [i**2/j**2 for i,j in zip(pointb, self.epsax)]
        self.cmd.do('print("The one is: {} (sum:{})")'.format(one, sum(one)))
        #Thus, if one is <= 1 it means the point is inside the Elipsoid.
        """

    def removeW(self, s="all"):
        """
        removing all overlap water SOL solvent  residues.
        the top boundary of surface can be auto detected.
        """
        try:
            atoms = self.cmd.get_model(s)
        except:
            self.hello("Can not find model({})".format(s))
            return 0

        self.xyzedges("surface")
        surface_top = self.zedgep
        self.xyzedges("protein")
        #We've got x y z edges, for the Elipsoid.

        self.dElipsoid()
        #We've got the Elipsoid three semi-axis, epsax.

        water_count = 0
        remove_count = 0

        for atom in atoms.atom:
            try:
                x = float(atom.coord[0])
                y = float(atom.coord[1])
                z = float(atom.coord[2])
                resn = atom.resn
                resi = atom.resi
                atid = atom.id
            except:
                x=y=z=0
                resi = rean = atid = ""
            if resn == "SOL":
                r = [i**2/j**2 for i,j in zip([x,y,z], self.epsax)]
                #if z <= surface_top+1 or r <= 1:
                if z <= surface_top+1:
                    self.cmd.remove("resn SOL and resi {}".format(resi))
                    remove_count+=1
                else:
                    water_count+=1
        self.hello("Removed {} elements.".format(remove_count))
        self.hello("We have got {} number of water elements.".format(water_count))
        """
            if resn != "PTS" and resn != "SUB" and resn != "SOL":
                #This is the protein.
                pass
            else:
                self.cmd.do('print("resn: {} atom id: {}")'.format(resn, atid))
        """

    def save_file(self):
        self.hello("Saving file...")
        self.cmd.save(self.filenm)
        self.hello("Saved: {}".format(self.filenm))

    def prota(self, ida="id 5", idb="id 1958"):
        """
        measure the eigenvectors of protein model.
        Target id 5 and 1958 is specified for LYZ protein."""
        self.hello("measure the protein and axises angle")
        CA = []
        CB = []
        try:
            target = ida
            atoms = self.cmd.get_model(target)
            #[ atom for atom in atoms]
            atom = atoms.atom[0]
            CA = np.array(atom.coord[:3])

            target = idb
            atoms = self.cmd.get_model(target)
            atom = atoms.atom[0]
            CB = np.array(atom.coord[:3])
        except IndexError:
            self.hello("No atom {} or invalid coord".format(target))
        except:
            self.hello("Can not find atom {}!".format(target))
            raise

        if CA != [] and CB != []:
            self.hello("CA: {}".format(CA))
            self.hello("CB: {}".format(CB))
            protein_v = CB - CA
            self.hello("protein vector: {}".format(protein_v))
            cos = np.dot(CA, CB)/(np.linalg.norm(CA)*np.linalg.norm(CB))
            deg = np.degrees(np.arccos(cos))
            self.hello("c vector angle: {}".format(deg))

            xi = [1, 0, 0]
            yi = [0, 1, 0]
            zi = [0, 0, 1]

            cosx = np.dot(protein_v, xi)/(np.linalg.norm(protein_v)*np.linalg.norm(xi))
            cosy = np.dot(protein_v, yi)/(np.linalg.norm(protein_v)*np.linalg.norm(yi))
            cosz = np.dot(protein_v, zi)/(np.linalg.norm(protein_v)*np.linalg.norm(zi))
            self.degx = np.degrees(np.arccos(cosx))
            self.degy = np.degrees(np.arccos(cosy))
            self.degz = np.degrees(np.arccos(cosz))
            self.hello("three axis angle is, X:{}, Y:{}, Z:{}".format(self.degx, self.degy, self.degz))
            #return (degx, degy, degz)

    def reshi(self):
        self.hello("First, select the protein & surface object.")
        self.hello("Find the edges of surface.")
        self.xyzedges("surface")
        ztop = self.zedgep
        self.hello("Second, loop all residues and calculate the minidistance between the surface and the residue.")
        self.resl = {}
        try:
            atoms = self.cmd.get_model("protein")
        except:
            self.hello("Could not find protein model.")
            return 0

        for atoma in atoms.atom:
            try:
                reid = int(atoma.resi)
                z  = float(atoma.coord[2])
            except ValueError:
                reid = -1
                z = 0.0
            except IndexError:
                z = 0.0
                self.hello("[reshi] Can not get Z value!.")
            if reid != -1:
                if reid not in self.resl or z < self.resl[reid]:
                    self.resl[reid] = z

        lx = []
        ly = []
        #ASSUMED: minh is initialy large enough, the coord field in PDB is just 8 columns
        minh = 999999999
        mini = 0
        for i in self.resl:
            #self.hello("resl[{}]: {}".format(i, self.resl[i]))
            lx += [i]
            ly += [self.resl[i]]
            if self.resl[i] < minh:
                minh = self.resl[i]
                mini = i
        nlx = np.array(lx)
        nly = np.array(ly)

        fig = plt.figure()
        fres = plt.plot(nlx, nly-ztop, 'k') 
        #plt.hlines(minh-ztop, nlx[0], nlx[-1], color='b')
        self.hello("found the mini distance #{}@{}".format(mini, minh))

        #plt.show()

        self.hello("Thrid, plot the graph.")
        plt.ylabel(u"The minimum distance to the surface of each residues (Å)")
        #plt.ylabel(u"Unit Å")
        plt.savefig("diagram_rh.pdf", dpi=300)

    def rmsd(self, a, b):
        if len(a) != len(b):
            m.hello("vector A and B must have same size!!!")
            return 1

        return np.sqrt(np.mean(np.square(a-b)))

    def avgex(self):
        """average all numbers for each X"""
        num_txt = 0
        max_len = 0
        LA = []
        for txt in self.args.txt:
            self.hello("Loading: {}...".format(txt))
            L = np.loadtxt(txt)
            LA.append(L)
            if max_len < len(L):
                max_len = len(L)
            num_txt += 1
        avgnum = []
        median = []
        for i in range(max_len):
            ct = 0
            sm = 0
            avgnp  = []
            for nbs in LA:
                if i < len(nbs):
                    sm+=nbs[i]
                    ct+=1
                    avgnp.append(nbs[i])
            avgnp = np.array(avgnp)
            median.append(np.median(avgnp))
            if ct == 0:
                #It should never happen.
                avgnum.append(0)
            else:
                avgnum.append(sm/ct)
        avgnum = np.array(avgnum)
        np.savetxt("average-num.txt.gz", avgnum)
        self.hello("Saved the average numbers as average-num.txt.gz")
        median = np.array(median)
        np.savetxt("median-num.txt.gz", median)
        self.hello("Saved the median numbers as median-num.txt.gz")

    def sumup(self):
        """sum up all the number input from text files.
        usually used to add the two energy curves in MD simulations."""
        LS = np.array([])
        for txt in self.args.txt:
            self.hello("Loading :{}...".format(txt))
            L = np.loadtxt(txt)
            try:
                LS = LS + L
            except ValueError:
                LS = L
                #The first time add up LS, the shape is different.
                #However, if the numbers list shape is different, then can't add up.
        np.savetxt("sumup.txt"+self.args.z, LS)
        self.hello("Saved sumup.txt"+self.args.z)


    def resh(self):
        "plot the residues minimum distance graph."
        filenm = "surface_s.pdb"
        filenm = self.args.pdb

        if self.load_file(filenm):
            #self.reshi()
            self.prota()
            self.hello("The target pdb structure angles: {}".format([self.degx, self.degy, self.degz]))
            resl = self.residst()
            resx = np.arange(1, len(resl)+1)

            fig = plt.figure()
            plt.plot(resx, resl, 'k')
            plt.xlabel("residue number")   #№
            plt.ylabel(u"minimum distance to surface of each residues (Å)")
            #plt.xticks(np.arange(min(self.ansx), max(self.ansx)+1, 1.0))
            #plt.xticks(self.ansx)
            plt.savefig("diagram_rxh.pdf", dpi=300)
            self.hello("diagram_rxh.pdf saved")

            L = np.column_stack((resx, resl))
            np.savetxt("diagram_rxh.txt", L, delimiter=" ", fmt="%03d %07.4f")
            self.hello("diagram_rxh.txt saved")
            L = np.loadtxt("diagram_rxh.txt", dtype=[('resn', int), ('dest', float)])
            L.sort(order=['dest'])
            #L.sort(axis=1)
            np.savetxt("diagram_sort.txt", L, delimiter=" ", fmt="%03d %07.4f")
            self.hello("diagram_sort.txt saved")

    def rmwtr(self):
        """remove all water residue under the top boundary of surface."""
        if self.load_file(self.args.pdb):
            self.hello("Load the pdb structure.")
            self.removeW("all")
            self.save_file()
        else:
            self.hello("Can't find the pdb structure.")

    def energy(self):
        """plot the energy curve from energy.xvg g_energy output.
        here pdb argument will be used as the input filename.
        Usually it should be energy.xvg. This is the default name.
        Features:
        1. load multi xvg/txt files.
        2. can concatenate xvg files into one single txt files per column.
        3. save the plot as a pdf file.
        4. supported compression on txt files.
        """
        plt.figure()
        at = mpl.offsetbox.AnchoredText(self.args.name,
                          prop=dict(size=self.args.ltxsize, weight="bold"), frameon=False,
                          loc=2,
                          )
        #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax = plt.gca()
        ax.add_artist(at)
        plt.locator_params(tight=None, nbins=self.args.nbins)
        if self.args.b:
            #TODO: to be fixed a bug, if the file name did not specified the format.
            # png is the default, and it will crash here. since dvi2png not found.
            #latex supporting ONLY available for pdf not viewing.
            #TODO: UPDATED: something in matplotlib, it needs dvipng to work.
            #even saving pdf and not showing the plotting it still reported dvipng not found.
            #Just installed it.
            plt.rc('text', usetex=True)
            #other text format are ONLY meaningful for pdf.
            plt.rc('font', family='serif', size=self.args.txtsize)
            plt.xlabel(self.args.xlabel, fontsize=self.args.txtsize)
            plt.ylabel(self.args.ylabel, fontsize=self.args.txtsize, multialignment='center')
            plt.xticks(size=self.args.txtsize)
            plt.yticks(size=self.args.txtsize)
        for bxt in self.args.bxt:
            self.energytxt(bxt, width=self.args.lin_width*self.args.bxt_width)
        for bvg in self.args.bvg:
            self.energyxvg(bvg, width=self.args.lin_width*self.args.bxt_width)
        for txt in self.args.txt:
            self.energytxt(txt, width=self.args.lin_width)
        for xvg in self.args.xvg:
            self.energyxvg(xvg, width=self.args.lin_width)
        #self.hello("Legends: {}".format(self.legend))
        c=1
        if self.args.save or self.args.savetxt:
            self.hello("saved the values into text file(s).")
            #Using the legend as the txtfile name
            #ASSUMED: the legend is correct, no invalid characters, no duplicates
            #ASSUMED: here the order is same as the parameters order,
            #but we are better to sort by the first column.
            for i in range(len(self.values)):
                if len(self.legend) > i:
                    legend = self.legend[i]
                else:
                    legend = "S{:02}".format(c)
                np.savetxt(legend+".txt"+self.args.z, self.values[i])
                c+=1
        #fig = plt.gcf()
        if self.args.save or self.args.savepdf:
            self.hello("saved the energy graph: {}".format(self.args.o))
            plt.tight_layout()
            plt.savefig(self.args.o, bbox_inches="tight", format="pdf", transparent=True)
            #fig.savefig(self.args.o, bbox_inches="tight")
        #Must show the plot in the end, otherwise after closing the window the pdb file will be empty.
        if not self.args.b:
            plt.show()



    def energytxt(self, txtfile, width=1.0):
        """A txtfile should contain only one array.
        Using the filename as the legend.
        """
        try:
            npvalues = np.loadtxt(txtfile)
        except IOError:
            self.hello("can not open file: {}".format(txtfile))
            return

        legends = [re.sub("\.txt(\.gz)?","",txtfile)]
        self.hello("Loaded numbers from {}".format(txtfile))
        npvt = np.vstack((np.arange(len(npvalues)), npvalues))
        self.legend += legends
        self.energyplot(npvt, self.legend, width)

        if self.args.savetxt:
            np.savetxt("out-"+txtfile, npvt, delimiter="\n")
            self.hello("saved numbers to out-{}".format(txtfile))

    def energyxvg(self, xvgfile, width=1.0):
        npvalues = []
        legends = []
        try:
            with open(xvgfile, "r") as xvg:
                for line in xvg:
                    m = re.search("legend\s*\"(.*?)\"", line)
                    if m:   #This is specified for energy.xvg
                        legends.append(m.group(1))
                        lgnotfd = False
                    if not re.match("^\s*$", line) and re.match("^[0-9+-.e\s]*$", line):
                        values = line.split()
                        fvalues = [ float(x.replace(',','')) for x in values  ] #must be float number matched in re.
                        npvalues.append(fvalues)
        except IOError:
            self.hello("Can not open file: {}".format(xvgfile))


        npvalues = np.array(npvalues)
        #self.hello(str(npvalues.T))
        self.hello("Loaded data {} from {}".format(legends, xvgfile))
        npvt = npvalues.T

        if self.args.vs:
            self.legend += legends
        else:
            self.legend = legends
        if self.args.save or self.args.savetxt:
            try:
                #We do not want the first row, which is incremental numbers
                if self.args.vs:
                    self.values = np.vstack((self.values, npvt[1:]))
                else:
                    self.values = np.hstack((self.values, npvt[1:]))
            except ValueError:
                self.values = np.array(npvt[1:])
            except:
                #TODO: exception case if nothing loaded from the xvg file.
                pass

        self.energyplot(npvt, self.legend, width)

    def energyplot(self, lonum, legends, width=1.0):
        legendss = []
        for i in range(1, len(lonum)):
            X = lonum[0] / 10**self.args.xten
            Y = lonum[i]
            if self.args.smooth:
                Y = self.quick_sliding_mean_0(Y, self.args.smooth)
                #np.savetxt("debug.txt", Y, delimiter="\n")
            if self.args.xuntil:
                idx = 0
                for x in X:
                    #ASSUMED: X is incremental numbers.
                    if x <= self.args.xuntil:
                        idx += 1
                X = X[:idx]
                Y = Y[:idx]
            if self.args.yoffset:
                Y+=self.args.yoffset
            plt.plot(X, Y, color=self.nxc(), linewidth=width)
            legendss.append("s{:02}".format(i))
            if self.args.no_legend:
                self.hello("Color[{}]: {{{}}}".format(self.crank,self.cthem[self.crank]))
            yl = Y.min()
            yh = Y.max()
            yr = np.abs(yh - yl)    
            yh += yr*0.1
            yl -= yr*0.1
            if yh > self.ylimp:
                self.ylimp = yh
            if yl < self.ylimn:
                self.ylimn = yl
            if self.args.ylim and len(self.args.ylim) == 2:
                plt.ylim(self.args.ylim[0], self.args.ylim[1])
            else:
                plt.ylim(self.ylimn, self.ylimp)
        if self.args.legends and len(self.args.legends) >= len(legends):
            legendl = self.args.legends
        elif legends:
            legendl = legends
        else:
            legendl = legendss
        #plt.legend(legendl, loc='upper left', ncol=7, frameon=False, labelspacing=0, columnspacing=0, handletextpad=0)
        numcol = 1
        bbanch = (1.0, 1.0)
        legloc = "upper right"
        if self.args.right_legend:
            #bbanch = (1.5, 1.0)
            legloc = "upper left"
        if self.args.top_legend:
            numcol = 7
        if not self.args.no_legend:
            plt.legend(legendl, loc=legloc, bbox_to_anchor=bbanch, ncol=numcol, frameon=False, \
                    labelspacing=0, columnspacing=0, handletextpad=0, prop={'size': self.args.ltxsize})


    def measure(self):
        """
        measure the protein features, including:
        mxlength - the maximum length of protein;
        xyzedges - the x/y/z edges and length of the surface and the whole system;
        prota    - measure the eigenvectors of protein;
        """
        if self.load_file(self.args.pdb):
            self.hello("Loaded the {} structure.".format(self.args.pdb))
            self.mxlength("protein")
            self.xyzedges("all")
            self.xyzedges("surface")
            ztops = self.zedgep
            self.xyzedges("protein")
            zbotp = self.zedgen
            self.hello("Minimum distance from protein to surface is: {} - {} = {}".format(zbotp, ztops, zbotp - ztops))
            self.prota("protein")
        else:
            self.hello("Can't find the pdb structure.")

    def sliding_mean(self, data_array, window=5):  
        data_array = np.array(data_array)  
        new_list = []  
        for i in range(len(data_array)):  
            indices = range(max(i - window + 1, 0),  
                    min(i + window + 1, len(data_array)))  
            avg = 0  
            for j in indices:  
                avg += data_array[j]  
            avg /= float(len(indices))  
            new_list.append(avg)  

        return np.array(new_list)  

    def quick_sliding_mean(self, data_array, window=10):
        """
        Smoothing data can help to understand the curves
        Here, the implementation running time is just O(N)
        python can not just return the data_array
        Maybe because of the pointer/reference problem.
        """
        data_array = np.array(data_array) #without this line, we can't have the correct result.
        pra = 0
        avg = 0
        asum = 0
        for i in range(len(data_array)):
            numb = data_array[i]
            if i < window:
                asum += numb
                avg = asum / (i+1)
            else:
                pra = last_numb/window
                avg = avg - pra + numb/window
            data_array[i] = avg    #this make the result incorrect, the i(=204834) is too large?
            last_numb = numb
        return data_array

    def quick_sliding_mean_1(self, data_array, window=10):
        """
        Smoothing data can help to understand the curves
        Here, the implementation running time is just O(N)
        Also used O(N) memory
        """
        #data_array = np.array(data_array)
        pra = 0
        avg = 0
        asum = 0
        for i in range(len(data_array)):
            numb = data_array[i]
            if i < window:
                asum += numb
                avg = asum / (i+1)
            else:
                pra = last_numb/window
                avg = avg - pra + numb/window
            data_array[i] = avg
            last_numb = numb
        return data_array

    def quick_sliding_mean_0(self, data_array, window=10):
        data_array = np.array(data_array)
        new_list = []
        pra = 0
        avg = 0
        asum = 0
        for i in range(len(data_array)):
            if i < window:
                asum += data_array[i]
                avg = asum / (i+1)
            else:
                pra = data_array[i-window]/window
                avg = avg - pra + data_array[i]/window
            new_list.append(avg)
        return np.array(new_list)

    def residst(self):
        self.xyzedges("surface")
        ztop = self.zedgep
        #self.hello("Got the surface z top: {}".format(ztop))

        resd = {}
        try:
            atoms = self.cmd.get_model("protein")
        except:
            #self.hello("Could not find protein model.")
            return 0

        for atoma in atoms.atom:
            try:
                reid = int(atoma.resi)
                z = float(atoma.coord[2])
            except ValueError:
                reid = -1
                z = 0.0
            except IndexError:
                z = 0.0
                #self.hello("No Z coord value")
            if reid != -1 and reid not in resd or resd[reid] > z-ztop:
                resd[reid] = z - ztop

        resl = []
        for i in range(len(resd)):
            j = i+1
            if j in resd:
                resl.append(resd[j])

        return np.array(resl)
        #self.hello(str(resl))
        #self.hello(str(resd))
        #self.hello(str(self.ansl))

    def fthree(self):
        self.ansf = "/Users/fong/workspace/thesism/Gromacs/testrun/PSO/pso_conf/conf08401.pdb"
        if self.load_file(self.ansf):
            self.hello("Load the PSO answer structure.")
            self.ansl = self.residst()
            self.prota()
            self.ansx = self.degx
            self.ansy = self.degy
            self.ansz = self.degz
        else:
            self.hello("Can't find the PSO answer structure.")



        if len(sys.argv) >= 2:
            dname = sys.argv[1]
        else:
            self.hello("please give a directory name.")
            sys.exit(1)

        dtype = [('id', "a15"), ('rmsd', float), ('anglex', float), ('angley', float), ('anglez', float)]
        #np.set_printoptions(threshold=np.nan)
        value = []
        ls = os.listdir(dname)
        count = 0
        for f in ls:
            ext = f.split(".")[-1]  #split always return non-empth list
            if ext == "pdb":
                self.hello("find this file: {}".format(f))
                if self.load_file(os.path.join(dname, f)):
                    resl = self.residst()
                    rmsd = self.rmsd(self.ansl, resl)
                    self.prota()
                    value.append((f, rmsd, self.degx, self.degy, self.degz))
                count += 1
            #if count == 10:
                #break

        #print(value)
        L = np.array(value, dtype=dtype)
        L.sort(order="rmsd")
        #print(L)
        np.savetxt("np.out.gz", L, delimiter=" ", fmt="%s")

        J = np.loadtxt("np.out.gz", dtype=dtype)
        self.hello("Top 10 result!")
        print(J[:10])
        #print(np.sort(J, order="rmsd"))
        self.hello("The PSO answer angle is {}:{}:{}".format(self.ansx, self.ansy, self.ansz))

    def ffive(self):
        self.ansf = "/Users/fong/workspace/thesism/Gromacs/testrun/PSO/pso_conf/conf08401.pdb"
        if self.load_file(self.ansf):
            self.hello("Load the PSO answer structure.")
            self.ansl = self.residst()
            self.prota()
            self.ansx = self.degx
            self.ansy = self.degy
            self.ansz = self.degz
        else:
            self.hello("Can't find the PSO answer structure.")



        dtype = [('id', "a15"), ('rmsd', float), ('anglex', float), ('angley', float), ('anglez', float)]
        J = np.loadtxt("np.out.gz", dtype=dtype)
        dtype = [('id', "a15"), ('rmsd', float), ('anglex', float), ('angley', float), ('anglez', float)]
        dtype2 = dtype + [('rmsang', float)]
        J2 = []
        la = np.array([self.ansx, self.ansy, self.ansz])
        for i in J:
            lb = np.array([i[2], i[3], i[4]])
            rmsang = self.rmsd(la, lb)
            J2.append((i[0], i[1], i[2], i[3], i[4], rmsang))

        J2 = np.array(J2, dtype=dtype2)
        J2.sort(order="rmsang")
        print(J2[:10])


    def open_file(self, filename=""):
        if len(sys.argv) >= 2:
            pdbfile = sys.argv[1]
        else:
            self.hello("Please give a pdb file name")
            return 0

        if self.load_file(pdbfile):
            return 1
        else:
            self.hello("Please give a pdb file name")
            return 0

    def ffour(self):
        self.ansf = "/Users/fong/workspace/thesism/Gromacs/testrun/PSO/pso_conf/conf08401.pdb"
        if self.load_file(self.ansf):
            self.hello("Load the PSO answer structure.")
            self.ansl = self.residst()
            self.ansx = np.arange(1, len(self.ansl)+1)
            self.prota()
            self.hello("The PSO answer pdb structure angles: {}".format([self.degx, self.degy, self.degz]))
            self.ansx = self.degx
            self.ansy = self.degy
            self.ansz = self.degz
        else:
            self.hello("Can't find the PSO answer structure.")



        if len(sys.argv) >= 2:
            filename = sys.argv[1]
        else:
            self.hello("please give a pdb file name.")
            sys.exit(1)



        if self.load_file(filename):
            self.prota()
            self.hello("The target pdb structure angles: {}".format([self.degx, self.degy, self.degz]))
            resl = self.residst()
            resx = np.arange(1, len(resl)+1)

            fig = plt.figure()
            plt.plot(resx, self.ansl, 'b')
            plt.plot(resx, resl, 'k')
            plt.title("The minimum distance between the surface of each residues.")
            plt.ylabel(u"Unit Å")
            #plt.xticks(np.arange(min(self.ansx), max(self.ansx)+1, 1.0))
            #plt.xticks(self.ansx)
            plt.savefig("diagram_rxh.pdf", dpi=300)

    def fplot(self):
        groupb = ["06968", "05100", "03648", "01383"]
        #energb = ["-1013", "-1007", "-918", "-907"]
        energb = ["pso2", "pso3", "pso6", "pso7"]
        groupa = ["08401", "09682", "03314", "01680"]
        #energa = ["-1060", "-954", "-890", "-878"]
        energa = ["pso1", "pso5", "pso8", "pso9"]
        groupc = ["03492"]
        #energc = ["-959"]
        energc = ["pso4"]
        select_group = groupa
        select_energ = energa
        select_file = "group_a.pdf"
        select_name = "each residues"

        self.hello("Plot the group residues height diagram")
        fig = plt.figure()
        plt.rc('font', family='serif', size=20)
        ENGL = []
        minl = +7777777
        for i in select_group:
            self.hello("load the data of {}".format(i))
            L = np.loadtxt("conf"+i+".txt", dtype=[('resn', int), ('dest', float)])
            ENGL.append(L)
            X = L["resn"]
            Y = L["dest"]
            plt.plot(X, Y)
            if minl > Y.min():
                minl = Y.min()

        N = np.arange(-9, len(X)+18, 1)
        M = np.zeros(len(N)) + minl
        #plt.plot(N, M, 'k')
        #plt.text(-9, minl, "{}".format(minl), ha="right", va="center")
        plt.legend(select_energ, loc='upper left', ncol=4, frameon=False, labelspacing=0, columnspacing=0, handletextpad=0)
        plt.xlabel("Residue number")   #№
        #plt.ylabel(u"Minimum distance to surface of {} (Å)".format(select_name))
        plt.ylabel(u"Minimum distance to surface (Å)")
        plt.xlim(-7, len(X)+7)
        #plt.ylim(0, Y.max()*1.2)
        plt.ylim(0,37)
        #plt.xticks(np.arange(min(X), max(X)+1, 5.0))
        plt.tight_layout()
        plt.savefig(select_file, dpi=1200)
        self.hello("saved the diagram")

    def helping(self):
        """print the help message."""
        self.hello("="*77)
        self.hello("list of actions:")
        #ASSUMED: loacts is always correct.
        for i in self.loacts:
            act = getattr(self, i, self.helping)
            self.hello("{}: {}".format(i, act.__doc__))
        self.hello("="*77)

    def main(self):
        act = getattr(self, self.args.act, self.helping)
        if callable(act):
            self.hello("Run the action ({})".format(str(act)))
            act()
            








if __name__ == "__main__":
    m=measure()
    m.hello()
    #m.fplot()
    #m.fneo()
    #m.fthree()
    #m.ffour()
    #m.ffive()
    m.main()

