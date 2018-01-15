#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI lolicon.jimmy@gmail.com

import sys, os, argparse, math, logging, json, re
import numpy as np

import pdb

from simpleLOG import logger

class simpleANS:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="This is simpleANS.py")
        self.parser.add_argument("--pdb", help="set the input protein pdb filename.", required=False, default="gbest.pdb")
        self.parser.add_argument("--jdb", help="the filename of json db file.", required=False, default="db.json")
        self.args = self.parser.parse_args()

        # pymol launching
        import pymol
        pymol.pymol_argv = ['pymol','-qc'] 
        pymol.finish_launching()
        self.cmd = pymol.cmd
        self.cmd.set("retain_order", 1)
        self.cmd.set("pdb_retain_ids", 1)

        logger.info("Starting to analyze the lowest energy orientation and search trajectory           ")



    def load_file(self, filenm):
        self.filenm = ""
        if not os.path.isfile(filenm):
            logger.error("No such ({}) pdb file!".format(filenm))
            return 0

        self.cmd.delete("*")    #Clean everything first
        try:
            self.cmd.load(filenm, "system")
            self.filenm = filenm
            self.cmd.select("surface", "resn PTS or resn SUB")
            self.cmd.select("protein", "not resn SOL and not resn PTS and not resn SUB")
        except:
            logger.error("Loading file ERROR!!!")
            return 0

        return 1

    def xyzedges(self, s="surface"):
        """measure the X/Y/Z edges and length."""
        #PyMol already implemented get_extent function to get the X Y Z min/max values. see: http://www.pymolwiki.org/index.php/Get_extent
        try:
            atoms = self.cmd.get_model(s)
        except:
            logger.error("Can not find model({})".format(s))
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

        self.centerp = [(j-i)/2 +i for i,j in [(self.xedgen, self.xedgep), (self.yedgen, self.yedgep), (self.zedgen, self.zedgep)]]


    def residst(self):
        self.xyzedges("surface")
        ztop = self.zedgep      #Ugly!

        resd = {}
        try:
            atoms = self.cmd.get_model("protein")
        except:
            logger.error("Could not find protein model.")
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
                logger.error("No Z coord value")
            if reid != -1 and reid not in resd or resd[reid] > z-ztop:
                resd[reid] = z - ztop

        resl = []
        for i in range(len(resd)):
            j = i+1
            if j in resd:
                resl.append(resd[j])

        return np.array(resl)

    def resh(self):
        """
        plot the residues minimum distance graph.
        """

        filenm = self.args.pdb

        if self.load_file(filenm):
            resl = self.residst()
            resx = np.arange(1, len(resl)+1)

            L = np.column_stack((resx, resl))
            np.savetxt("gbest.txt", L, fmt="%03d, %02.9f", header="residue, distance (angstrom)")
            logger.info("Final gbest residue min-distance profile:     gbest.txt ")
            L = np.loadtxt("gbest.txt", dtype=[('resn', int), ('dest', float)], delimiter=",")
            L.sort(order=['dest'])
            np.savetxt("gbest_sorted.txt", L, fmt="%03d, %02.9f", header="residue, distance (angstrom)")
            logger.info("Sorted by the distance of each residue:       gbest_sorted.txt ")


    def loaddb(self, name=""):
        try:
            with open(name, "r") as openfile:
                self.jsdb = json.load(openfile)
        except (IOError, TypeError):
            logger.error("Can't open the file: {}".format(name))
            return 1

        return 0

    def energyP(self, target=-1):
        """
        Just plot the gbest's energy values and location as function of iterations

        Its design can be imporved, it mixed up tow tasks in this function.
        Locate the target bird's values and find the gbest bird.
        I think these should be two separated tasks.
        The current version is try to save time, find the gbest bird with its values on-the-fly
        """
        try:
            target = int(target)
        except ValueError:
            logger.error("Please gives the target ID, should be integer!")
            return 1

        if self.loaddb(self.args.jdb) != 0:
            logger.error("Please check the db file existed or not!")
            return 1

        try:
            birds = self.jsdb["birds"]
        except KeyError:
            logger.error("jsdb does not have birds!")
            return 1

        if not isinstance(birds, list):
            logger.error("there are no birds in it!")
            return 1

        lasti = 0
        I = []
        E = []
        X = []
        Y = []
        Z = []
        A = []
        B = []
        C = []
        P = []
        EE = []
        maxe = -999999999
        mine = +999999999
        Pnp = np.array([])

        laste = 0
        lastxr = lastyr = lastzr = 0
        lastxt = lastyt = lastzt = 0

        for bird in birds:
            try:
                n = int(bird["bird"])
                i = int(bird["iteration"])
                e = float(bird["energy"]) * -1
                g = bool(bird["gbest"])
                v = bird["velocity"]
                p = bird["position"]
                xr = float(p[0])
                yr = float(p[1])
                zr = float(p[2])
                xt = float(p[3])
                yt = float(p[4])
                zt = float(p[5])
            except (TypeError, IndexError, ValueError, KeyError):
                logger.error("There are alien in the birds!")
                return 1

            if i == 0:
                P.append(e)
            if i == 1:
                Pnp = np.array(P)
            if P[n] > e:
                P[n] = e
            if mine > e:
                mine = e
            if maxe < e:
                maxe = e
            if (target == -1 and g) or (target != -1 and target == n):
                laste = e           #will case ERROR if this statement never runs
                lastxr = xr
                lastyr = yr
                lastzr = zr
                lastxt = xt
                lastyt = yt
                lastzt = zt
            if (target != -1 and target == n) or (target == -1 and i > lasti):
                if target == -1:    #fixed incorrect iteration number for non-gbest target
                    ino = lasti     #for gbest, the data is added after the iteration
                else:
                    ino = i         #for non-gbest, the data is added on-the-fly
                I.append(ino)
                E.append(laste)
                A.append(lastxr)
                B.append(lastyr)
                C.append(lastzr)
                X.append(lastxt)
                Y.append(lastyt)
                Z.append(lastzt)
                lasti = i
                maxe = -999999999
                mine = +999999999
                if i == 0:          #TODO: can get EE for non-gbest target
                    EE.append(0)    #because the data is added on-the-fly
                elif len(Pnp) > 0:               
                    #it's OK if added after the iteration
                    EE.append(Pnp.std())


        #gbest must be added after the iteration
        #because there are more than one gbest records in one iteration
        if target == -1:
            I.append(lasti)
            E.append(laste)
            A.append(lastxr)
            B.append(lastyr)
            C.append(lastzr)
            X.append(lastxt)
            Y.append(lastyt)
            Z.append(lastzt)
            lasti = i
            maxe = -999999999
            mine = +999999999
            if len(Pnp) > 0:
                EE.append(Pnp.std())

        I = np.array(I)
        E = np.array(E)
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        A = np.array(A)
        B = np.array(B)
        C = np.array(C)
        EE = np.array(EE)


        #save all the numbers into txt
        EN = np.arange(1, len(E)+1)
        ES = np.column_stack((EN, E))
        np.savetxt("gbest_energy.txt", ES, fmt="%03d, %2.9e", header="iteration, ProtPOS score (kJ/mol)")
        logger.info("Gbest energy evolution:                       gbest_energy.txt")
        EV = np.column_stack((EN, A, B, C, X, Y, Z))
        np.savetxt("gbest_vector.txt", EV, fmt="%03d, %2.9e, %2.9e, %2.9e, %2.9e, %2.9e, %2.9e", \
            header="iteration, (orientation:) rotation-X, rotation-Y, rotation-Z, translation-X, translation-Y, translation-Z")
        logger.info("Gbest orientation evolution:                  gbest_vector.txt")






if __name__ == "__main__":
#two main function:
    ans = simpleANS()
#   1. resh from measure, to generate the minimum distance profile
    ans.resh()   
#   2. gbest energy from plotp, to generate the gbest energy and position 
    ans.energyP()
