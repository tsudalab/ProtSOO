#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI jimmycfngai@gmail.com

from __future__ import print_function       # for python 2.7
from subprocess import call
import sys, shutil, os, argparse, json, re
import random, math, errno

import simpleMOVE

import numpy as np

from simpleLOG import logger
"""
        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M:%S',
        filename='debug.log',
        filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(levelname)8s : %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
logger = logging.getLogger(__name__)
# add the handler to the root logger
logger.addHandler(console)
#TODO: to change the logger level from cmd args.
#TODO: and disable debug.log file log
"""

from measure import measure


class measurem(measure):
    """override self.cmd object, I can reuse the methods in measure."""
    def __init__(self, pymol_cmd):
        """override this"""
        self.cmd = pymol_cmd


class simplePSO(simpleMOVE.simpleMOVE):

    def __init__(self):
        super(simplePSO, self).__init__()
        self.measure = measurem(self.cmd)
        logger.info("Initialized simplePSO object ")

    def initcmd(self):
        logger.info("Initialized command line arguments")
        self.parser = argparse.ArgumentParser(description="This is PSO testing")
        self.parser.add_argument("--outdir", help="set the output configuration directory.", required=False, default="pso_conf")
        self.parser.add_argument("--jsdbf", help="set the json db file path.", default="db.json")
        self.parser.add_argument("--proteinf", help="set the input protein pdb filename.", default="protein.pdb", metavar="protein.pdb")
        self.parser.add_argument("--surfacef", help="set the input surface pdb filename.", default="surface.pdb", metavar="surface.pdb")
        self.parser.add_argument("--n", help="set the number of birds.", default=200, required=False, type=int)
        self.parser.add_argument("--r", help="set the total iteration number.", default=10, required=False, type=int)
        self.parser.add_argument("--w", help="set the weight for updating velocity.", default=0.721, required=False, type=float)
        self.parser.add_argument("--c1", help="set the parameter C1.", default=1.193, required=False, type=float)
        self.parser.add_argument("--c2", help="set the parameter C2.", default=1.193, required=False, type=float)
        self.parser.add_argument("--minx", help="set the minimum range of X distance (Unit:Å)", default=1.0, required=False, type=float)
        self.parser.add_argument("--maxx", help="set the maximum range of X distance (Unit:Å)", default=5.5, required=False, type=float)
        self.parser.add_argument("--miny", help="set the minimum range of Y distance (Unit:Å)", default=1.0, required=False, type=float)
        self.parser.add_argument("--maxy", help="set the maximum range of Y distance (Unit:Å)", default=5.5, required=False, type=float)
        self.parser.add_argument("--minz", help="set the minimum range of Z distance (Unit:Å)", default=1.0, required=False, type=float)
        self.parser.add_argument("--maxz", help="set the maximum range of Z distance (Unit:Å)", default=5.5, required=False, type=float)
        self.parser.add_argument("--maxitr", help="set the maximum number of iteration.", default=0, required=False, type=int)
        self.parser.add_argument("--emdir", help="set the directory for EM, energy minimization using Gromacs.", default="EM", required=False, type=str)
        self.parser.add_argument("--keep-pdb", help="keep all of the pdf files generated during the searching.", default=False, action="store_true", required=False)
        self.initpg = self.parser.add_mutually_exclusive_group()
        self.initpg.add_argument("--resi", help="searching only the prefered residue surface", type=int, nargs="+")
        self.initpg.add_argument("--offset", help="set the initial searching orientation", nargs=6, type=float)
        self.initpg.add_argument("--init", help="use the exactly same molecule conformation as input", action="store_true", default=False)
        self.args = self.parser.parse_args()
        logger.debug("args: {}".format(self.args))


    def initvar(self):
        super(simplePSO, self).initvar()
        n = int(self.args.n)

        self.dimen = 6
        self.v = np.zeros((n, self.dimen))+1   # X, Y, Z rotation angle and Z distance
        self.x = np.zeros((n, self.dimen))
        self.gbests = 0.0                      # score is not vectors
        self.pbests = [0.0]*n
        self.gbestx = np.array([0.0]*self.dimen)
        self.pbestx = np.array([self.gbestx]*n)
        self.gbestx = np.zeros(self.dimen)
        self.pbestx = np.zeros((n, self.dimen))
        self.gbestn = 99
        self.coff   = 20                        # cut-off for checking dimension size
        self.gbestf = ""
        #ASSUMED: x range a and b can not be all zero.
        #ASSUMED: range a is always less than b
        self.xrangea = [0.0, 0.0, 0.0, self.args.minx, self.args.miny, self.args.minz]
        self.xrangeb = [360.0, 360.0, 360.0, self.args.maxx, self.args.maxy, self.args.maxz]
        #TODO:this should let user input
        self.maxv = [036.0, 036.0, 036.0, (self.args.maxx-self.args.minx)*0.1, \
        (self.args.maxy-self.args.miny)*0.1, (self.args.maxz-self.args.minz)*0.1]

        logger.debug("initialized PSO variables: {}".format(self.__dict__))
        self.tid = 1        # pdb file id, incremental

    def reslst(self):
        """compute the list of residues distance"""
        #ASSUMED all residue id are continuous number
        lst = self.measure.residst()
        lst = np.array(lst)
        lst_min = lst.min()
        lst = lst - lst.min()
        return lst

    def comp_reslst(self, resl):
        """Compare the residue list with-in adsorption site.
        It will obtain the list of residues height internally."""
        one_in = all_in = False
        in_count = 0
        reslst = self.reslst()
        logger.debug("residue height list: {}".format(str(reslst)))
        for i in resl:
            if i >= 0 and i < len(reslst):
                i_height = reslst[i]
            if i_height <= 5:       #TODO: hard-code number
                one_in = True
                in_count += 1
        if in_count == len(resl):
            all_in = True

        return one_in


    def preprocessing(self):
        """Check the inputs, esidue list existing.
        Compute the radius of adsorption side input list.
        This function should be called after load the pdb molecules files.
        """
        n = int(self.args.n)

        if self.args.resi:
            resr = self.face2surface(self.args.resi)
            self.xrangea[0] = self.xrangea[1] = -1*resr
            self.xrangeb[0] = self.xrangeb[1] = +1*resr
        if self.args.offset:
            self.initorient(self.args.offset)

        logger.debug("RAndomly assign the initial birds position and velocity.")
        self.v = self.v * np.random.rand(n, self.dimen)
        self.x = self.x + np.random.rand(n, self.dimen)
        #print ("preprocessing:",self.x)
        for i in range(6):
            self.x[:,i] *= (self.xrangeb[i] - self.xrangea[i])
            self.x[:,i] += self.xrangea[i]

        logger.debug("Measuring the dimension of protein and surface molecules")
        r = self.measure.mxlength("protein")    #sqroot always return positive
        logger.debug("Protein maximum length: {}".format(r))
        self.measure.xyzedges("surface")
        rx = np.absolute(self.measure.xedgep - self.measure.xedgen)
        ry = np.absolute(self.measure.yedgep - self.measure.yedgen)
        rz = np.absolute(self.measure.zedgep - self.measure.zedgen)
        logger.debug("Surface dimension: {}, {}, {}".format(rx, ry, rz))
        #c=2nm = 20
        #rx > r+c
        #ry > r+c
        #simbox z > rz + r+c+1
        text = ""
        with open("set-up.sh", "r") as setupf:
            text = setupf.read()
        setupf.close()

        simsiz = 0.0
        m = re.search('sysboxs="(\d+) (\d+) (\d+)"', text)
        if m:
            simsix = float(m.group(1)) * 10.0
            simsiy = float(m.group(1)) * 10.0
            simsiz = float(m.group(3)) * 10.0
        logger.debug("simulation box height: {}".format(simsiz))
        c = self.coff
        if not rx >= r+c or not ry >= r+c or not simsiz >= rz+r+c+1 or not simsix >= r+c or not simsix >= r+c:
            logger.error("System dimension mismatched!!")
            logger.error("Please check the surface dimension or the simulation box size is larger enough")
            sys.exit("System dimension mismatched")

        if self.args.resi or self.args.offset or self.args.init:
            self.save_init()

    def opendb(self):
        """
        json db schema: {
            "N": int,
            "R": int,
            "bests": float,
            "bestb": int,
            "besti": int,
            "bestf": str,
            "birds": [ bird ]
        }

        bird: {
            "iteration": int,
            "bird": int,
            "energy": float,
            "postition": [float, float, float, float],
            "velocity": [float, float, float, float],
            "gbest": bool,
            "fpath": str
        }
        """
        try:
            logger.debug("Trying to load json db file.")
            with open(self.args.jsdbf, "r") as openfile:
                self.jsdb = json.load(openfile)
                logger.info("Loaded json db file ")
        except IOError:
            self.jsdb = json.loads("{}")        # create an empty json
            self.jsdb["N"] = self.args.n
            self.jsdb["R"] = self.args.r
            self.jsdb["bests"] = 0.0
            self.jsdb["bestb"] = 0
            self.jsdb["besti"] = 0
            self.jsdb["birds"] = []
            logger.info("Could not find previous json db, initialized a new one ")

    def savejdb(self):
        with open(self.args.jsdbf, "w") as openfile:
            json.dump(self.jsdb, openfile)
        logger.debug("Saved the jsdb file: {}.".format(self.args.jsdbf))


    def myPSO(self):
        if len(self.pbestx) != self.args.n or len(self.pbests) != self.args.n or \
            len(self.x) != self.args.n or len(self.v) != self.args.n:
            logger.warning("intial X and V does not initialized correctly !!")
            return -1
        for i in self.x + self.v:
            if len(i) != self.dimen:
                logger.warning("X or V does not have enought number of elements {}".format(self.dimen))
        logger.debug("checking completed! Variables オールグリーン！")



        c = self.args.r
        if self.args.maxitr:
            c = self.args.maxitr
        total = 0
        logger.info("{} birds have been initialized, PSO searching start!".format(self.args.n))
        while c > 0:
            logger.info("[===PSO===] iteration number: {}, r={} is remaining".format(total, c))
            for i in range(self.args.n):
                # scoringone, generate the pdb system and then calculate the energy
                x = self.x[i]
                x=[188,190,200,1,2,3]
                #print (x.shape())
                print ("myPSO:",x)
                v = self.v[i]
                for idx, step in enumerate(x[:3]):
                    self.protate(self.xyzcv[idx], step)
                for idx, step in enumerate(x[3:]):
                    self.ptransl(self.xyzcv[idx], step)


                #Drop invalid surface conformation if user input the list of adsorption site
                if self.args.resi and not self.comp_reslst(self.args.resi):
                    logger.debug(
                    "[===PSO=== count({})] Drop this conformation because not included in users input adsorption site".format(total))
                    logger.debug("The position is: {}".format(x))
                    fpath = "N/A"
                    s = 0.0


                else:
                    fpath = self.savers()
                    #print ("fpath:",fpath)
                    logger.debug("[===PSO=== count({})] The bird {:02} arrived, saved as pdb file: {}".format(total, i, fpath))
                    logger.debug("The position is: {}".format(x))
                    print (format(x))

                    jsons = self.scoringone(fpath, mdir=self.args.emdir)
                    #print ("jsons:",jsons)
                    s = jsons["energy"] * -1
                    print ("energy:",s)
                #logger.debug("[===PSO=== count({})] With energy value: {}".format(total, s))
                #jsona = {}
                #jsona["iteration"] = total
                #jsona["bird"] = i
                #jsona["energy"] = s
                #jsona["position"] = x.tolist()
                #jsona["position"] = x
                #jsona["velocity"] = v.tolist()
                #jsona["gbest"] = False
                #jsona["fpath"] = fpath
                if s > self.pbests[i]:
                    self.pbests[i] = s
                    self.pbestx[i] = list(self.x[i])
                if s > self.gbests:
                    self.gbests = s
                    self.gbestx = list(self.x[i])
                    self.gbestn = i
                    self.gbestf = fpath
                    jsona["gbest"] = True
                    if not self.args.maxitr:
                        c = self.args.r         #reset the c if a better solution was found
                    logger.debug("Updated global best: bird #{:03}@score:{}".format(self.gbestn, self.gbests))

                    shutil.copyfile(fpath, "gbest.pdb")     #TODO gbest.pdb file name is hard-code!
                    logger.debug("Copied the conformation to gbest.pdb")
                if not self.args.keep_pdb and os.path.isfile(fpath):
                    os.remove(fpath)
                    logger.debug("Removed the pdb file to save space.")
                try:
                    self.jsdb["birds"].append(jsona)
                except KeyError:
                    self.jsdb["birds"] = [jsona]
                logger.debug("Updated gbest and pbest scoring and position.")
                # update velocity v and position x
                self.v[i] = self.args.w*self.v[i] + self.args.c1*random.random()*(self.pbestx[i] - self.x[i]) +\
                        self.args.c2*random.random()*(self.gbestx - self.x[i])
                if len(self.maxv) >= len(self.v[i]):
                    for k, v in enumerate(self.maxv):
                        if np.abs(self.v[i][k]) > v:
                            velocity = self.v[i][k]
                            velocity = (velocity / np.abs(velocity)) * v
                            self.v[i][k] = velocity
                            logger.debug("reset the velocity index:{} value:{}".format(k, self.v[i][k]))
                else:
                    logger.warning("xrangea/xrangeb/vmax array does not have enought elements!")
                    logger.warning("Cant update position/velocity x/v correctly!")
                self.x[i] += self.v[i]
                logger.debug("position x before normalized: {}".format(self.x[i]))
                if len(self.x[i]) > len(self.xrangea) or len(self.x[i]) > len(self.xrangeb):
                    logger.warning("xrangea or xrangeb array does not have enought elements!")
                    logger.warning("Cant update position/velocity x/v correctly!")
                else:
                    #MOD operations, let x within the range.
                    px = list(self.x[i])
                    self.x[i] = [ (px[x]-self.xrangea[x])%(self.xrangeb[x] - self.xrangea[x]) + self.xrangea[x] for x in range(len(px)) ]
                    # special z translation
                    # Since numpy does not allow to copy the array without the last element, z translation
                    # Here we replace it by the backup px.
                    ftz = px[-1]
                    bnd = 0
                    if px[-1] < self.xrangea[-1]:
                        bnd = self.xrangea[-1]
                    if px[-1] > self.xrangeb[-1]:
                        bnd = self.xrangeb[-1]
                    if bnd:
                        #This is another boundary condition:
                        #ftz = bnd - math.fmod((px[-1] - bnd), (self.xrangeb[-1] - self.xrangea[-1]))
                        ftz = bnd   #this stop at the boundary.
                    self.x[i][-1] = ftz
                logger.debug("updated velocity({}): {}".format(i, self.v[i]))
                logger.debug("updated position({}): {}".format(i, self.x[i]))
            c-=1
            total+=1

        if self.gbests == 0:
            logger.info("No favorable orientation is found")
            sys.exit(1)
        logger.info("Finally, PSO stop after {} number of iterations ".format(total))
        logger.info("Found the best scoring result ")
        logger.info("Bird ID: {:03}".format(self.gbestn))
        #ASSUMED list length is corrent, will not out of range
        logger.info("Rotation (deg): x={0} y={1} z={2}".format(*self.gbestx))
        logger.info("Translation (Ang): x={3} y={4} z={5}".format(*self.gbestx))
        logger.info("Score (kJ/mol): {}".format(self.gbests*-1))
        logger.info("Output files:")
        logger.info("Search history file: {}".format(self.args.jsdbf))
        logger.info("Final gbest structure: {}".format("gbest.pdb"))

    def savers(self):
        """save a pdb file and restore the protein to original position from backup"""
        fpath = self.save_file(self.tid)
        self.tid+=1
        logger.debug("self.tid increased({})".format(self.tid))
        self.cmd.create("protein", "protein_bak")
        logger.debug("restored the protein to original position.")

        return fpath


    def save_file(self, state=0):
        """save pse pymol session file if state equal to zero,
        save pdb stand format if state not zero"""
        if not os.path.isdir(self.args.outdir):
            os.makedirs(self.args.outdir)       # ASSUMED outdir is a correct path, NO .. parent path

        fpath = ""
        if state == 0:
            logger.debug("PASSED, don't save the session file. it is useless.")
        else:
            self.cmd.create("system", "protein")
            self.cmd.order("system", location="top")
            fpath = os.path.join(self.args.outdir, "conf{:05}.pdb".format(self.tid))
            self.cmd.save(fpath, "system surface")
            logger.debug("saved pdb file: {}".format(fpath))

        return fpath

    def save_init(self):
        """save the system initial conformation for debug propose"""
        if not os.path.isdir(self.args.outdir):
            os.makedirs(self.args.outdir)       # ASSUMED outdir is a correct path, NO .. parent path
        fpath = os.path.join(self.args.outdir, "conf_init.pdb")
        self.cmd.create("system", "protein_bak")
        self.cmd.order("system", location="top")
        self.cmd.save(fpath, "system surface")
        logger.debug("saved the initial conformation as pdb file (conf_init.pdb)")


    def main(self):
        logger.debug("Protein Surface adsorption PSO search algorithm")
        logger.debug("はじまるよ～♪ ")

        self.load_file()
        self.load_post(self.args.init)
        #ASSUMED all files has been alread loaded in the following function without checking
        self.preprocessing()
        self.myPSO()
        self.savejdb()



if __name__ == "__main__":
    simplePSO().main()
