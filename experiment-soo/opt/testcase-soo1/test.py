
from __future__ import print_function       # for python 2.7
from subprocess import call
import sys, shutil, os, argparse, json, re
import random, math, errno

import simpleMOVE

import numpy as np

from simpleLOG import logger
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

        self.parser.add_argument("--maxitr", help="set the maximum number of iteration.", default=0, required=False, type=int)
        self.parser.add_argument("--emdir", help="set the directory for EM, energy minimization using Gromacs.", default="EM", required=False, type=str)
        self.parser.add_argument("--keep-pdb", help="keep all of the pdf files generated during the searching.", default=False, action="store_true", required=False)
        self.initpg = self.parser.add_mutually_exclusive_group()
        self.initpg.add_argument("--resi", help="searching only the prefered residue surface", type=int, nargs="+")
        self.initpg.add_argument("--offset", help="set the initial searching orientation", nargs=6, type=float)
        self.initpg.add_argument("--init", help="use the exactly same molecule conformation as input", action="store_true", default=False)
        self.args = self.parser.parse_args()
        logger.debug("args: {}".format(self.args))




    def myPSO(self):

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

    def main(self):
        #ilogger.debug("Protein Surface adsorption PSO search algorithm")

        self.load_file()
        #self.load_post(self.args.init)
        #ASSUMED all files has been alread loaded in the following function without checking
        #self.preprocessing()
    #    self.myPSO()
        #self.savejdb()





if __name__ == "__main__":
    simplePSO().main()
