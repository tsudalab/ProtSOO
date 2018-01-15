#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI lolicon.jimmy@gmail.com

from __future__ import print_function       # for python 2.7
from subprocess import call
import sys, os, argparse, json, re
import random, math, errno

import numpy as np

from simpleLOG import logger
"""
logging.basicConfig(level=logging.DEBUG,
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



class simpleMOVE(object):

    def initcmd(self):
        # initial the cmd args
        self.parser = argparse.ArgumentParser(description="This is MC testing")
        self.parser.add_argument("--total", help="total number of sample", default=10000, type=int)
        self.parser.add_argument("--steps", help="set the rotation sampling steps.", default=100, type=int)
        self.parser.add_argument("--stepd", help="set the moving down steps number.", default=1, type=int)
        group = self.parser.add_mutually_exclusive_group()
        group.add_argument("-k", help="the number unit of sample is 10^3 (1K)", action="store_true")
        group.add_argument("-m", help="the number unit of sample is 10^6 (1M)", action="store_true")
        self.parser.add_argument("--minlv", help="set the minimum energy level",  type=int)
        self.parser.add_argument("--maxlv", help="set the maximum energy level",  type=int)
        self.parser.add_argument("--act", help="run something main/etc...",  default="", \
                choices=['sampling', 'testing', 'scoring', 'simulate'])
        self.parser.add_argument("--log", help="set the logging level to DEBUG/INFO/WARNING", default="WARNING", \
                choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], type=str)
        self.parser.add_argument("--proteinf", help="set the input protein pdb filename.", default="protein.pdb", metavar="protein.pdb")
        self.parser.add_argument("--surfacef", help="set the input surface pdb filename.", default="surface.pdb", metavar="surface.pdb")
        self.parser.add_argument("--topolf", help="set the input topol top filename.", default="topol.top", metavar="topol.top")
        self.parser.add_argument("--systemf", help="set the output system pdb filename.", default="system.pdb")
        self.parser.add_argument("--gromacs", help="set the executable path of gromacs command.", default="gromacs")
        self.parser.add_argument("--init-height", help="set the initial height, unit is 1Å", type=float, metavar="N")
        self.parser.add_argument("--init-xr", help="set the initial X axis rotation angle, unit is degrees", type=float, default=0.0)
        self.parser.add_argument("--init-yr", help="set the initial Y axis rotation angle, unit is degrees", type=float, default=0.0)
        self.parser.add_argument("--init-xt", help="set the initial X axis translation distance, unit is 1Å", type=float, default=0.0)
        self.parser.add_argument("--init-yt", help="set the initial Y axis translation distance, unit is 1Å", type=float, default=0.0)
        self.parser.add_argument("--init-zr", help="set the initial Z axis rotation angle, unit is degree", type=float, default=0.0)
        self.parser.add_argument("--random-zdist-range", help="set the random z distance range.", type=int, default=0)
        self.parser.add_argument("--no-plot", help="skip the plot progress.", action="store_true")
        self.parser.add_argument("--no-plot-re", help="skip the residues plot progress.", action="store_true")
        self.parser.add_argument("--jsdbf", help="set the json db file path.", default="db.json")
        self.args = self.parser.parse_args()

    def initpmo(self):
        """
        initial the PyMOL environment
        """
        # pymol environment TODO: only for MacOS, e.g. pymol is installed from MacPort
        self.moddir='/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/'
        sys.path.insert(0, self.moddir)
        os.environ['PYMOL_PATH'] = os.path.join(self.moddir, 'pymol/pymol_path')

        # pymol launching
        try:
            import pymol
        except:
            #logger.error("PyMOL environment was not installed correctly!!")
            sys.exit(1)

        #pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
        pymol.pymol_argv = ['pymol','-qc']
        pymol.finish_launching()
        self.cmd = pymol.cmd
        self.cmd.set("retain_order", 1)
        self.cmd.set("pdb_retain_ids", 1)
        #logger.info("PyMOL environment initialized   ")

    def __init__(self):
        self.initcmd()
        self.initpmo()
        self.initvar()
        self.opendb()
        #logger.info("Initialized simpleMOVE objects")

    def initvar(self):
        # initial the variables
        self.systemfn = "system"
        self.topol = []
        self.svcount = 0
        self.atom_cutoff = 7                # Stop downward when closer
        self.miniZdistR = 0                 # random distance range is 1 A
        self.bash = "/bin/bash"             # bash path
        self.runsh = "score.sh"             # score.sh script
        self.runshp = self.runsh            # default is same as runs.h
        self.energyx = "energy.xvg"         # energy.xvg file name
        self.mdir = "EM"                    # EM directory
        self.ctlvs = 10                     # contour map levels number
        self.slogn = "score.log"            # output log of score.sh
        self.xyzcv = ["x", "y", "z"]
        self.xyzid = {"x":0, "y":1, "z":2}
        self.testname = "yes, I am simpleMOVE."

        logger.debug("searching {}...".format(self.runsh))
        for sysp in sys.path:
            for root, dirs, files in os.walk(sysp):
                for f in files:
                    if f == self.runsh:
                        self.runshp = os.path.abspath(os.path.join(root, f))
                        logger.debug("{} was found!".format(self.runshp))
                        # it will return the last one found


        #logger.debug("initialized simpleMOVE variables.")

    def opendb(self):
        try:
            with open(self.args.jsdbf, "r") as openfile:
                self.jsdb = json.load(openfile)
                #logger.info("Loaded json db file.")
        except IOError:
            self.jsdb = json.loads("{}")        # create an empty json
            self.jsdb["mindst"] = 0
            self.jsdb["confn"] = 0
            self.jsdb["confd"] = ""
            self.jsdb["confs"] = []
            #logger.info("Can not find previous json db, initial a new one.")
        except AttributeError:
            logger.warning("Seems there are no jsdbf setting")
        self.tmdb = {}

        #TODO: print out all system conf variables




    def savefile(self, diry=None, pdb=False):
        self.jsdb["confn"] = 0      #TODO may need to clear the db?
        self.jsdb["confd"] = ""
        self.jsdb["confs"] = []

        if not diry:
            diry = "."
        if not os.path.isdir(diry):
            os.mkdir(diry)
        self.cmd.save("{}/system{:05}.pse".format(diry, self.svcount))
        logger.info("Saved state file file.")
        if pdb:
            try:
                os.mkdir("{}/conf{:05}".format(diry, self.svcount))
                logger.debug("Creating dir: conf{:05}".format(self.svcount))
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir("{}/conf{:05}".format(diry, self.svcount)):
                    logger.debug("Already exist, pass.")
                    pass
                else:
                    logger.warning("Can not create directory, skip...")
                    return
            finally:
                for i in range(self.cmd.count_states("mov")):
                    self.cmd.create("system", "mov", i+1, 1)
                    self.cmd.order("system", location="top")
                    self.cmd.save("{}/conf{:05}/system{:05}.pdb".format(diry, self.svcount, i), "system surface")
                    logger.debug("Saved system{:05}.pdb (state: {}) under directory: conf{:05}.".format(i, i+1, self.svcount))

                    pdbfs = {}
                    pdbfs["xi"] = i % self.args.steps
                    pdbfs["yi"] = i / self.args.steps
                    pdbfs["id"] = "system{:05}".format(i)
                    self.jsdb["confs"].append(pdbfs)
                    logger.debug("append ({}) confs, xi: {} yi: {}".format(pdbfs["id"], pdbfs["xi"], pdbfs["yi"]))

            self.jsdb["confn"] = self.cmd.count_states("mov")
            self.jsdb["steps"] = self.args.steps
            self.jsdb["confd"] = "conf{:05}".format(self.svcount)
            logger.debug("Saved tmp data, {} number of confs under {} directory.".format(self.jsdb["confn"], self.jsdb["confd"]))
        self.svcount += 1


    def gtha(self, atoms):
        d = {}
        for at in atoms.atom:
            logger.info("ATOM DEFINITION: "+at.chain+" "\
                                     +at.resn+" "\
                                     +str(at.resi)+" "\
                                     +str(at.name)+" "\
                                     +str(at.index)+" "\
                                     +str(at.b)+" "\
                                     +str(at.coord[0])+" "\
                                     +str(at.coord[1])+" "\
                                     +str(at.coord[2]))
            d[at.resn] = at.resi
        return d

    def testing(self):
        self.load_file()
        self.load_top()



    def savejdb(self):
        #ASSUMED: it should have jsdbf argument.
        with open(self.args.jsdbf, "w") as openfile:
            json.dump(self.jsdb, openfile)
        logger.info("Saved the jsdb file: {}.".format(self.args.jsdbf))

    def load_file(self):
        logger.debug("loading "+self.args.proteinf)
        self.cmd.load(self.args.proteinf, "protein")
        logger.debug("loading "+self.args.surfacef)
        self.cmd.load(self.args.surfacef, "surface")
        logger.info("loaded protein and surface pdb files.")
        #TODO: load will raise error


    def load_post(self, skip=False):
        #ASSUMED: it is called after load_file
        if not skip:
            self.cmd.center("surface")
            surfacec = self.cmd.get_position()
            logger.debug("Center of the surface: " + ", ".join(format(x, "f") for x in surfacec))
            self.cmd.center("protein")
            proteinc = self.cmd.get_position()
            logger.debug("Center of the protein: " + ", ".join(format(x, "f") for x in proteinc))
            delta = [i-j for i, j in zip(proteinc, surfacec)]
            logger.debug("move delta: "+", ".join(format(x, "f") for x in delta))
            self.cmd.translate(delta, "surface", camera=0)
            logger.debug("Translate Z axis of protein-surface.")
            miniZ = self.mindst(state=0, z=True)
            self.cmd.translate([0, 0, miniZ-0], "surface", camera=0)    #let the distance to be zero, translate later
            logger.debug("The minimal Z distance now is : {}".format(self.mindst(state=0, z=True)))

        #logger.info("The initial structure has been created ")
        self.cmd.create("protein_bak", "protein")
        logger.debug("backup the original initial protein position")

    def face2surface(self, resi):
        """
        rotate the protein let the residue number resi face to the surface.
        input argument resi should be a list which contains at least one element.
        """
        #ASSUMED the list of residue (resi) is valid
        self.cmd.center("protein")
        proteinc = np.array(self.cmd.get_position())


        xi = [1, 0, 0]
        yi = [0, 1, 0]

        #TODO: check resi in range
        #TODO: should project the vector to x-z and y-z plane
        self.cmd.center("protein and resi {}".format("+".join(str(x) for x in resi)))
        residuec = np.array(self.cmd.get_position())
        rvt = residuec - proteinc

        cosx = np.dot(rvt, xi)/(np.linalg.norm(rvt)*np.linalg.norm(xi))
        cosy = np.dot(rvt, yi)/(np.linalg.norm(rvt)*np.linalg.norm(yi))
        degx = np.degrees(np.arccos(cosx))
        degy = np.degrees(np.arccos(cosy))

        logger.info("Rotate the residue group to face to surface.")
        logger.info("degx: {}, degy: {}".format(degx, degy))
        self.protate("y", -90-degx)
        self.protate("x", -90-degy)

        self.cmd.create("protein_bak", "protein")
        logger.debug("backup the initial protein position after face to surface")

        #fpath = os.path.join(self.args.outdir, "bak_init.pdb")
        #self.cmd.save(fpath, "protein surface")
        #logger.debug("backup the initial orientation pdb file: {}".format(fpath))


        #Compute the radius
        #cosine = adjacent over hypotenuse
        hypoten = 0
        for i in resi:
            self.cmd.center("protein and resi {}".format(str(i)))
            ric = np.array(self.cmd.get_position())
            length = np.linalg.norm(ric - proteinc)
            if length > hypoten:
                hypoten = length
        adjacent = np.linalg.norm(rvt)
        cosTheta = adjacent / (hypoten + 5)        #hard code! 5 is the adsorption site cut-off
        #It should even works if only one residue in the list
        Theta = np.degrees(np.arccos(cosTheta))
        logger.debug("Adjacent({}) over hypoten({}) +5".format(adjacent, hypoten))
        logger.debug("Computed cosTheta: {}, and Theta: {}".format(cosTheta, Theta))

        return Theta

    def initorient(self, initp):
        """
        set the initial orientation
        """
        xyz = ["x", "y", "z"]
        for i in range(3):
            self.protate(xyz[i], initp[i])
        #preform rotation and translation separately, try to avoid any accuracy problem
        for i in range(3):
            self.ptransl(xyz[i], initp[i+3])

        self.cmd.create("protein_bak", "protein")
        logger.debug("backup the initial protein position after initial the starting orientation")


    def load_top(self):
        logger.info("loading dihedral angles from gromacs top file.")
        f = open(self.args.topolf)
        logger.info("opened file: "+self.args.topolf)
        m = "no"                #mode, true when goes into dihedrals section
        logger.info("reading dihedrals ...")
        anglen = 0
        for l in f.readlines():
            if len(l) > 0 and l[0] == ";":
                continue        #skip all comments
            if l.strip() == "[ dihedrals ]":
                m = "dh"        #this is dihedrals section
            if m == "dh" and l.strip() == "":
                m = "no"        #ending by empty line
            if m == "dh":
                r = l.split()   #split heading and ending space chs
                logger.debug("line: "+", ".join(r))
                if len(r) >= 4 and r[4] == "1":     #we need funct 1 only.
                    #filter out all backbone atoms
                    is_bkbon = 0
                    for i in r[:4]:
                        logger.debug("get_model protein and id "+i)
                        atoms = self.cmd.get_model("protein and id "+i)
                        for idx, a in enumerate(atoms.atom):
                            logger.debug("({})~ [ID: {}] chain {} res {} ({}) atom {} ({})".format(idx, i, a.chain, a.resn, a.resi, a.name, a.index))
                    if not is_bkbon:
                        self.topol.append(r[:4])    #we only need non-backbone
                        anglen += 1
                    logger.debug("record was just added: "+" ".join(map(str, self.topol[-1])))
        logger.info("finished reading top !!")
        random.shuffle(self.topol)      #TODO: randomly insert on the fly.

    def mindst(self, state=-1, z=True, surfc=None):
        '''Find the minimal distance bewteen the surface and the protein
        if the mov model not yet exist, please give state=0 to choose the protein model
        z=False may not work in current version'''
        pdown = float("+inf")
        adown = [0.0, 0.0, 0.0]
        rdown = "" #resi of atom
        ndown = "" #name of atom
        if surfc:
            surf  = surfc
        else:
            surf  = self.surfaceh()
        stop  = surf["stop"]
        atop  = surf["atop"]
        rtop  = surf["rtop"]
        ntop  = surf["ntop"]
        #loop all atoms to find the top of surface and the bottom of protein

        try:
            int(state)
        except ValueError:
            logger.warning("given state is not integer!!")
            return -1

        #TODO: extract surface height to another function
        try:
            atoms = self.cmd.get_model("surface")
        except:
            logger.warning("Can not find model: surface")
            return 0

        for atom in atoms.atom:
            if atom.coord[2] > stop:
                stop = atom.coord[2]
                atop = atom.coord[:3]
                rtop = str(atom.resi)
                ntop = str(atom.name)
        logger.debug("(state: {})sTOP is : {}.".format(self.cmd.get_state(),stop))

        model = ""
        get_state = self.cmd.get_state()
        if state == 0:
            model = "protein"
            get_state = 0
        elif state != -1:
            model = "mov"
            get_state = state
        else:
            model = "mov"
        try:
            logger.debug("Trying to get model: {} at state({}).".format(model, get_state))
            atoms = self.cmd.get_model(model, get_state)
        except:
            logger.warning("Can not find model({}): in state({})".format(model, get_state))
            return 0

        for atom in atoms.atom:
            if atom.coord[2] < pdown:
                pdown = atom.coord[2]
                adown = atom.coord[:3]
                rdown = str(atom.resi)
                ndown = str(atom.name)
        logger.debug("(state: {})pdown is : {}.".format(self.cmd.get_state(), pdown))

        logger.debug("Calculate the destince bewteen two atoms.")
        deltasq = [(i-j)**2 for i,j in zip(atop, adown)]
        distanc = math.sqrt(sum(deltasq))
        distanz = pdown - stop
        logger.debug("The distance is : {}".format(distanc))
        logger.debug("Z distance is : {}".format(distanz))

        r = distanc
        if z == True:
            r = distanz
        return r

    def surfacez(self):
        '''Find the highest Z value of the atoms in the surface'''
        stop  = float("-inf")
        atop  = [0.0, 0.0, 0.0]
        rtop  = "" #resi of atom
        ntop  = "" #name of atom

        try:
            logger.debug("Try getting model: surface")
            atoms = self.cmd.get_model("surface")
        except:
            logger.warning("Can not find model: surface")
            return 0

        for atom in atoms.atom:
            if atom.coord[2] > stop:
                stop = atom.coord[2]
                atop = atom.coord[:3]
                rtop = str(atom.resi)
                ntop = str(atom.name)

        logger.debug("Found stop is: {}".format(stop))
        return stop

    def proteinz(self, state=0):
        '''Find the lowest Z value of the atoms in the protein'''
        pdown = float("+inf")
        adown = [0.0, 0.0, 0.0]
        rdown = "" #resi of atom
        ndown = "" #name of atom

        try:
            int(state)
        except ValueError:
            logger.warning("given state is not integer!!")
            return -1

        model = ""
        get_state = self.cmd.get_state()
        if state == 0:
            model = "protein"
            get_state = 0
        elif state != -1:
            model = "mov"
            get_state = state
        else:
            model = "mov"
        try:
            logger.debug("Try getting model: {}".format(model))
            atoms = self.cmd.get_model(model, get_state)
        except:
            logger.warning("Can not find model({}): in state({})".format(model, get_state))
            return 0

        for atom in atoms.atom:
            if atom.coord[2] < pdown:
                pdown = atom.coord[2]
                adown = atom.coord[:3]
                rdown = str(atom.resi)
                ndown = str(atom.name)

        logger.debug("Found pdown is: {}".format(pdown))
        return pdown

    def pdown(self):
        s = self.args.stepd
        s = 10      #roughly defined, to be changed.
        distanc = self.mindst(z=True, state=0)
        if s > distanc / 2:
            s = distanc / 2     #avoid move too far
        if s >= 1:
            #avoid move too short
            self.cmd.translate([0,0,s], "surface", state=0, camera=0)
            logger.debug("[pdown] Move the surface {} units.".format(s))

        #check crashing
        distanc = self.mindst(z=True, state=0)

        if distanc >= self.atom_cutoff:
            return distanc
        else:
            return 0    #Before Crashing

    def covaxis(self, axis):
        """convert the axis from number 1 2 3 to x y z"""
        if axis == 1 or axis == 2 or axis == 3:
            return self.xyzcv[axis]
        else:
            return axis

    def pcheck(self, axis, step, action):
        logger.debug("checking protein move action {}".format(action))
        if action != "ptransl" and action != "protate":
            logger.warning("move protein action should only be ptransl or protate!")
            return -1

        if axis != "x" and axis != "y" and axis != "z":
            logger.warning("[{}] The axis should be either ''x'' or ''y'' or ''z''.".format(action))
            return -1

        try:
            float(step)
        except ValueError:
            logger.warning("[{}] The step should only be float number.".format(action))
            return -1

        logger.debug("The axis and step are all right, continue...")
        return 0

    def ptransl(self, axis, step):
        if self.pcheck(axis, step, "ptransl"): return -1

        self.cmd.center("protein")
        center = self.cmd.get_position()
        logger.debug("[ptransl] protein position BEFORE translation: {}".format(center))
        delta = [0.0,0.0,0.0]
        delta[self.xyzid[axis]] = step
        logger.debug("[ptransl] translate protein delta: {}".format(delta))
        self.cmd.translate(delta, "protein", state=0, camera=0)
        self.cmd.center("protein")
        center = self.cmd.get_position()
        logger.debug("[ptransl] protein position AFTER translation: {}".format(center))

    def protate(self, axis, step):
        if self.pcheck(axis, step, "protate"): return -1

        za = self.proteinz()
        logger.debug("[protate] The lowest Z B4 rotate is : {}".format(za))

        self.cmd.center("protein")
        origin = self.cmd.get_position()
        logger.debug("[protate] define the origin: {}".format(origin))
        self.cmd.rotate(axis, step, "protein", camera=0, origin=origin)
        logger.debug("[protate] Rotate {} axis {} units ".format(axis, step))
        logger.debug("rotate {}, {}, protein, camera=0, origin={}".format(axis, step, origin))

        zb = self.proteinz()
        logger.debug("[protate] The lowest Z AF rotate is : {}".format(zb))

        delta = za - zb
        logger.debug("[protate] distA: {} - distB: {} = delta: {}".format(za, zb, delta))

        self.cmd.translate([0,0,delta], "protein", state=0, camera=0)
        logger.debug("[protate] translated the protein {}(delta) units.".format(delta))
        logger.debug("[protate] now the Z value is : {}".format(self.proteinz()))

    def protation(self):
        s = self.args.steps
        step = 360.0 / s
        sum = 1

        for i in range(s):

            for j in range(s):

                #save state
                self.cmd.create("mov", "protein", self.cmd.get_state(), sum)
                logger.debug("[Current state ({})]Saved new state (count: {}), after rotate x".format(self.cmd.get_state(), i*s+j+1))
                logger.debug("[protation] for loop i: {}, j: {}, ia: {}, ja: {}".format(i,j, i*step, j*step))

                logger.debug("{:=^70}".format("next round"))
                self.protate("x", step)
                sum+=1

            logger.debug("{:=^70}".format("Y AXIS"))
            self.protate("y", step)


        logger.debug("{:=^70}".format("N TURN"))
        self.cmd.ending()
        logger.debug("[protation] Goes to the last state ({}).".format(self.cmd.get_state()))
        self.cmd.create("final", "protein", self.cmd.get_state())
        logger.debug("Saved the final structure.")
        self.savefile(pdb=True)
        logger.info("[protation] Saved file (count: {}), rotation sampling is done".format(self.svcount))



    def simulate(self):
        '''run sampling and then scoring'''
        self.sampling()
        self.scoring()

    def scoring(self):
        logger.info("Now we call gromacs to calculate the score.")

        for con in self.jsdb["confs"]:
            f = con["id"] + ".pdb"
            logger.debug("processing file: {}".format(f))
            fn = os.path.join(self.jsdb["confd"], f)
            logger.debug("scoringone the file: {}".format(fn))
            r = self.scoringone(fn)
            try:
                con["stepN"] = r["stepN"]
                con["coul"] = r["coul"]
                con["ljsr"] = r["ljsr"]
                con["file"] = r["file"]
                logger.debug("updated the configuration db, {}:{}".format(con["coul"], con["ljsr"]))
            except:
                logger.warning("scoringone return incorrect data!!!")


        self.savejdb()


    def scoringone(self, filename, mdir=None):
        if mdir == None:
            mdir = self.mdir

        logger.debug("This is scoring ONE, the filename: {}".format(filename))
        for filei in [filename, self.runshp]:
            if not os.path.isfile(filei):
                logger.warning("The file({}) does NOT exist! Please Check!".format(filei))
                return

        last = ""
        line = ""   #in case the energy file is empty.
        rpath = os.path.realpath(filename)
        slogf = open(self.slogn, "a")
        call([self.bash, self.runshp, rpath], stdout=slogf, stderr=slogf)
        #ASSUMED the script will run and exit correctly.
        #AND output absolutely right answer.
        try:
            with open(os.path.join(mdir, self.energyx), "r") as openfile:
                # just read the last line
                for line in openfile:   pass
                last = line
        except IOError:
            logger.critical("SCORINGONE ERROR !!!")
            logger.critical("Unable to read the file: {}!".format(\
                    os.path.join(mdir, self.energyx)))
            last = ""
            sys.exit(1)
        logger.debug("The last line: {}".format(last))
        jsona = {}
        lasts = last.split()
        if len(lasts) >= 3:
            jsona["file"] = filename
            jsona["lasta"] = lasts[0:3]
            try:
                jsona["stepN"] = int(float(lasts[0]))
                jsona["coul"] = float(lasts[1])
                jsona["ljsr"] = float(lasts[2])
            except ValueError:
                jsona["stepN"] = 0
                jsona["coul"] = 0.0
                jsona["ljsr"] = 0.0
                logger.critical("SCORINGONE ERROR !!!")
                logger.critical("Unable to read the values from energy file!")
                sys.exit(1)
            jsona["energy"] = jsona["coul"] + jsona["ljsr"]
            logger.debug("scoringone energy: {}".format(jsona["energy"]))
        else:
            jsona["stepN"] = 0
            jsona["lasta"] = ["0","1","2"]
            jsona["energy"] = 0.0
            logger.warning("scoringone, spliting last line of energy file failed!")
        logger.debug("splitted last line, should be the energy values: [{}]".format(",".join(jsona["lasta"])))

        return jsona

    def sampling(self):
        r = random.random()
        logger.info("the random number is :"+ str(r))

        self.cmd.do('print "Hajimaru yo ~♪"')
        self.load_file()    # Load surface object and protein object

        logger.debug("SECRET thing, rotate X/Y slightly... {}/{} degrees...".format(self.args.init_xr, self.args.init_yr))
        self.protate("x", self.args.init_xr)
        self.protate("y", self.args.init_yr)
        logger.info("ONE more thing, translate X/Y {}/{} and rotate Z {}".format(self.args.init_xt, self.args.init_yt, self.args.init_zr))
        self.protate("z", self.args.init_zr)
        self.ptransl("x", self.args.init_xt)
        self.ptransl("y", self.args.init_yt)
        self.ptransl("z", self.args.init_height)
        logger.info("FINALLY, The minimal Z distance after adjustment is : {}".format(self.jsdb["mindst"]))
        self.savefile()
        logger.info("Saved file (count: {}), after rotation and translation".format(self.svcount))


        # Ready to start sampling
        r = 0   # count the number of runs
        c = 1   # when c is true, continue.

        while r < self.args.stepd and c != 0:
            self.protation()    #rotate the protein
            c = 0               #now we only run once protation
            r+=1
        if c != 0:
            self.savefile(pdb=True)
            logger.info("Saved file (count: {}), final configuration.".format(self.svcount))
        else:
            logger.info("The final configuration should be the last count ({}).".format(self.svcount))
        self.jsdb["mindst"] = self.mindst(state=0, z=True)
        logger.debug("updated the mini distance Z: {}".format(self.jsdb["mindst"]))
        self.jsdb["svcount"] = self.svcount

        with open(self.args.jsdbf, "w") as openfile:
            json.dump(self.jsdb, openfile)
        logger.info("Saved the jsdb file, after sampling all confs!")

    def plot(self):
        '''output the db.json to csv'''
        if self.args.no_plot:
            logger.info("skipped the plot progress...")
            return 0

        try:
            a = self.jsdb["confs"]
            steps = self.jsdb["steps"]
            mindst = self.jsdb["mindst"]
        except KeyError:
            logger.error("Can not find the confs/steps from the json db file. It might be corrupted.")
            return 1


        count = 0

        x = []
        y = []
        z = []
        xt = []
        yt = []
        zt = []
        ang = 360.0 / steps
        lx = []
        ly = []
        ld = []
        minres = "UNKNOW"
        mineng = 7777777
        mincnt = 0
        maxeng = -7777777
        maxres = "UNKNOW"
        minxra = 0
        minyra = 0

        logger.debug("preparing X Y Z")
        for i in a:
            try:
                enegy = float(i["coul"]) + float(i["ljsr"])
                resnm = i["id"]
                xi = int(i["xi"])
                yi = int(i["yi"])
            except (ValueError, KeyError):
                enegy = 7777777
                resnm = "UNKNOW"
                xi = 0
                yi = 0
            if enegy < mineng:
                mineng = enegy
                minres = resnm
                mincnt = count
                minxra = xi
                minyra = yi
                logger.debug("Found a lower energy residue: {}@{}".format(minres, mineng))
            if enegy > maxeng:
                maxeng = enegy
                logger.debug("Found a higher energy residue: {}@{}".format(maxres, maxeng))

            xti = float(xi * ang - self.args.init_xr)
            yti = float(yi * ang - self.args.init_yr)
            xt.append( xti )
            yt.append( yti )
            zt.append(enegy)
            lx.append(count)
            ly.append(enegy)
            ld.append(resnm)
            count += 1
            if (count % steps) == 0:
                x.append(xt)
                xt = []
                y.append(yt)
                yt = []
                z.append(zt)
                zt = []
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        logger.debug("get X, Y and Z, then plot the graph.")

        figa = plt.figure()
        cmap = plt.cm.get_cmap("RdBu_r")
        if self.args.minlv != None and self.args.maxlv != None:
            a = self.args.minlv
            b = self.args.maxlv
            i = (b - a) / self.ctlvs    #ASSUMED: the maxlv is always bigger than minlv
            levels = range(a, b, i)
            logger.debug("The contour map levels: {}".format(levels))
        else:
            levels=None
            logger.debug("Use the default contour map levels.")
        plt.plot(112, -65, 'ko')   #the PSO ans
        #plot the minimum location
        plt.plot( float(minxra * ang - self.args.init_xr), float(minyra * ang - self.args.init_yr), 'kx')
        Da = plt.contourf(x, y, z, cmap=cmap, levels=levels, norm=mpl.colors.SymLogNorm(011))
        plt.title(u"contour diagram\ndistance={}Å".format(mindst))
        plt.xlabel("X rotation angle")
        plt.ylabel("Y rotation angle")
        cbar = plt.colorbar(Da)
        cbar.ax.set_ylabel("energy level")
        plt.savefig("diagram_0a.pdf")
        logger.debug("plot contour diagram and save as pdf file.")

        figb = plt.figure()
        plt.title(u"energy line\ndistance={}Å".format(mindst))
        plt.xlabel("Iteration Number")
        plt.ylabel("Energy Value")
        Db = plt.plot(lx, ly, 'k')
        plt.plot(mincnt, mineng, 'bo')
        plt.plot([mincnt, mincnt*1.1], [mineng, mineng], 'k')
        plt.text(mincnt*1.13, mineng, "id: {}\nen: {}".format(minres, mineng), verticalalignment="center", horizontalalignment="left")
        logger.debug("plot energy line only.")
        if self.args.minlv != None:
            plt.ylim(ymin=self.args.minlv)
            logger.debug("set the y-axis minimum range.")
        if self.args.maxlv != None:
            plt.ylim(ymax=self.args.maxlv)
            logger.debug("set the y-axis maximum range.")
        plt.savefig("diagram_0k.pdf")
        ly = np.array(ly)
        ld = np.array(ld)
        L = np.column_stack((ly, ld))
        np.savetxt("energy.txt.gz", L, delimiter=" ", fmt="%11s %11s")
        logger.debug("plot energy line diagram and save as pdf file.")

        fige = plt.figure()
        plt.title(u"normalized energy line\ndistance={}Å".format(mindst))
        plt.ylabel("Energy Value")
        plt.xlabel("residues")
        plt.axis("off")
        plt.grid("on")
        plt.xticks([])
        plt.yticks([])
        sly = np.sort(ly)
        sry = sly[::-1]
        nly = (sry - mineng) / (maxeng - mineng)
        De = plt.plot(range(len(nly)), nly)
        plt.text(0,0, "id: {}@{}".format(minres, mineng))
        logger.debug("plot energy histogram.")
        plt.savefig("diagram_0e.pdf")
        logger.debug("plot energy histogram diagram and saved as pdf file.")

        figh = plt.figure()
        plt.title(u"normalized energy line\ndistance={}Å".format(mindst))
        Dh = plt.hist(ly, 100)
        plt.xlabel("The lowest configuration is id: {}@eng: {}".format(minres, mineng))
        plt.savefig("diagram_0h.pdf")
        logger.debug("plot another histogram diagram and saved as pdf file.")


        if self.args.no_plot_re:
            logger.info("skipped the residues plot progress...")
            return 0
        logger.debug("plot residues configuration graph, go through all configurations...")
        self.cmd.load("system00003.pse")
        for con in a:
            cid = con["id"]
            try:
                nid =  re.findall('\d+', cid)[0]
                nid = int(nid)
            except:
                logger.warning("function plot - can not find digit from cid.")
                nid = 0
            logger.debug("find the digit: {} from cid.".format(nid))
            cen = float(con["coul"]) + float(con["ljsr"])
            logger.debug("processing conf: {}, the energy value is: {}".format(cid, cen))
            logger.debug("Now create the system state for processing...")
            self.cmd.create("system", "mov", nid+1, 1)
            atoms = self.cmd.get_model("system")
            resds = atoms.get_residues()
            xr = []         # x axis: residues id
            yr = []         # y axis: mini distance
            mindt = [999999.99] * len(resds)
            surfc = self.surfaceh()
            szv = surfc["stop"]
            for atom in atoms.atom:
                rid = int(atom.resi) - 1
                azv = atom.coord[2]
                dst = azv - szv
                if dst < mindt[rid]:
                    mindt[rid] = dst

            #save the graph under jsdb["confd"] directory
            figs = plt.figure()
            Ds = plt.plot(range(1, len(mindt)+1), mindt, 'r+', \
                    range(1, len(mindt)+1), mindt, 'k')
            plt.title(u"residues configuraions diagram\ndistance={}Å; energy={}kj".format(self.jsdb["mindst"], cen))
            plt.savefig(os.path.join(self.jsdb["confd"], cid+".pdf"))
            logger.debug("plot a residues diagram and saved as pdf.")


    def surfaceh(self):
        logger.debug("calculating the sruface top position.")
        surf = {}
        stop  = float("-inf")
        atop  = [0.0, 0.0, 0.0]
        rtop  = "" #resi of atom
        ntop  = "" #name of atom

        try:
            atoms = self.cmd.get_model("surface")
        except:
            logger.warning("Can not find model: surface")
            return 0

        for atom in atoms.atom:
            if atom.coord[2] > stop:
                stop = atom.coord[2]
                atop = atom.coord[:3]
                rtop = str(atom.resi)
                ntop = str(atom.name)
        logger.debug("(state: {})sTOP is : {}.".format(self.cmd.get_state(),stop))
        surf["stop"] = stop
        surf["atop"] = atop
        surf["rtop"] = rtop
        surf["ntop"] = ntop
        logger.debug("return surface dict.")
        return surf


    def cleaning(self):
        flist = ["db.json", "debug.log"]
        dlist = ["conf"]

    def helping(self):
        return


    def main(self):
        # TODO: assume all files exist.
        # checking...
        # self.args.proteinf
        # self.args.surfacef
        # self.args.topolf

        logger.info("Starting the main function...")
        #act = getattr(self, self.args.act, self.helping)
        #print("printing self.args:")
        #print(self.args)
        #if callable(act):
            #logger.info("Run the action ({})".format(str(act)))
            #act()


if __name__ == "__main__":
    simpleMOVE().main()
