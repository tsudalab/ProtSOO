#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Written by: Jimmy NGAI jimmycfngai@gmail.com

from __future__ import print_function       # for python 2.7
from subprocess import call
import sys, os, argparse, math, re
from collections import Counter
import pprint
import pdb

# These are the "Color Blind 10" colors as RGB.  
#blind10=[(255,128,14),(171,171,171),(95,158,209),(89,89,89),(0,107,164),(255,188,121),(207,207,207),(200,82,0),(162,200,236),(137,137,137)]
blind10=[(0,107,164),(200,82,0),(255,128,14),(137,137,137),(171,171,171),(162,200,236),(89,89,89),(255,188,121),(95,158,209),(207,207,207)]
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

from measure import measure

from simpleLOG import logger



class dbcluster:
    def __init__(self):
        self.loacts = ['helping', 'cluster']
        self.parser = argparse.ArgumentParser(description="This is measure.py")
        self.parser.add_argument("act", help="run something main/etc...", choices=self.loacts)
        self.parser.add_argument("--gbs", help="load the gbest txt files.", nargs="*", default=[], type=str)
        self.parser.add_argument("--eps", help="The maximum distance between two samples for them to be considered as in the same neighborhood.", default=0.3, type=float)
        self.parser.add_argument("--minn", help="MIN", default=1, type=int)
        self.parser.add_argument("--out-dir", help="output directory", type=str, default="cluster")
        self.parser.add_argument("--no-plot", help="skip the plotting", action="store_true", default=False)
        self.args = self.parser.parse_args()

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

        self.crank = -1
        self.cthem = blind10
        self.clenf = len(self.cthem)
        self.curc = ""  #current color code 

        self.ylimp = -7777777
        self.ylimn = +7777777


        """
        initialized variables, data structure:

        index: *gbest_id*           #for each gbest result
        self.labels = [ "gbest_filename", ... ]   #gbest filename list
        self.Hnums = [ 
                        [ 1.0, 2.0, 3.0, ... ],   #the residue height
                        ...
                    ]
        self.dbresult = [ *cluster_id*, ... ]

        self.clusters = { 
                        *cluster_id* : {
                            coden : "01",                   #codename, ordered by lowest energy
                            gbest : [ *gbest_id*, ... ],    #list of gbest id
                            englt : [ -100.0, ... ],        #list of gbest energy
                            adsit : [ 1, 2, 3, ...],        #common adsorption site
                            adslt : [ 1, 2, 3, ...],        #all adsorption site list
                            adsle : [ -100.00, ...],        #the energy values of each in adslt
                            fpath : "path/to/cluster.pdf"   #cluster plotting file
                            elwid : 1,                      #lowest energy gbest id
                            englw : -900,                   #loewst energy value
                            engmd : -800                    #cluster median energy value
                            },
                        ...
                        }

        #Clusters id ordered by lowest energy
        index: *byenergy_id*
        self.cs_eid = [ *cluster_id*, ... ]
        """

        self.Hnums = []
        self.labels = []
        self.clusters = {}
        self.clusters_byenergy = []
        self.parametbar = ["pso50", "pso50-200", "pso100", "pso150", "pso200", "pso250"]


        self.pp = pprint.PrettyPrinter(indent=2)
        


    def nxc(self):
        """get the index of next color"""
        self.crank += 1
        self.crank %= self.clenf
        #convert into matplotlib accepted format
        r = [ x/255.0 for x in self.cthem[self.crank] ]
        return r

    def hello(self, msg="はじまるよ～♪ "):
        logger.debug("[hello] {}".format(msg))
        logger.debug("hello function is deprecated.")






    def readgbest(self, f):
        try:
            L = np.loadtxt(f, delimiter=',')
            H = L[:,1]
        except:
            self.hello("Can not load the file: {}".format(f))
            return

        if len(self.Hnums) > 0:
            self.Hnums = np.vstack((self.Hnums,H))
        else:
            self.Hnums = H
        self.labels.append(f)

    def get_gbest(self):
        """given a list of gbest data filename
        In each files, there are two numbers per lines.
        the first number is the residue number, and the second one is the distance to surface.
        Load the number and create a matrix of all the vectors, with a corresponding filename array.
        """
        for f in self.args.gbs:
            self.readgbest(f)


    def avgenglst(self):
        energy = []
        stdenergy = []
        lowenergyv = []
        lowenergyp = []
        medianengy = []
        for idx, cst in enumerate(self.clusters):
            if len(cst) == 0:
                #pass if no noise group
                energy.append(0)
                lowenergyv.append(0)
                lowenergyp.append(0)
                medianengy.append(0)
                continue
            englst = self.avgenergy(cst)
            energy.append(np.average(englst))
            stdenergy.append(np.std(englst))
            #energy.append(self.avgenergy(cst))
            #self.hello("DEBUG {:02d}th cluster lowest energy conformation is: {}".format(idx, englst.min()))
            lowenergyv.append(englst.min())
            lowenergyp.append(englst.argmin())
            medianengy.append(np.median(englst))

        #TODO: ugly! should not be created in this way!
        self.lowest_eng = lowenergyv[:]
        self.low_energy = lowenergyp[:]
        self.avg_energy = energy[:]
        self.median_eng = medianengy[:]
        self.std_energy = stdenergy[:]
        return energy
        
    def avgenergy(self, paras):
        """
        paras: list of gbest.txt file path

        load the list of gbest final energy values (in np array)
        the average value can be computed through numpy methods
        """
        psoeng = [ re.sub("gbest.txt", "gbest_energy.txt", x) for x in paras ]
        engval = []
        for f in psoeng:
            try:
                v = np.loadtxt(f, delimiter=',')
            except:
                logger.error("Can NOT open the numbers in textfile: {}".format(f))
                v = []

            if len(v) > 0:
                engval.append(v[-1])
            else:
                logger.error("Can NOT open the numbers in textfile: {}".format(f))

        #engval = np.array(engval)
        return np.array(engval)
        #return np.average(engval)

    def loadadsite(self, gbestf, cutoff=5.0):
        """
        gbestf is gbest.txt file path
        read it and filter out all residues within cut-off
        """
        gbests = re.sub("gbest.txt", "gbest_sorted.txt", gbestf)
        try:
            L = np.loadtxt(gbests, delimiter=",")
            logger.debug("Loaded gbest: {}".format(gbests))
        except:
            logger.error("Can NOT read the numbers in {}".format(gbestf))
            L = []

        R = []
        if len(L) > 0 and len(L[0]) > 1:
            #The residues height was stored in second column
            min_h = L[:,1].min() * 1.0
            n_cutoff = min_h + cutoff
            #filter out with-in cut-off
            R = L[L[:,1] < n_cutoff][:,0]
            logger.debug("filter cutoff={}: {}".format(n_cutoff, str(R)))

        return R

    def loadenergy(self, gbestf):
        """
        gbestf if the gbest.txt file path.
        need to convert to gbest_energy.txt
        return the gbest final energy.
        """
        gbeste = re.sub("gbest.txt", "gbest_energy.txt", gbestf)
        try:
            v = np.loadtxt(gbeste, delimiter=",")
            laste = v[-1][-1]
        except:
            logger.error("Can NOT read the energy values of {}".format(gbeste))
            laste = 0

        return int(laste)

    def loadreseng(self, gbestf, rid):
        """
        call score_res.sh to compute the energy between the residue rid and surface
        """
        gbestp = re.sub("gbest.txt", "gbest.pdb", gbestf)
        reshp = ""
        slogn = "score_res.log"
        scoringf="score_res.sh"
        bash = "/bin/bash"
        mdir = "EM"
        energyx = "energy.xvg"
        logger.debug("searching {}...".format(scoringf))
        for sysp in sys.path:
            for root, dirs, files in os.walk(sysp):
                for f in files:
                    if f == scoringf:
                        runshp = os.path.abspath(os.path.join(root, f))
                        logger.debug("{} was found!".format(runshp))
                        # it will return the last one found
        logger.debug("Going to compute the energy of residue: {}".format(rid))
        for filei in [gbestp, runshp]:
            if not os.path.isfile(filei):
                logger.warning("The file({}) does NOT exist! Please Check!".format(filei))
                return

        last = ""
        line = ""   #in case the energy file is empty.
        rpath = os.path.realpath(gbestp)
        slogf = open(slogn, "w")    #Ugly: simply overwritten existing log
        #ASSUMED
        call([bash, runshp, rpath, str(int(rid))], stdout=slogf, stderr=slogf)
        #ASSUMED the script will run and exit correctly.
        #AND output absolutely right answer.
        try:
            with open(os.path.join(mdir, energyx), "r") as openfile:
                # just read the last line 
                for line in openfile:   pass
                last = line
        except IOError:
            logger.critical("SCORINGONE ERROR !!!")
            logger.critical("Unable to read the file: {}!".format(\
                    os.path.join(mdir, energyx)))
            last = ""
            sys.exit(1)
        logger.debug("The last line: {}".format(last))
        jsona = {}
        lasts = last.split()
        if len(lasts) >= 3:
            try:
                stepN = int(float(lasts[0]))
                coul  = float(lasts[1])
                ljsr  = float(lasts[2])
            except:
                logger.critical("Can't read the energy!!")
                stepN = coul = ljsr = 0.0

        logger.debug("Read the energy of residue ({}), Coul: {}, LJ: {}".format(rid, coul, ljsr))
        with open("residues_energy.txt", "a") as fres:
            fres.write("{}, {}, {}, {}\n".format(gbestp, rid, coul, ljsr))



    def clustersub(self, clusterlist, substring):
        #clusterbar = {}
        count = 0.0
        for item in clusterlist:
            if substring in item:
                count += 1
            #clusterbar[substring] = count
        return count

    def clusterbar(self):
        self.hello("cluster id is according to the size.")
        clusterbar = {}
        for n in self.parametbar:
            clusterbar[n] = []
        for nth, idx in enumerate(self.rankings):
            if idx < len(self.clusters):
                cst = self.clusters[idx]
            else:
                cst = []
            for pso in self.parametbar:
                count = self.clustersub(cst, pso+"/")
                clusterbar[pso].append(count)
            #self.hello("cluster[{:02d}] = {}".format(nth+1, clusterbar))
        self.hello(str(clusterbar))

    def ranking(self):
        self.rankings = []
        self.rankinge = []
        average_energy = self.avgenglst()
        for cst in range(len(self.clusters)-1):
            csize = len(self.clusters[cst])
            energy = average_energy[cst]

            #ranking by size first then by average energy
            tgidx = 0
            while tgidx < len(self.rankings) and self.rankings[tgidx] < len(self.clusters) \
                and len(self.clusters[self.rankings[tgidx]]) > csize:
                tgidx += 1
            while tgidx < len(self.rankings) and self.rankings[tgidx] < len(self.clusters) \
                and len(self.clusters[self.rankings[tgidx]]) == csize    \
                and tgidx < len(self.rankings) and self.rankings[tgidx] < len(average_energy) \
                and average_energy[self.rankings[tgidx]] < energy:
                #remember here the energy is negative
                tgidx += 1
            self.rankings.insert(tgidx, cst)

            #ranking by average energy first then by size
            tgidx = 0
            while tgidx < len(self.rankinge) and self.rankinge[tgidx] < len(average_energy) \
                and average_energy[self.rankinge[tgidx]] < energy:
                #remember here the energy is negative
                tgidx += 1
            while tgidx < len(self.rankinge) and self.rankinge[tgidx] < len(average_energy) \
                and average_energy[self.rankinge[tgidx]] == energy   \
                and tgidx < len(self.rankinge) and self.rankinge[tgidx] < len(self.clusters) \
                and len(self.clusters[self.rankinge[tgidx]]) > csize:
                tgidx += 1
            self.rankinge.insert(tgidx, cst)

            #self.hello("DEBUG rankings: {}".format(self.rankings))
            #self.hello("DEBUG rankinge: {}".format(self.rankinge))


    def printoutcluster(self):
        """
        example output:

Assume homogeneous surface, performing clustering analysis based on the residue min-distance profiles:

=== Sorted by cluster size ===
Cluster      Size      Average E     Median E      Lowest E          Highest E           Min-distance profile
01               35       -842.84335    -855.28589   -1002.338455   highest energy   give the path to the profile of the cluster plot
02               10       -788.97523    -806.48693   -1015.307678   highest energy   give the path to the profile of the cluster plot
...

==== Cluster 01 ====
#  ProtPOS Score        Predicted PDB
1  -1002.338455            ./200/trun-PSO-072112-test03/gbest.pdb
2  - 800.3949830           ./250/trun-PSO-072112-test03/gbest.pdb
3  - 700.3949830           ./210/trun-PSO-072112-test03/gbest.pdb
4  - 600.3949830           ./220/trun-PSO-072112-test03/gbest.pdb
5  - 500.3949830           ./240/trun-PSO-072112-test03/gbest.pdb
...

Contacting Residues Within 5 Angstrom from Surface:
ResID  Freq   %          Avg Coulomb       Avg LJ
128      35     100.00    -123.232312        -232.929392   <-- i.e. occurrence frequency in all cluster members, the percentage, average coulomb (of all members) and LJ energy 
20        25       71.43    -123.232312        -232.929392  <-- =25/35*100 = 71.42857... round to 2 decimal places 
18        10       28.57    -123.232312        -232.929392  <-- =10/35*100 = 28.571428... round to 2 decimal places


==== Cluster 02 ====
...

        """
        logger.info("Assume homogeneous surface, performing clustering analysis based on the residue min-distance profiles:")
        #logger.info("=== Sorted by cluster size ===")
        logger.info("=== Clusters summary ===\n")
        labels = ["Cluster id", "Size", "Average E", "Median E", "Lowest E", "Highest E", "Min-distance profile"]
        labelp = ["{:<12s}"] + ["{:^7s}"] + ["{:^12s}"] * 4 + ["{}"]
        labelw = ["{:<12s}"] + ["{:^7d}"] + ["{:^12.3f}"] * 4 + ["{}"]
        outstr=""
        for i,j in zip(labelp, labels):
            outstr+=(i.format(j))
        logger.info(outstr)
        #for c in self.clusters:
        for c in self.cs_eid:
            cs = self.clusters[c]
            enl = np.array(cs["englt"])
            csl = [cs["coden"], enl.size, np.average(enl), np.median(enl), enl.min(), enl.max(), cs["fpath"]]
            outstr = ""
            for i,j in zip(labelw, csl):
                outstr += i.format(j)
            logger.info(outstr)

        logger.info("\nEnergy (E): in kJ/mol")
        logger.info("noise: not grouped into any cluster\n")


        for c in self.cs_eid:
            cs = self.clusters[c]
            logger.info("\n=== Cluster {} ===\n".format(cs["coden"]))
            englt = cs["englt"]
            csize = len(englt) * 1.0
            gbest = cs["gbest"]     #this is gbest id
            gbest_pdb = [re.sub("\.txt", ".pdb", self.labels[x]) for x in gbest]
            cs_list = zip(englt, gbest_pdb)
            cs_sorted = sorted(cs_list, key=lambda x: x[0], reverse=False)
            idslt = range(1, len(englt)+1)
            logger.info("{:<2s} {:^15s} {}".format("#", "ProtPOS Score", "Predicted PDB"))
            #for i,j,k in zip(idslt, englt, gbest_pdb):
            for i,j in zip(idslt, cs_sorted):
                #ASSUMED j has two elements
                logger.info("{:<2d} {:^15.5f} {}".format(i, j[0], j[1]))
            logger.info("\nContacting Residues Within 5 Angstrom from Surface:")
            logger.info(",".join(["{:.0f}".format(x) for x in cs["adsit"]]))
            #logger.info("Common adsorption site: {}".format(str(cs["adsit"])))

            #logger.info("Contacting Residues Within 5 Angstrom from Surface:")
            logger.info("\nResID  Freq         %          ")
            cc = Counter(cs["adslt"])
            cc_sorted = sorted(cc, key=lambda x: cc[x], reverse=True)
            for i in cc_sorted:
                logger.info("  {:03.0f}   {:>3d}   {:06.2f}%".format(i, cc[i], cc[i]/csize*100))


    def plotcluster(self, cid, paras, bold="", justprint=False):
        """
        cid: cluster plotting filename
        paras: list of gbest.txt file path
        bold: blod line for the lowest energy profile
        """
        outdir = re.sub("/+", "", self.args.out_dir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        #filename=os.path.join(outdir, "{}-cluster{:02d}.pdf".format(byrank, cid))
        filename=os.path.join(outdir, "cluster-{}.pdf".format(str(cid)))
        bak_argv = sys.argv
        sys.argv = ['measure.py', 'energy', '--no-legend', '-b', '--savepdf', 
        '--txtsize', '24',
        '--xlabel', 'residue number', '--ylabel', 'minimum distance to surface',
        '--ylim', '0', '35', '-o', filename, '--xvg']
        #try:
            #paras = paras[:]    #remember to make a copy, otherwise it will modify the original clusters list
            #paras.remove(bold)
        #except ValueError:
            #pass
        paras = [ x for x in paras if x != bold ]
        sys.argv.extend(paras)
        if bold != "":
            sys.argv.extend(["--bxt-width", "02", "--bvg", bold])
        logger.debug("DEBUG calling: `{} \"{}\"`".format(sys.argv[0], '" "'.join(sys.argv[1:])))

        if justprint:
            #just print the debug message.
            return

        try:
            bak_out = sys.stdout
            bak_err = sys.stderr
            null_out = open(os.devnull, 'w')
            sys.stdout = null_out
            sys.stderr = null_out
            measure().main()
        except:
            self.hello("({})something wrong!".format(paras))
            raise
        finally:
            sys.argv = bak_argv
            sys.stdout = bak_out
            sys.stderr = bak_err
            null_out.close()

        return filename

    def cluster(self):
        """
        cluster by dbscae
        it will return a list with the corresponding cluster id
        we convert it into a dict index by cluster id

        """
        self.get_gbest()
        X = StandardScaler().fit_transform(self.Hnums)
        db = DBSCAN(eps=self.args.eps, min_samples=self.args.minn).fit(X)
        self.dbresult = np.array(db.labels_)
        logger.debug("DBSCAN labels: {}".format(",".join(map(str,self.dbresult))))
        """
        self.clusters = [ [] for x in range(labels.max()+2) ]
        for idx, val in enumerate(labels):
            if idx < len(self.labels) and val < len(self.clusters):
                self.clusters[val].append(self.labels[idx])
        self.ranking()
        if not self.args.no_plot and len(self.clusters) > 0:
            idx = 0
            while idx < len(self.rankings) and idx < len(self.rankinge):
                rankids = self.rankings[idx]
                rankide = self.rankinge[idx]
                loweids = self.low_energy[rankids]
                loweide = self.low_energy[rankide]
                s_clusters = self.clusters[rankids]
                e_clusters = self.clusters[rankide]
                #self.hello("DEBUG the lowest energy by size is {}".format(s_clusters[loweids]))
                #self.hello("DEBUG the lowest energy by energy is {}".format(e_clusters[loweide]))
                #self.plotcluster(idx+1, s_clusters, "bysize", s_clusters[loweids])
                self.plotcluster("{:02d}".format(idx+1), e_clusters, "byenergy", e_clusters[loweide])
                idx+=1
            self.hello("Plot the noise group")
            self.plotcluster("noise", self.clusters[-1], "noise")
        """
        for idx, cid in enumerate(self.dbresult):
            #Here, idx is *gbest_id*
            if cid not in self.clusters:
                self.clusters[cid] = {}
                self.clusters[cid]["coden"] = "{:02d}".format(cid+1)
                self.clusters[cid]["gbest"] = []
                self.clusters[cid]["englt"] = []
                self.clusters[cid]["adsit"] = []
                self.clusters[cid]["adslt"] = []

                if cid == -1:
                    self.clusters[cid]["coden"] = "noise"

            gbestf = self.labels[idx]
            energy = self.loadenergy(gbestf)

            adsita = self.loadadsite(gbestf)
            adsitb = self.clusters[cid]["adsit"]
            logger.debug("site a:{}".format(adsita))
            logger.debug("site b:{}".format(adsitb))
            if len(adsitb) > 0:
                adsitc = np.intersect1d(adsita, adsitb)
                logger.debug("site A and B intersection.")
            else:
                adsitc = adsita
            logger.debug("site c:{}".format(adsitc))
            self.clusters[cid]["adslt"] += list(adsita)


            self.clusters[cid]["gbest"].append(idx)
            self.clusters[cid]["englt"].append(energy)
            self.clusters[cid]["adsit"] = adsitc

        #ASSUMED: lambda should correct
        self.cs_eid = sorted(self.clusters, key=lambda x: np.mean(self.clusters[x]["englt"]), reverse=False)
        if -1 in self.cs_eid:
            self.cs_eid.remove(-1)
            self.cs_eid.append(-1)
        logger.debug("cluster sorted by energy: {}".format(str(self.cs_eid)))
        for i, c in enumerate(self.cs_eid):
            if c != -1:
                self.clusters[c]["coden"] = "{:02d}".format(i+1)

        if not self.args.no_plot and len(self.clusters) > 0:
            cidlt = self.clusters.keys()
            for c in cidlt:
                list_of_gbesti = self.clusters[c]["gbest"]
                list_of_gbestf = [ self.labels[x] for x in list_of_gbesti ]
                logger.debug("gbesti: {}".format(str(list_of_gbesti)))
                logger.debug("gbestf: {}".format(str(list_of_gbestf)))
                filename = "{}".format(self.clusters[c]["coden"])
                fpath = self.plotcluster(filename, list_of_gbestf, bold="")
                self.clusters[c]["fpath"] = fpath

        logger.debug("clusters: \n{}".format(self.pp.pformat(self.clusters)))
        """
        for idx, cst in enumerate(self.clusters):
            self.hello("DEBUG cluster id {:02d} average energy {} lowest energy {} size {}".format(idx, self.avg_energy[idx], self.lowest_eng[idx], len(cst)))
            self.hello("{}".format(cst))
        self.clusterbar()
        self.hello("{:=^70}".format("FINALLY"))
        self.hello("{:=^70}".format("by size"))
        for idx, cst in enumerate(self.rankings):
            cstlst = self.clusters[cst]
            loweng = self.low_energy[cst]
            #self.hello("The {:02d}th cluster (size={}): {}".format(idx+1, len(self.clusters[cst]), self.clusters[cst]))
            self.hello("The {:02d}th cluster (size={}; average energy={}; median energy={}; lowest energy={}; energy std={}): \n{}".format(\
            idx+1, len(cstlst), self.avg_energy[cst], self.median_eng[cst], self.lowest_eng[cst], self.std_energy[cst], cstlst[loweng]))
        self.hello("{:=^70}".format("by energy"))
        for idx, cst in enumerate(self.rankinge):
            cstlst = self.clusters[cst]
            loweng = self.low_energy[cst]
            #self.hello("The {:02d}th cluster (energy={}): {}".format(idx+1, self.avg_energy[cst], self.clusters[cst]))
            self.hello("The {:02d}th cluster (size={}; energy={}): {}".format(idx+1, len(cstlst), self.avg_energy[cst], cstlst[loweng]))
        """
        self.printoutcluster()

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
        logger.debug("command line arguments:")
        logger.debug(str(self.args))
        act = getattr(self, self.args.act, self.helping)
        if callable(act):
            logger.debug("Run the action ({})".format(str(act)))
            act()
            








if __name__ == "__main__":
    m=dbcluster()
    logger.debug("Hello! This is main function of dbcluster")
    m.hello()
    #m.fplot()
    #m.fneo()
    #m.fthree()
    #m.ffour()
    #m.ffive()
    m.main()

