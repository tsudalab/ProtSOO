
from __future__ import print_function       # for python 2.7
from subprocess import call
import sys, shutil, os, argparse, json, re
import random, math, errno

import simpleMOVE

import numpy as np

from simpleLOG import logger
from measure import measure
from math import floor
import time


class measurem(measure):
    """override self.cmd object, I can reuse the methods in measure."""
    def __init__(self, pymol_cmd):
        """override this"""
        self.cmd = pymol_cmd


class simplePSO(simpleMOVE.simpleMOVE):

    def __init__(self):
        super(simplePSO, self).__init__()
        self.measure = measurem(self.cmd)
        #logger.info("Initialized simplePSO object ")
    def initcmd(self):
        #logger.info("Initialized command line arguments")
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

    def initvar(self):
        super(simplePSO, self).initvar()

        logger.debug("initialized PSO variables: {}".format(self.__dict__))
        self.tid = 1

    def myPSO(self,x):

        #x = self.x[i]
        #for i in range(10):
            #x=[200,200,200,1,1,1]


            #print (x.shape())
        #print ("myPSO:",x)
        #v = self.v[i]
        #time9=time.time()
        #x=[219.50617283950612, 178.84773662551441, 64.93827160493828, -2.532235939643347, 0.7736625514403288, 1.1635802469135803]
        for idx, step in enumerate(x[:3]):
            self.protate(self.xyzcv[idx], step)
        for idx, step in enumerate(x[3:]):
            self.ptransl(self.xyzcv[idx], step)
        fpath = self.savers()
        jsons = self.scoringone(fpath, mdir=self.args.emdir)
            #print ("jsons:",jsons)
        s = jsons["energy"] * -1
        #timeo=time.time()-time9
        #print (timeo)
        #print ("energy:",s)
        jsona = {}
        #jsona["iteration"] = total
        #jsona["bird"] = i
        jsona["energy"] = s
        #jsona["position"] = x.tolist()
        #jsona["velocity"] = v.tolist()
        #jsona["gbest"] = False
        jsona["fpath"] = fpath
        jsona["gbest"] = True
        shutil.copyfile(fpath, "gbest.pdb")

        return s

    def main(self):

        self.load_file()
        self.load_post(self.args.init)
        #self.myPSO()
        #SOO()

def rescaling(lmin,lmax,real):
    x_sc=real*(lmax-lmin)+lmin

    return x_sc
def split8(x,y):
    element=(y-x)/3.0
    return element
def split(x1left,x1right,x2left,x2right):
    left_mid=(x1left+x1right)/2.0
    right_mid=(x2left+x2right)/2.0
    return left_mid,right_mid
def midpoint(p1, p2):
    mid=(p1+p2)/2.0

    return mid
def split6(x1l,x1r,x2l,x2r,x3l,x3r,x4l,x4r,x5l,x5r,x6l,x6r):
    lx1=(x1l+x1r)/2.0
    lx2=(x2l+x2r)/2.0
    lx3=(x3l+x3r)/2.0
    lx4=(x4l+x4r)/2.0
    lx5=(x5l+x5r)/2.0
    lx6=(x6l+x6r)/2.0

    return lx1,lx2,lx3,lx4,lx5,lx6
def split3(x1left,x1right,x2left,x2right,x3left,x3right):
    x1mid=(x1left+x1right)/2.0
    x2mid=(x2left+x2right)/2.0
    x3mid=(x3left+x3right)/2.0

    return x1mid,x2mid,x3mid
def tree_depth(node):
    dep=[]
    #print "node len:",len(node)
    for i in range(len(node)):
        dep.append(node[i].depth)
    dindex=np.argmax(dep)

    return dep[dindex]

"""This is the test function you can search this function and replace it with your own one"""
def test_function(x1,x2,x3,x4,x5,x6):
    y=(x1**2+x2*3)*np.sin(x3*x4)+10*x5+np.cos(x6)
    return y

def kernel_calculation(idtr,xtr,xte):
    ktr=np.zeros((idtr,idtr))
    ktr_new=np.matrix(ktr)
    kte=np.zeros((1,idtr))
    kte_new=np.matrix(kte)
    #kte_train=np.zeros((idtr,idtr))
    #kte_test=np.zeros((idtr_l,idtr_l))
    #kte_kte=np.zeros((idtr,idtr_l))
    ### training kernel matrix calculation
    for i in range(0,idtr):
        for j in range(i,idtr):
            #ktr_new[i,j]=exp(-0.5*0.03*LA.norm(xtr[:,i]-xtr[:,j])**2)
            ktr_new[i,j]=0.5*(1+np.sqrt((5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.05)+
            5*LA.norm(xtr[:,i]-xtr[:,j])**2/(3*0.05))*np.exp(-(np.sqrt(5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.05))
            #print ktr[i,j]
            #ktr_new[j,i]=ktr_new[i,j]

    for i in range(1):
        for j in range(0,idtr):
            #kte_new[i,j]=exp(-0.5*0.03*(xte[:,i]-xtr[:,j])**2)
            kte_new[i,j]=0.5*(1+np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.05)+
            5*LA.norm(xte[:,i]-xtr[:,j])**2/(3*0.05))*np.exp(-np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.05))

    return ktr_new,kte_new

class Node:
    def __init__(self, values=None, x1center=None, x2center=None, x3center=None, x4center=None,x5center=None,x6center=None,depth=None, x1leftmin=None, x1rightmax=None,x2leftmin=None,x2rightmax=None, x3leftmin=None, x3rightmax=None,x4leftmin=None, x4rightmax=None, x5leftmin=None, x5rightmax=None, x6leftmin=None,x6rightmax=None,parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.x1center=x1center
        self.x2center=x2center
        self.x3center=x3center
        self.x4center=x4center
        self.x5center=x5center
        self.x6center=x6center
        self.depth=depth
        self.values=values
        self.x1leftmin=x1leftmin
        self.x1rightmax=x1rightmax
        self.x2leftmin=x2leftmin
        self.x2rightmax=x2rightmax
        self.x3leftmin=x3leftmin
        self.x3rightmax=x3rightmax
        self.x4leftmin=x4leftmin
        self.x4rightmax=x4rightmax
        self.x5leftmin=x5leftmin
        self.x5rightmax=x5rightmax
        self.x6leftmin=x6leftmin
        self.x6rightmax=x6rightmax

    def Selectnode(self,node,h):
        leaves_of_depth=[]
        index=[]
        value=[]
        for i in range(len(node)):
            if node[i].depth==h:
                leaves_of_depth.append(node[i])
                index.append(i)

        for i in range(len(leaves_of_depth)):
            value.append(leaves_of_depth[i].values)
            #print value
        new_index=np.argmax(value)
        h1=node[index[new_index]].x1center
        h2=node[index[new_index]].x2center
        h3=node[index[new_index]].x3center
        h4=node[index[new_index]].x4center
        h5=node[index[new_index]].x5center
        h6=node[index[new_index]].x6center

        #print ("domain1:",rescaling(0.0,360.0,h1))
        #print ("domain2:",rescaling(0.0,360.0,h2))
        #print ("domain3:",rescaling(0.0,360.0,h3))
        #print ("domain4:",rescaling(-3.0,3.0,h4))
        #print ("domain5:",rescaling(-3.0,3.0,h5))
        #print ("domain6:",rescaling(1.0,5.5,h6))

        #print ("domain1:",h1)
        #print ("domain2:",h2)
        #print ("domain3:",h3)
        #print ("domain4:",h4)
        #print ("domain5:",h5)
        #print ("domain6:",h6)

        fh=node[index[new_index]].values
        #print "function value:",fh
        #gl=function(node[index[new_index]].center)
        #print gl


        return node[index[new_index]],index[new_index],fh
    def Addnode(self,values,x1center,x2center,x3center, x4center,x5center,x6center,depth,x1leftmin,x1rightmax,x2leftmin,x2rightmax,x3leftmin,x3rightmax,x4leftmin,x4rightmax,x5leftmin,x5rightmax,x6leftmin,x6rightmax):
        n = Node(values=values,x1center=x1center,x2center=x2center,x3center=x3center,x4center=x4center,x5center=x5center,x6center=x6center,depth=depth,x1leftmin=x1leftmin,x1rightmax=x1rightmax,
        x2leftmin=x2leftmin,x2rightmax=x2rightmax, x3leftmin=x3leftmin,x3rightmax=x3rightmax,x4leftmin=x4leftmin,x4rightmax=x4rightmax,
        x5leftmin=x5leftmin,x5rightmax=x5rightmax,x6leftmin=x6leftmin,x6rightmax=x6rightmax,parent = self)
        self.childNodes.append(n)
        return n



def SOO(center1,center2,center3,center4,center5,center6,fe,x1min,x1max,x2min,
        x2max,x3min,x3max,x4min,x4max,x5min,x5max,x6min,x6max):
    """set g(0,0)=f(x(0,0))"""
    ini_f=float("-inf")
    #g_function=get_energy_ads(1.5,1.5,2.5)
    #g_function=simplePSO().myPSO([center1,center2,center3,center4,center5,center6])
    g_function=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
    #print g_function
    """set f^=g(0,0)"""
    f_value=g_function

    """initialize the tree"""
    rootnode = Node(values=g_function,x1center=center1,x2center=center2, x3center=center3,x4center=center4,x5center=center5,x6center=center6,depth=0,x1leftmin=x1min, x1rightmax=x1max,x2leftmin=x2min,x2rightmax=x2max, x3leftmin=x3min,x3rightmax=x3max,x4leftmin=x4min,x4rightmax=x4max,x5leftmin=x5min,x5rightmax=x5max,x6leftmin=x6min,x6rightmax=x6max)


    """set t=1,n=1,N=1 and D_t={(x00,g00)}"""
    t=1 ## represent the number of sampled data
    #n=1 ## represent the iteration of the algorithm
    N=1 ## represent the number of GP calculation
    #D_x_train=np.array([center1,center2,center3,center4,center5,center6])
    #D_y_train=np.array([g_function])

    current_node=[]
    current_node.append(rootnode)
    #print current_node[0].values
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]
    final_result=[]
    axis=[]
    best_axis=[]
    t=1
    h_tree=0
    #ini_f=float("-inf")
    f_eva=0
    while f_eva<=10000:
        #print n
        v_max=float("-inf")
        h_max=np.sqrt(t)
        h_tree=tree_depth(current_node)
        #print "h_tree",h_tree
        loop=min(h_max,h_tree)
        #print "min:",loop
        h=0
        #print "iteration:",t

        for h in range(int(loop)+1):
            #print "tree depth:",h
            #print len(current_node)
            #print "tree depth:", h
            #print len(current_node)
            check_leaves=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            #print (check_leaves)
            if h in check_leaves:
                node,index,f_value = node.Selectnode(current_node,h)
                current_node.pop(index)
                #print "g_value and v_max", g_value,v_max
                if f_value>=v_max:
                    #print "cao"
                    center=[]
                    most_width=[]
                    widthx1=node.x1rightmax-node.x1leftmin
                    most_width.append(widthx1)
                    widthx2=node.x2rightmax-node.x2leftmin
                    most_width.append(widthx2)
                    widthx3=node.x3rightmax-node.x3leftmin
                    #print widthx1,widthx2,widthx3
                    most_width.append(widthx3)
                    widthx4=node.x4rightmax-node.x4leftmin
                    most_width.append(widthx4)
                    #print "node",node.x5rightmax,node.x5leftmin
                    widthx5=node.x5rightmax-node.x5leftmin
                    #print "widthx5",widthx5
                    most_width.append(widthx5)
                    widthx6=node.x6rightmax-node.x6leftmin
                    most_width.append(widthx6)
                    width_index=np.argmax(most_width)
                    #print ("index",width_index)
                    if width_index==5:
                        ele=split8(node.x6leftmin,node.x6rightmax)
                        e2=node.x6leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,e2)

                        centerx21,centerx22,centerx23,centerx24,centerx25,centerx26=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,e2,e3)

                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,e3,node.x6rightmax)


                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centerx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)

                        for i in range(3):
                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])

                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                                #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,e2))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,e2,e3))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,e3,node.x6rightmax))

                    if width_index==4:
                        ele=split8(node.x5leftmin,node.x5rightmax)
                        e2=node.x5leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,e2,node.x6leftmin,node.x6rightmax)

                        centerx21,centerx22,centerx23,centerx24,centerx25,centerx26=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,e2,e3,node.x6leftmin,node.x6rightmax)

                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,e3,node.x5rightmax,node.x6leftmin,node.x6rightmax)


                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centerx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)

                        for i in range(3):
                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])

                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                            #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,e2,node.x6leftmin,node.x6rightmax))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,e2,e3,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,e3,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                    if width_index==3:
                        ele=split8(node.x4leftmin,node.x4rightmax)
                        e2=node.x4leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,e2,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)

                        centerx21,centerx22,centerx23,centerx24,centerx25,centerx26=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,e2,e3,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)

                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,e3,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)


                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centerx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)

                        for i in range(3):
                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])

                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                            #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                node.x4leftmin,e2,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                e2,e3,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,
                                e3,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                    if width_index==2:
                        ele=split8(node.x3leftmin,node.x3rightmax)
                        e2=node.x3leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)

                        centerx21,centerx22,centerx23,centerx24,centerx25,centerx26=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)

                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)


                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centerx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)

                        for i in range(3):

                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])

                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                            #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax,
                                node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                    if width_index==1:
                        ele=split8(node.x2leftmin,node.x2rightmax)
                        e2=node.x2leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)
                        centerx21,centerx22,centerx23,centrx24,centerx25,centerx26=split6(node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)
                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)


                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centrx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)
                        #node_function_value=[]

                        for i in range(3):
                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])

                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                            #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax
                                ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax
                                ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax
                                ,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                    if width_index==0:
                        ele=split8(node.x1leftmin,node.x1rightmax)
                        e2=node.x1leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13,centerx14,centerx15,centerx16=split6(node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)
                        centerx21,centerx22,centerx23,centerx24,centerx25,centerx26=split6(e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)
                        centerx31,centerx32,centerx33,centerx34,centerx35,centerx36=split6(e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax)

                        leftnode=[centerx11,centerx12,centerx13,centerx14,centerx15,centerx16]
                        midnode=[centerx21,centerx22,centerx23,centerx24,centerx25,centerx26]
                        rightnode=[centerx31,centerx32,centerx33,centerx34,centerx35,centerx36]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)
                        #node_function_value=[]

                        for i in range(3):
                            value=simplePSO().myPSO([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            #value=simplePSO().myPSO([rescaling(0.0,360.0,center1),rescaling(0.0,360.0,center2),rescaling(0.0,360.0,center3),
                            #rescaling(-3.0,3.0,center4),rescaling(-3.0,3.0,center5),rescaling(1.0,5.5,center6)])
                            #value=simplePSO().myPSO([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            #value=get_energy_ads(rescaling(0.0,1.0,center[i][0]),rescaling(0.0,4.0,center[i][1]),rescaling(1.2,4.2,center[i][2]),
                            #rescaling(0.0,90.0,center[i][3]),rescaling(0.0,90.0,center[i][4]),rescaling(0.0,360.0,center[i][5]))
                            #print ([center[i][0],center[i][1],center[i][2],center[i][3],center[i][4],center[i][5]])
                            print ([rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])])
                            axis.append([[rescaling(0.0,360.0,center[i][0]),rescaling(0.0,360.0,center[i][1]),rescaling(0.0,360.0,center[i][2]),
                            rescaling(-3.0,3.0,center[i][3]),rescaling(-3.0,3.0,center[i][4]),rescaling(1.0,5.5,center[i][5])]])
                            function_evalution.append(value)
                            kl=np.argmax(function_evalution)
                            final_result.append(function_evalution[kl])
                            best_axis.append(axis[kl])
                            #print "centers:",center[i][0],center[i][1],center[i][2]
                            print ("function value:",value)
                            f_eva=f_eva+1
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(value,center[0][0],center[0][1],center[0][2],center[0][3],center[0][4],center[0][5],node.depth+1,node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                            if i==1:
                                current_node.append(node.Addnode(value,center[1][0],center[1][1],center[1][2],center[1][3],center[1][4],center[1][5],node.depth+1,e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            if i==2:
                                current_node.append(node.Addnode(value,center[2][0],center[2][1],center[2][2],center[2][3],center[2][4],center[2][5],node.depth+1,e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax,node.x4leftmin,node.x4rightmax,node.x5leftmin,node.x5rightmax,node.x6leftmin,node.x6rightmax))
                                #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))




                    t=t+1
                    v_max=f_value
                    #print (f_eva)

            if f_eva>=fe:
                break

        if f_eva>=fe:
            break

    return final_result,best_axis




if __name__ == "__main__":
    simplePSO().main()
    #simplePSO().load_post(args.init)
    """the highper-parameters of current version is set by hand"""



    """6 domain [0,1] [0,4] [1.2,4.2] [0,90] [0,90] [0,360]"""
    """center of each domain value"""
    x1=0.5
    x2=0.5
    x3=0.5
    x4=0.5
    x5=0.5
    x6=0.5
    """search range of each domain"""
    x1min=0.0
    x1max=1.0
    x2min=0.0
    x2max=1.0
    x3min=0.0
    x3max=1.0
    x4min=0.0
    x4max=1.0
    x5min=0.0
    x5max=1.0
    x6min=0.0
    x6max=1.0
    """number of function evaluation"""
    fe=10000
    d,b=SOO(x1,x2,x3,x4,x5,x6,fe,x1min,x1max,x2min,x2max,x3min,x3max,x4min,x4max,x5min,x5max
            ,x6min,x6max)
    #print (d)
    #print (b)
    with open('energy.txt', 'w') as f:
        for s in d:
            f.write(str(s) + '\n')
    with open('axis.txt', 'w') as f:
        for a in b:
            f.write(str(a) + '\n')
    #simplePSO().myPSO()
