#! usr/bin/env python

# File: traithull.py
# Author: Dylan Schwilk
# Copyright 2003 Dylan W. Schwilk
# www.pricklysoft.org

####################################################################
# GNU
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version.
#   
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
####################################################################

"""
Traithull.py -

Provides an interface to the Qhull program that allows ecologists to easily
calculate the convex hull volume (CHV) metric of functional diversity and do
tests against null models. For a description of the method see: Cornwell, W.K,
D.W. Schwilk and D.D. Ackerly. 2006. A trait-based test for habitat filtering:
convex hull volume. Ecology 87: 1465--1471.
"""

__version__ = "1.4"
__author__ = '''Dylan Schwilk'''
__needs__ = '2.3'

import os.path
import random
import math
import stats # requires the included GNU stats module by Gary Strangman
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
traithull_logger = logging.getLogger('bibstuff_logger')

# globals
qhull_path = ""

def G2(vals) :
    """Calculates the adjusted and corrected kurtosis (G2) as
    defined on pg 115 of Sokal and Rohlf"""
    n = len(vals)
    if n < 4 : return -99  # hack
    return ((n+1)*n*stats.kurtosis(vals)) / ((n-1) * (n-2)) - ((3*(n-1)*(n-1))/((n-2)*(n-3)))

def TraitMatrix(lines):
    '''reads trait matrix. Returns tuple (dict, list of trait names)'''
    result = {}
    #print trait_names
    for line in lines[1:] : # skip header
        species =  line.split()
        result[species[0]] = map(float,species[1:])
    return (result, lines[0].split()[1:])

def PlotMatrix(lines):
    '''Reads in species x plot data'''
    result = {}
    plot_names = lines[0].split()[1:]
    for line in lines[1:] :
        species = line.split()
        species_name = species[0]
        for plot in range(len(plot_names)) :
            if int(species[plot+1]) :
                if result.has_key(plot_names[plot]) :
                    result[plot_names[plot]].append(species_name)
                else :
                     result[plot_names[plot]] = [species_name,]
    #print "\nplot: ",result
    return result

def HullVolume(species_list, trait_matrix):
    '''Use qhull to produce volume measure'''

    # if single dimension, return range
    dim = len(trait_matrix[species_list[0]])
    if dim < 2 :
        vals = map(lambda x: trait_matrix[x][0], species_list)
        return max(vals) - min(vals)
          
    num_points = len(species_list)

    if num_points < dim + 1 :
        return 0
       
    instring = "%d %d " % (dim , num_points )
    for name in species_list :
        for trait in trait_matrix[name] :
            instring += " %f" % trait
        instring += '\n'

    input, output =  os.popen4(os.path.join(qhull_path, "qconvex") + " FA")
    input.write(instring)
    input.close()
    
    data = output.read()
    traithull_logger.debug( "Total volume: %s" % data[data.find("Total volume:") + 14 :])
    ind = data.find("volume:")
    if ind == -1 :
        print data
        return 0
    return float(data[ind + 7 :] )
     

def IndividualTraitStats(species_list, trait_matrix, trait_names):
    '''Get variance and kurtosis for each trait. returns dict with
    keys trait_names and data a tuple of range, var, kurtosis (G2),
    nndist_mean, nndist_var
    '''
    result = {}
    for trait in trait_names :
        temp_matrix = {}
        vals = []
        for spec in species_list :
            val = trait_matrix[spec][trait_names.index(trait)]
            vals.append(val)
            temp_matrix[spec] = [val]
        distances = NearestNeighborDistances(species_list, temp_matrix)
        range = max(vals) - min(vals)
        result[trait] = (range, stats.var(vals), G2(vals),
                         stats.mean(distances), stats.var(distances))  
    return result
                        

def EuclideanDistance(species_a, species_b, trait_matrix):
    "Return euclidian distance between a and b in trait space"
    return math.sqrt(EuclideanDistanceSquare(species_a, species_b, trait_matrix))


def EuclideanDistanceSquare(species_a, species_b, trait_matrix):
  "Return square of euclidian distance between a and b in trait space"
  sum = 0
  dim = len(trait_matrix[species_a])
  for trait in range(0,dim):
      v = trait_matrix[species_a][trait] - trait_matrix[species_b][trait]
      sum += v*v
  return sum

def NeighborDistances(species_list, trait_matrix, dist_fun = EuclideanDistance ):
    "return list of all pairwise distances.  This is a list of lists"
    result = []
    for i in species_list :
        l = []
        for j in species_list:
            if i != j :
                l.append(dist_fun(i,j, trait_matrix))
        result.append(l)
    return result    

def NearestNeighborDistances(species_list, trait_matrix,
                             dist_fun = EuclideanDistance):
    "return list of nearest neighbor distances"
    return(map(min, NeighborDistances(species_list, trait_matrix, dist_fun)))

def AussieDistances(species_list, trait_matrix,
                    dist_fun = EuclideanDistanceSquare):
    "return list of nearest neighbor distances"
    return(map(sum, NeighborDistances(species_list, trait_matrix, dist_fun)))


def PrintResultRow(rowname, species_set, trait_matrix, trait_names,
                   options, include_indiv=0) :
    "Print one row of output"
    print "%s\t%d\t%f" % (rowname, len(species_set), HullVolume(species_set, trait_matrix)),

    #individual trait stats     
    if include_indiv :
        trait_stats = IndividualTraitStats(species_set, trait_matrix, trait_names)
        for trait in trait_names :
            print "\t%f\t%f\t%f\t%f\t%f" % trait_stats[trait],

    #Nearest-neighbor distances        
    if options.do_dist :
        distances = NearestNeighborDistances(species_set, trait_matrix)
        print "\t%f\t%f" % (stats.mean(distances), stats.var(distances)),

    if options.do_aussie :
        distances = AussieDistances(species_set, trait_matrix,EuclideanDistanceSquare)
        print "\t%f" % sum(distances), 

    # finish row
    print '\n',    

        
def main():
    '''Command line version of tool'''
    
    from optparse import OptionParser
    import sys
    global qhull_path
    
    usage = "usage: %prog [options] [trait_file]"
    parser = OptionParser(usage=usage, version ="%prog " + __version__)
    parser.add_option("-p", "--plotfile", action="store", type="string", \
                      dest="plotfile", default = "", help="File containing species occurance by plot")
    parser.add_option("-q", "--qhull", action="store", type="string", \
                      dest="qpath", default = qhull_path, help="path to qconvex executable")
    parser.add_option("-r", "--randsample", action="store", type="int", \
                      dest="replicates", default = 1, help="Number of random samples per richness")
    parser.add_option("-d", "--distance", action = "store_true", dest="do_dist", default = 0, \
                      help = "Output mean and variance of nearest-neighbor distances")
    parser.add_option("-a", "--Aussie", action = "store_true", dest="do_aussie", default = 0, \
                      help = "Output Aussie (Walker) fun. div. index")
    parser.add_option("-i", "--individual", action = "store_true", dest="do_indiv", default = 0, \
                      help = "Do each treat in matrix individually (1-dimensional version)")
    parser.add_option("-t", "--total", action="store_true", \
                      dest="do_total", default=0, help="Output total species pool results")
    parser.add_option("-v", "--verbose", action="store_true", \
                      dest="verbose", default=0, help="Verbose output")

    # get options
    (options, args) = parser.parse_args()
    if options.verbose:
	traithull_logger.setLevel(logging.INFO)
    qhull_path = options.qpath

    if len(args) > 0 :
        try :
           species_file =  open(args[0]).readlines()
        except:
            traithull_logger.error("Bad or missing species input file: %s" % args[0])
            sys.exit(1)
    else :
        species_file = sys.stdin.readlines()

    tMatrix, trait_names = TraitMatrix(species_file)
  
    print "Species_set\tRichness\tVolume",
    if options.do_indiv :
        for trait in trait_names :
            print "\t%s_range\t%s_var\t%s_kurtosis\t%s_nnmean\t%s_nnvar" % (trait, trait,trait,trait,trait),
    if options.do_dist :
        print "\tMeanNNeighbor\tNeighborVar"
    if options.do_aussie :
        print "\tAussie"
    print '\n',
   
    
    if options.plotfile :
        try:
            pMatrix = PlotMatrix(open(options.plotfile).readlines())
        except:
            traithull_logger.error("Bad or missing plot input file: %s" % options.plotfile)
            sys.exit(1)
        for plot in pMatrix.keys() :
           PrintResultRow(plot, pMatrix[plot], tMatrix, trait_names, options, options.do_indiv)
    elif options.replicates :
        dim = len(tMatrix[tMatrix.keys()[0]])    
        for r in range(dim+1, len(tMatrix.keys())):
            for i in range(options.replicates):
                plot = random.sample(tMatrix.keys(), r)
                PrintResultRow("Random", plot, tMatrix, trait_names, options, options.do_indiv)
    # output totals
    if options.do_total :
        PrintResultRow("Total", tMatrix.keys(), tMatrix, options)
            
if __name__ == '__main__':
    main()
