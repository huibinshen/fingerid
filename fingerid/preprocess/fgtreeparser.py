"""
=============================
Parser for fragmentation tree
=============================

Parse fragmenation tree (.dot file) to FragTree object defined in fgtree.py.

"""
import re
import sys
import commands
import os
import numpy

from fgtree import FragTree

class FragTreeParser:
    
    def parse_dir(self, dir_path="NULL"):
        """
        Parse directory of dot tree file into list of FragTree instance   

        """
        tree_list = []
        dir_path = os.path.abspath(dir_path)
        # if user's path is not having a "/"                                   
        if dir_path[-1] != "/":
            dir_path = dir_path + "/"
        # invoke parse file for every file in the dir_path directory           
        files = commands.getoutput("ls %s" % dir_path).split()
        count = 0
        for f in files:
            tree = self.parse_file(dir_path + f)
            tree_list.append(tree)
            count = count + 1
        return tree_list

    def parse_file(self, f_path="NULL"):
        """ parse file into FragTree instance """
        if f_path == "NULL":
            raise Exception("ERROR: please specify fgtree file path")
        return self._parse_dot_file(f_path)

    def _parse_dot_file(self, f_path):
#        print "Parse dot file:", f_path
        fgtree = FragTree()
        fgtree.path = f_path
        nodes = {}
        edges = []

        data = open(f_path).read()
        for line in data.split('\n'):
            if not line or line.find("{")!=-1 or line.find("}")!=-1:
                continue
            if line.find('->') == -1: # is a node
                node = FragTree.Node()
                node.id = line[:line.find(' ')]
                node.label = line[line.find("\"")+1:line.find("\\n")]
                node.mass = float(line[line.find("\\n")+2:line.find(" Da")])
                node.inten = float(line[line.find(",")+2:line.find(" %")])
                nodes[node.id] = node
                if numpy.isinf(node.inten):
                    node.inten = 100
            else: # is a edge
                edge = FragTree.Edge()
                s = line.find("\"")
                t = line.find("\"",s+1)
                edge.label = line[s+1:t]
                i = line.find(" -> ")
                j = line.find(" ",i+4)
                edge.A = line[:i]
                edge.B = line[i+4:j]
                edges.append(edge)

        fgtree.metlin_id = f_path[f_path.find("pos")+3:f_path.find(".")]
        fgtree.nodes = nodes
        fgtree.edges = edges
        return fgtree

