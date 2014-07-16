"""
=================================
The class for fragmentation tree.
=================================

"""
import re

class FragTree:
    def __init__(self):
        self.metlin_id = "" # the metlin id corresponds to the fgtree
        self.edges = [] # a list of Edge instance
        self.nodes = {} # id as key, node as value
        self.leaves = [] # a list of ids which are leaves
        self.root = None # a node instance
        self.path = "" # the file location for the fgtree
        

    class Node:
        def __init__(self):
            self.id = ""
            self.label = "" # chemical formula
            self.mass = 0
            self.inten = 0
            self.children = [] # a list of Node instance
            self.par = None # the parent node
        
        def __str__(self):
            return self.label + ' ' + str(self.mass) + ' ' + str(self.inten)
    
    class Edge:
        def __init__(self):
            self.label = "" # chemical formula loss
            self.A = "" # node id
            self.B = "" # node id
            
        def __str__(self):
            return self.A + '->' + self.B + ' ' + self.label

    def link_nodes(self):
        """ add children and parent attibuties to every node"""
        for e in self.edges:
            n1 = self.nodes[e.A]
            n2 = self.nodes[e.B]
            if n2 in n1.children:
                continue
            n1.children.append(n2)
            n2.par = n1

    def get_leaves_and_root(self):
        """ set leaves and root of the tree """
        for nid, node in self.nodes.items():
            if len(node.children) == 0:
                self.leaves.append(nid)
            if node.par is None:
                self.root = node

    def get_post_order(self):
        """return a list of nodes id in post order traversal """
        root = self.root
        queue = []
        self._post_order(root, queue)
        return queue

    def _post_order(self,node, queue):
        for c in node.children:
            self._post_order(c,queue)
        queue.append(node.id)

    def get_edges_dict(self):
        """ Return a dict with v1->v2 as keys and edge loss as values"""
        e_dict = {}
        for e in self.edges:
            e_dict[e.A+'->'+e.B] = e.label
        return e_dict

    def get_root_loss(self):
        """ return a list of losses from the root, losses are mergered"""
        self.link_nodes()
        self.get_leaves_and_root()
        order = self.get_post_order()
        edges = self.get_edges_dict()
        root_loss = []
        for nid in order:
            node = self.nodes[nid]
            t_list = []
            while node.par is not None:
                par = node.par
                t_loss = edges[par.id+'->'+node.id]
                t_list.append(t_loss)
                node = par
            # merge chemical formula in t_list 
            atoms_dict = {'C':0,'H':0,'O':0,'N':0,'S':0,'P':0,'I':0,'Br':0,'F':0,'Cl':0}
            for loss in t_list:
                l = re.findall(r'([A-Z][a-z]*)(\d*)', loss)
                for a,c in l:
                    if not c:
                        c = 1
                    if a not in atoms_dict:
                        continue
                    atoms_dict[a] = atoms_dict[a] + int(c)
            m_loss = 'C'+str(atoms_dict['C'])+'H'+str(atoms_dict['H'])+'O'+str(atoms_dict['O'])+'N'+str(atoms_dict['N'])+'P'+str(atoms_dict['P'])+'S'+str(atoms_dict['S'])+'I'+str(atoms_dict['I'])+'F'+str(atoms_dict['F'])+'Cl'+str(atoms_dict['Cl'])+'Br'+str(atoms_dict['Br'])
            root_loss.append((nid,m_loss))
        self.root_loss = root_loss
        return root_loss


                
                
