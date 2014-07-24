"""
=====================================
Compute kernel for fragmentation tree
=====================================

"""
from mskernel import Kernel
import numpy

class FragTreeKernel(Kernel):
    
    
    def compute_train_kernel(self, trees, k_type, sm=0.00001, si=100000):
        """
        Parameters:
        -----------
        trees: list, a list of fgtree objects.
        k_type: str, kernel type, must be one of the following types:
         "NB", "NI", "LB", "LI", "LC", "RLB", "RLI", "CSC", "CPK", "CP2","CSC".
        sm: float, std in mass dimention, only used by CPK
        si: float,  std in intensit dimention only used by CPK.

        Returns:
        --------
        km, numpy 2d array, kernel matrix.
        
        """
        n = len(trees)
#        nodes_dict = {}
#        loss_dict = {}
#        edge_dict = {}
#        for t in trees:
#            edges = t.edges
#            nodes = t.nodes
#            for n_id, node in nodes.items():
#                nodes_dict[node.label] = 0
#            for e in edges:#
#                loss_dict[e.label] = 0
#            for e in edges:
#                e_name = nodes[e.A].label+'->'+nodes[e.B].label
#                edge_dict[e_name] = 0
        #print "Has %d nodes features" % len(nodes_dict)
        #print "Has %d losses features" % len(loss_dict)
        #print "Has %d edges features" % len(edge_dict)

#        nodes_list = nodes_dict.keys()
#        loss_list = loss_dict.keys()

        # takes all the nodes in the trees as features
        if k_type == 'NB':
            # binary valued feature vector
            km = numpy.zeros((n,n))
#            m = len(nodes_dict)
            for i in range(n):
                for j in range(i,n):
                    T1_nodes = [node.label for node in trees[i].nodes.values()]
                    T2_nodes = [node.label for node in trees[j].nodes.values()]
#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for node_label in nodes_list:
#                        if node_label in T1_nodes:
#                            T1_v[k] = 1
#                            ks_1.append(k)
#                        if node_label in T2_nodes:
#                            ks_2.append(k)
#                            T2_v[k] = 1
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
                    count = 0
                    for t1node in T1_nodes:
                        if t1node in T2_nodes:
                            count +=1
                    km[i,j] = count
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type=='NC':
            # count existence for every chemical formula, 
            km = numpy.zeros((n,n))
#            m = len(nodes_dict)
            for i in range(n):
                for j in range(i,n):
                    T1_labels= [node.label for node in trees[i].nodes.values()]
                    T2_labels= [node.label for node in trees[j].nodes.values()]
                    T1_dict = {}
                    T2_dict = {}
                    for l in T1_labels:
                        T1_dict[l] = T1_dict.get(l,0) + 1
                    for l in T2_labels:
                        T2_dict[l] = T2_dict.get(l,0) + 1
                    tsum = 0
                    for t1node, t1count in T1_dict.items():
                        if t1node in T2_dict:
                            tsum = tsum + t1count*T2_dict[t1node]
                    km[i,j] = tsum
                    km[j,i] = km[i,j]
                    #T1_v = numpy.zeros(m)
                    #T2_v = numpy.zeros(m)
                    #k = 0
                    #for node_label in nodes_list:
                    #    if node_label in T1_labels:
                    #        count = 0
                    #        for t1_label in T1_labels:
                    #            if t1_label == node_label:
                    #                count = count +1
                    #        T1_v[k] = count
                    #    if node_label in T2_labels:
                    #        count = 0
                    #        for t2_label in T2_labels:
                    #            if t2_label == node_label:
                    #                count = count +1
                    #        T2_v[k] = count
                    #    k = k+1
                    #km[i,j] = sum(T1_v * T2_v)
                    #km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'NI':
            # add the intensities for every chemical formula
            km = numpy.zeros((n,n))
            #m = len(nodes_dict)
            for i in range(n):
                for j in range(i,n):
                    T1_nodes = [node for node in trees[i].nodes.values()]
                    T2_nodes = [node for node in trees[j].nodes.values()]

                    T1_dict = {}
                    T2_dict = {}
                    for t1node in T1_nodes:
                        t1label = t1node.label
                        t1inten = t1node.inten
                        T1_dict[t1label] = t1inten
                    for t2node in T2_nodes:
                        t2label = t2node.label
                        t2inten = t2node.inten
                        T2_dict[t2label] = t2inten

                    tsum = 0
                    for t1node, t1inten in T1_dict.items():
                        if t1node in T2_dict:
                            tsum = tsum + t1inten*T2_dict[t1node]
                    km[i,j] = tsum
                    km[j,i] = km[i,j]

#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for node_label in nodes_list:
#                        if node_label in T1_labels:
#                            inten = 0
#                            for t1_node in T1_nodes:
#                                if t1_node.label == node_label:
#                                    inten = inten + t1_node.inten
#                            T1_v[k] = inten
#                        if node_label in T2_labels:
#                            inten = 0
#                            for t2_node in T2_nodes:
#                                if t2_node.label == node_label:
#                                    inten = inten + t2_node.inten
#                            T2_v[k] = inten
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
#                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        # take the edge loss as features
        elif k_type == 'LB':
            km = numpy.zeros((n,n))
#            m = len(loss_dict)
            for i in range(n):
                for j in range(i,n):
                    T1_edges = trees[i].edges
                    T2_edges = trees[j].edges
                    T1_loss = {}
                    T2_loss = {}
                    for e in T1_edges:
                        T1_loss[e.label] = 1
                    for e in T2_edges:
                        T2_loss[e.label] = 1
                    count = 0
                    for t1loss in T1_loss.keys():
                        if t1loss in T2_loss:
                            count += 1
                    km[i,j] = count
                    km[j,i] = km[i,j]
#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for loss_label in loss_list:
#                        if loss_label in T1_loss:
#                            T1_v[k] = 1
#                        if loss_label in T2_loss:
#                            T2_v[k] = 1
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
#                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'LC':
            km = numpy.zeros((n,n))
#            m = len(loss_dict)
            #print loss_list
            for i in range(n):
                for j in range(i,n):
                    T1_edges = trees[i].edges
                    T2_edges = trees[j].edges
                    T1_loss = {}
                    T2_loss = {}
                    for e in T1_edges:
                        T1_loss[e.label] = T1_loss.get(e.label,0) + 1
                    for e in T2_edges:
                        T2_loss[e.label] = T2_loss.get(e.label,0) + 1
                    tsum = 0
                    for t1loss, t1count in T1_loss.items():
                        if t1loss in T2_loss:
                            tsum = tsum + t1count*T2_loss[t1loss]
                    km[i,j] = tsum
                    km[j,i] = km[i,j]
#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0                    
#                    for loss_label in loss_list:
#                        count = 0
#                        for e in T1_edges:
#                            if e.label == loss_label:
#                                count = count +1
#                        T1_v[k] = count
 
#                        count = 0
#                        for e in T2_edges:
#                            if e.label == loss_label:
#                                count = count +1
#                        T2_v[k] = count
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
#                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'LI':
            km = numpy.zeros((n,n))
#            m = len(loss_dict)
            for i in range(n):
                for j in range(i,n):
                    T1_edges = trees[i].edges
                    T2_edges = trees[j].edges
                    T1_nodes = trees[i].nodes
                    T2_nodes = trees[j].nodes

                    T1_loss = {}
                    T2_loss = {}
                    for e in T1_edges:
                        T1_loss[e.label] = T1_loss.get(e.label,0) + T1_nodes[e.B].inten
                    for e in T2_edges:
                        T2_loss[e.label] = T2_loss.get(e.label,0) + T2_nodes[e.B].inten
                    tsum = 0
                    for t1loss, t1inten in T1_loss.items():
                        if t1loss in T2_loss:
                            tsum = tsum + t1inten*T2_loss[t1loss]
                    km[i,j] = tsum
                    km[j,i] = km[i,j]

#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for loss_label in loss_list:
#                        inten = 0
#                        for e in T1_edges:
#                            if e.label == loss_label:
                                #inten+= T1_nodes[e.A].inten+T1_nodes[e.B].inte
#                                inten += T1_nodes[e.B].inten
#                        T1_v[k] = inten

#                        inten = 0
#                        for e in T2_edges:
#                            if e.label == loss_label:
#                                inten += T2_nodes[e.B].inten
#                                #inten+=T2_nodes[e.A].inten+T2_nodes[e.B].inten
#                        T2_v[k] = inten
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
#                    km[j,i] = km[i,j]
#                    print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'CPC':
            km = numpy.zeros((n,n))
            # nodes id in post order for all the trees
            visit = [] 
            # list of dict with v1->v2 as key and edge label of value
            edges_dict = []
            for t in trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                visit.append(t.get_post_order())
                edges_dict.append(t.get_edges_dict())

            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]
                    T1_visit = visit[i]
                    T2_visit = visit[j]
                    edict1 = edges_dict[i]
                    edict2 = edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                    # compare two nodes in two trees

                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            score = score + (1 + D[ii1,ii2])
                                D[ind1,ind2] = score
                    km[i,j] = numpy.sum(D)
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])
        elif k_type == 'CSC':
            km = numpy.zeros((n,n))
            # nodes id in post order for all the trees
            visit = [] 
            # list of dict with v1->v2 as key and edge label of value
            edges_dict = []
            for t in trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                visit.append(t.get_post_order())
                edges_dict.append(t.get_edges_dict())

            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]
                    T1_visit = visit[i]
                    T2_visit = visit[j]
                    edict1 = edges_dict[i]
                    edict2 = edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                        
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                
                                D[ind1,ind2] = 1
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 1
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            score = score * (2 + D[ii1,ii2] -1)
                                D[ind1,ind2] = score
                    km[i,j] = numpy.sum(D)
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])
        elif k_type == 'CP2':
            km = numpy.zeros((n,n))
            # nodes id in post order for all the trees
            visit = [] 
            # list of dict with v1->v2 as key and edge label of value
            edges_dict = []
            for t in trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                visit.append(t.get_post_order())
                edges_dict.append(t.get_edges_dict())

            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]
                    T1_visit = visit[i]
                    T2_visit = visit[j]
                    edict1 = edges_dict[i]
                    edict2 = edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                        
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            # first half match
                                            chlchl1 = c1.children
                                            chlchl2 = c2.children
                                            if (not chlchl1) or (not chlchl2):
                                                # c1 or c2 don't have children
                                                pass
                                            else: 
                                                for cc1 in chlchl1:
                                                    for cc2 in chlchl2:
                                                        eee1 = c1.id+'->'+cc1.id
                                                        eee2 = c2.id+'->'+cc2.id
                                                        if edict1[eee1] == edict2[eee2]:# id1->c1->cc1 == id2->c2->cc2
                                                            ii1 = int(cc1.id[1:])-1
                                                            ii2 = int(cc2.id[1:])-1
                                                            score = score + (1 + D[ii1,ii2])
                                                            
                                D[ind1,ind2] = score
                    km[i,j] = numpy.sum(D)
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'CPK':
            """ Instead of counting common path, the value in the DP table
                is the gaussian product of two comparing nodes"""
            km = numpy.zeros((n,n))
            # nodes id in post order for all the trees
            visit = [] 
            # list of dict with v1->v2 as key and edge label of value
            edges_dict = []
            for t in trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                visit.append(t.get_post_order())
                edges_dict.append(t.get_edges_dict())

            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]
                    T1_visit = visit[i]
                    T2_visit = visit[j]
                    edict1 = edges_dict[i]
                    edict2 = edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                    # compare two nodes in two trees
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            s = self._gaussproduct(c1.mass,c2.mass,c1.inten,c2.inten,sm,si)
                                            score = score + (s + D[ii1,ii2])
                                D[ind1,ind2] = score
                    km[i,j] = numpy.sum(D)
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'RLB':
            km = numpy.zeros((n,n))
#            all_r_loss = {}
            for t in trees:
                r_loss = t.get_root_loss()
#                for nid, loss in r_loss:
#                    all_r_loss[loss] = 0
            #print len(all_r_loss)
#            all_r_loss = all_r_loss.keys()

            # binary valued feature vector
            km = numpy.zeros((n,n))
#            m = len(all_r_loss)

            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]

                    T1_rloss = T1.root_loss
                    T2_rloss = T2.root_loss

                    T1_dict = {}
                    T2_dict = {}

                    for nid, loss in T1_rloss:
                        T1_dict[loss] = T1.nodes[nid].inten
                    for nid, loss in T2_rloss:
                        T2_dict[loss] = T2.nodes[nid].inten

#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for rloss in all_r_loss:
#                        for nid,loss in T1_rloss:
#                            if loss == rloss:
#                                T1_v[k] = 1
#                        for nid,loss in T2_rloss:
#                            if loss == rloss:
#                                T2_v[k] = 1
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
                    count = 0
                    for t1loss,it1nten in T1_dict.items():
                        if t1loss in T2_dict:
                            count +=1
                    km[i,j] = count
                    km[j,i] = km[i,j]
                    #print "Computing kernel for %d, %d" % (i,j)

        elif k_type == 'RLI':
            km = numpy.zeros((n,n))
#            all_r_loss = {}
            for t in trees:
                r_loss = t.get_root_loss()
#                for nid, loss in r_loss:
#                    all_r_loss[loss] = 0
#            #print len(all_r_loss)
#            all_r_loss = all_r_loss.keys()

            # binary valued feature vector
            km = numpy.zeros((n,n))
#            m = len(all_r_loss)
            for i in range(n):
                for j in range(i,n):
                    T1 = trees[i]
                    T2 = trees[j]

                    T1_rloss = T1.root_loss
                    T2_rloss = T2.root_loss

                    T1_dict = {}
                    T2_dict = {}

                    for nid, loss in T1_rloss:
                        T1_dict[loss] = T1.nodes[nid].inten
                    for nid, loss in T2_rloss:
                        T2_dict[loss] = T2.nodes[nid].inten

                    tsum = 0
                    for t1loss, t1inten in T1_dict.items():
                        if t1loss in T2_dict:
                            tsum = tsum + t1inten*T2_dict[t1loss]
                    km[i,j] = tsum
                    km[j,i] = km[i,j]
#                    T1_v = numpy.zeros(m)
#                    T2_v = numpy.zeros(m)
#                    k = 0
#                    for rloss in all_r_loss:
#                        for nid,loss in T1_rloss:
#                            if loss == rloss:
#                                T1_v[k] = T1.nodes[nid].inten
#                        for nid,loss in T2_rloss:
#                            if loss == rloss:
#                                T2_v[k] = T2.nodes[nid].inten
#                        k = k+1
#                    km[i,j] = sum(T1_v * T2_v)
#                    km[j,i] = km[i,j]
#                    print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        else:
            raise Exception("Error in tree kernel type")
        return self._normalize_km(km)


    def compute_test_kernel(self, test_trees, train_trees, k_type, sm=0.00001, si=10000):
        """
        Parameters:
        -----------
        test_trees: list, a list of test fgtree objects.
        train_trees: list, a list of train fgtree object
        k_type: str, kernel type, must be one of the following types:
          "NB", "NI", "LB", "LI", "LC", "RLB", "RLI", "CSC", "CPK", "CP2","CSC".       sm: float, std in mass dimention, only used by CPK
        si; float,  std in intensit dimention only used by CPK

        Returns:
        --------
        km, numpy 2d array (n_test * n_train), test kernel matrix.
        
        """
        n_test = len(test_trees)
        n_train = len(train_trees)

        if k_type == 'NB':
            km = numpy.zeros((n_test,n_train))
            for i in range(n_test):
                T1_nodes = [node.label for node in test_trees[i].nodes.values()]
                for j in range(n_train):
                    T2_nodes = [node.label for node in train_trees[j].nodes.values()]
                    count = 0
                    for t1node in T1_nodes:
                        if t1node in T2_nodes:
                            count +=1
                    km[i,j] = count / numpy.sqrt(len(T1_nodes)*len(T2_nodes))
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'NI':
            km = numpy.zeros((n_test,n_train))
            for i in range(n_test):
                T1_nodes = [node for node in test_trees[i].nodes.values()]
                T1_dict = {}
                t1sum = 0
                for t1node in T1_nodes:
                    t1label = t1node.label
                    t1inten = t1node.inten
                    T1_dict[t1label] = t1inten                
                    t1sum = t1sum + t1inten*t1inten

                for j in range(n_train):
                    T2_nodes = [node for node in train_trees[j].nodes.values()]
                    T2_dict = {}
                    t2sum = 0
                    for t2node in T2_nodes:
                        t2label = t2node.label
                        t2inten = t2node.inten
                        T2_dict[t2label] = t2inten
                        t2sum = t2sum + t2inten*t2inten
                    tsum = 0
                    for t1node, t1inten in T1_dict.items():
                        if t1node in T2_dict:
                            tsum = tsum + t1inten*T2_dict[t1node]
                    km[i,j] = tsum/numpy.sqrt(t1sum*t2sum)

        elif k_type == 'LB':
            km = numpy.zeros((n_test,n_train))
            for i in range(n_test):
                T1_edges = test_trees[i].edges
                T1_loss = {}
                for e in T1_edges:
                    T1_loss[e.label] = 1
                for j in range(n_train):
                    T2_edges = train_trees[j].edges
                    T2_loss = {}
                    for e in T2_edges:
                        T2_loss[e.label] = 1
                    count = 0
                    for t1loss in T1_loss.keys():
                        if t1loss in T2_loss:
                            count += 1
                    if len(T1_loss) == 0 or len(T2_loss) == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = count/numpy.sqrt(len(T1_loss)*len(T2_loss))
#                    print i,j
#                    print count, len(T1_edges), len(T2_edges)
#                    raw_input()

        elif k_type == 'LC':
            km = numpy.zeros((n_test,n_train))
            for i in range(n_test):
                T1_edges = test_trees[i].edges
                T1_loss = {}
                for e in T1_edges:
                    T1_loss[e.label] = T1_loss.get(e.label,0) + 1
                t1sum = 0
                for t1loss, t1count in T1_loss.items():
                    t1sum = t1sum+ t1count*t1count

                for j in range(n_train):
                    T2_edges = train_trees[j].edges
                    T2_loss = {}
                    for e in T2_edges:
                        T2_loss[e.label] = T2_loss.get(e.label,0) + 1
                    tsum = 0
                    for t1loss, t1count in T1_loss.items():
                        if t1loss in T2_loss:
                            tsum = tsum + t1count*T2_loss[t1loss]
                    t2sum = 0
                    for t2loss, t2count in T2_loss.items():
                        t2sum = t2sum + t2count*t2count

                    if len(T1_loss) == 0 or len(T2_loss) == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = tsum/numpy.sqrt(t1sum*t2sum)

        elif k_type == 'LI':
            km = numpy.zeros((n_test,n_train))
            for i in range(n_test):
                T1_edges = test_trees[i].edges
                T1_loss = {}
                T1_nodes = test_trees[i].nodes
                for e in T1_edges:
                    T1_loss[e.label] = T1_loss.get(e.label,0) + T1_nodes[e.B].inten
                t1sum = 0
                for t1loss, t1inten in T1_loss.items():
                    t1sum = t1sum+ t1inten*t1inten

                for j in range(n_train):
                    T2_edges = train_trees[j].edges
                    T2_nodes = train_trees[j].nodes                
                    T2_loss = {}
                    for e in T2_edges:
                        T2_loss[e.label] = T2_loss.get(e.label,0) + T2_nodes[e.B].inten
                    tsum = 0
                    for t1loss, t1inten in T1_loss.items():
                        if t1loss in T2_loss:
                            tsum = tsum + t1inten*T2_loss[t1loss]
                    t2sum = 0
                    for t2loss, t2inten in T2_loss.items():
                        t2sum = t2sum+ t2inten*t2inten
                    if len(T1_loss) == 0 or len(T2_loss) == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = tsum/numpy.sqrt(t1sum*t2sum)

        elif k_type == 'CPC':
            km = numpy.zeros((n_test,n_train))
            test_visit = [] 
            test_edges_dict = []
            train_visit = [] 
            train_edges_dict = []
            for t in test_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                test_visit.append(t.get_post_order())
                test_edges_dict.append(t.get_edges_dict())
            for t in train_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                train_visit.append(t.get_post_order())
                train_edges_dict.append(t.get_edges_dict())
            for i in range(n_test):
                for j in range(n_train):
                    T1 = test_trees[i]
                    T2 = train_trees[j]
                    T1_visit = test_visit[i]
                    T2_visit = train_visit[j]
                    edict1 = test_edges_dict[i]
                    edict2 = train_edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                    D1 = numpy.zeros(n1)-1
                    D2 = numpy.zeros(n2)-1
                    # compare two nodes in two trees
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            score = score + (1 + D[ii1,ii2])
                                D[ind1,ind2] = score
                    # compute normalizer 
                    for id1 in T1_visit:
                        ind1 = int(id1[1:])-1
                        if id1 in T1.leaves:
                            D1[ind1] = 0
                        else:
                            chl1 = T1.nodes[id1].children
                            score = 0
                            for c1 in chl1:
                                ii1 = int(c1.id[1:])-1
                                score = score + (1 + D1[ii1])
                    for id2 in T2_visit:
                        ind2 = int(id2[1:])-1
                        if id2 in T2.leaves:
                            D2[ind2] = 0
                        else:
                            chl2 = T2.nodes[id2].children
                            score = 0
                            for c2 in chl2:
                                ii2 = int(c2.id[1:])-1
                                score = score + (1 + D2[ii2])
                    sum1 = numpy.sum(D1)
                    sum2 = numpy.sum(D2)
                    if sum1 == 0 or sum2 == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = numpy.sum(D)/numpy.sqrt(sum1*sum2)
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])
        elif k_type == 'CSC':
            km = numpy.zeros((n_test,n_train))
            test_visit = [] 
            test_edges_dict = []
            train_visit = [] 
            train_edges_dict = []
            for t in test_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                test_visit.append(t.get_post_order())
                test_edges_dict.append(t.get_edges_dict())
            for t in train_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                train_visit.append(t.get_post_order())
                train_edges_dict.append(t.get_edges_dict())
            for i in range(n_test):
                for j in range(n_train):
                    T1 = test_trees[i]
                    T2 = train_trees[j]
                    T1_visit = test_visit[i]
                    T2_visit = train_visit[j]
                    edict1 = test_edges_dict[i]
                    edict2 = train_edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2)) -1
                    D1 = numpy.zeros((n1,n1))  -1
                    D2 = numpy.zeros((n2,n2)) -1
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 1
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 1
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            score = score * (2 + D[ii1,ii2] -1)
                                D[ind1,ind2] = score
                    # compute normalizer
                    for id1 in T1_visit:
                        for id2 in T1_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T1.leaves:
                                D1[ind1,ind2] = 1
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T1.nodes[id2].children
                                score = 1
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict1[ee2]:
                                            score = score * (2 +D1[ii1,ii2] -1)
                                D1[ind1,ind2] = score

                    for id1 in T2_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T2.leaves or id2 in T2.leaves:
                                D2[ind1,ind2] = 1
                            else:
                                chl1 = T2.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 1
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict2[ee1] == edict2[ee2]:
                                            score = score * (2 +D2[ii1,ii2] -1)
                                D2[ind1,ind2] = score
                    sum1 = numpy.sum(D1)
                    sum2 = numpy.sum(D2)
                    if sum1 == 0 or sum2 == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = numpy.sum(D)/numpy.sqrt(sum1*sum2)
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])
        elif k_type == 'CP2':
            km = numpy.zeros((n_test,n_train))
            test_visit = [] 
            test_edges_dict = []
            train_visit = [] 
            train_edges_dict = []
            for t in test_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                test_visit.append(t.get_post_order())
                test_edges_dict.append(t.get_edges_dict())
            for t in train_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                train_visit.append(t.get_post_order())
                train_edges_dict.append(t.get_edges_dict())
            for i in range(n_test):
                for j in range(n_train):
                    T1 = test_trees[i]
                    T2 = train_trees[j]
                    T1_visit = test_visit[i]
                    T2_visit = train_visit[j]
                    edict1 = test_edges_dict[i]
                    edict2 = train_edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                    D1 = numpy.zeros(n1)-1
                    D2 = numpy.zeros(n2)-1
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            # first half match
                                            chlchl1 = c1.children
                                            chlchl2 = c2.children
                                            if (not chlchl1) or (not chlchl2):
                                                # c1 or c2 don't have children
                                                pass
                                            else: 
                                                for cc1 in chlchl1:
                                                    for cc2 in chlchl2:
                                                        eee1 = c1.id+'->'+cc1.id
                                                        eee2 = c2.id+'->'+cc2.id
                                                        if edict1[eee1] == edict2[eee2]:# id1->c1->cc1 == id2->c2->cc2
                                                            ii1 = int(cc1.id[1:])-1
                                                            ii2 = int(cc2.id[1:])-1
                                                            score = score + (1 + D[ii1,ii2])
                                                            
                                D[ind1,ind2] = score

                    for id1 in T1_visit:
                        ind1 = int(id1[1:])-1
                        if id1 in T1.leaves:
                            D1[ind1] = 0
                        else:
                            chl1 = T1.nodes[id1].children
                            score = 0
                            for c1 in chl1:
                                chlch1 = c1.children
                                if not chlch1:
                                    pass
                                else:
                                    for cc1 in chlch1:
                                        ii1 = int(cc1.id[1:])-1
                                        score = score + (1 + D1[ii1])
                    for id2 in T2_visit:
                        ind2 = int(id2[1:])-1
                        if id2 in T2.leaves:
                            D2[ind2] = 0
                        else:
                            chl2 = T2.nodes[id2].children
                            score = 0
                            for c2 in chl2:
                                chlch2 = c2.children
                                if not chlch2:
                                    pass
                                else:
                                    for cc2 in chlch2:
                                        ii2 = int(cc2.id[1:])-1
                                        score = score + (1 + D2[ii2])
                    sum1 = numpy.sum(D1)
                    sum2 = numpy.sum(D2)
                    if sum1 == 0 or sum2 == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = numpy.sum(D)/numpy.sqrt(sum1*sum2)
                    #print "Computing kernel for %d, %d (%s)" % (i,j,km[i,j])

        elif k_type == 'CPK':
            """ Instead of counting common path, the value in the DP table
                is the gaussian product of two comparing nodes"""
            km = numpy.zeros((n_test,n_train))
            test_visit = [] 
            test_edges_dict = []
            train_visit = [] 
            train_edges_dict = []
            for t in test_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                test_visit.append(t.get_post_order())
                test_edges_dict.append(t.get_edges_dict())
            for t in train_trees:
                t.link_nodes() # add children and parent info
                t.get_leaves_and_root() # find leaves and root for trees
                train_visit.append(t.get_post_order())
                train_edges_dict.append(t.get_edges_dict())
            for i in range(n_test):
                for j in range(n_train):
                    T1 = test_trees[i]
                    T2 = train_trees[j]
                    T1_visit = test_visit[i]
                    T2_visit = train_visit[j]
                    edict1 = test_edges_dict[i]
                    edict2 = train_edges_dict[j]
                    n1 = len(T1.nodes)
                    n2 = len(T2.nodes)
                    D = numpy.zeros((n1,n2))-1
                    D1 = numpy.zeros((n1,n1))-1
                    D2 = numpy.zeros((n2,n2))-1
                    # compare two nodes in two trees
                    for id1 in T1_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T2.leaves:
                                D[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict2[ee2]:
                                            s = self._gaussproduct(c1.mass,c2.mass,c1.inten,c2.inten,sm,si)
                                            score = score + (s + D[ii1,ii2])
                                D[ind1,ind2] = score
                    # compute normalizer
                    for id1 in T1_visit:
                        for id2 in T1_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T1.leaves or id2 in T1.leaves:
                                D1[ind1,ind2] = 0
                            else:
                                chl1 = T1.nodes[id1].children
                                chl2 = T1.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict1[ee1] == edict1[ee2]:
                                            s = self._gaussproduct(c1.mass,c2.mass,c1.inten,c2.inten,sm,si)
                                            score = score + (s + D1[ii1,ii2])
                                D1[ind1,ind2] = score

                    for id1 in T2_visit:
                        for id2 in T2_visit:
                            ind1 = int(id1[1:])-1
                            ind2 = int(id2[1:])-1
                            if id1 in T2.leaves or id2 in T2.leaves:
                                D2[ind1,ind2] = 0
                            else:
                                chl1 = T2.nodes[id1].children
                                chl2 = T2.nodes[id2].children
                                score = 0
                                for c1 in chl1:
                                    for c2 in chl2:
                                        ii1 = int(c1.id[1:])-1
                                        ii2 = int(c2.id[1:])-1
                                        ee1 = id1+'->'+c1.id
                                        ee2 = id2+'->'+c2.id
                                        if edict2[ee1] == edict2[ee2]:
                                            s = self._gaussproduct(c1.mass,c2.mass,c1.inten,c2.inten,sm,si)
                                            score = score + (s + D2[ii1,ii2])
                                D2[ind1,ind2] = score


                    sum1 = numpy.sum(D1)
                    sum2 = numpy.sum(D2)
                    
                    if sum1 == 0 or sum2 == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = numpy.sum(D)/numpy.sqrt(sum1*sum2)

        elif k_type == 'RLB':
            km = numpy.zeros((n_test,n_train))
            for t in test_trees:
                r_loss = t.get_root_loss()
            for t in train_trees:
                r_loss = t.get_root_loss()

            for i in range(n_test):
                T1 = test_trees[i]
                T1_rloss = T1.root_loss
                T1_dict = {}
                for nid, loss in T1_rloss:
                    T1_dict[loss] = 1

                for j in range(n_train):
                    T2 = train_trees[j]
                    T2_rloss = T2.root_loss
                    T2_dict = {}
                    for nid, loss in T2_rloss:
                        T2_dict[loss] = 1
                    
                    count = 0
                    for t1loss, dummy  in T1_dict.items():
                        if t1loss in T2_dict:
                            count +=1
                    if len(T1_rloss) == 0 or len(T2_rloss) == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = count/numpy.sqrt(len(T1_rloss)*len(T2_rloss))
                    #print "Computing kernel for %d, %d" % (i,j)

        elif k_type == 'RLI':
            km = numpy.zeros((n_test,n_train))
            for t in test_trees:
                r_loss = t.get_root_loss()
            for t in train_trees:
                r_loss = t.get_root_loss()
            for i in range(n_test):
                T1 = test_trees[i]
                T1_rloss = T1.root_loss
                T1_dict = {}
                for nid, loss in T1_rloss:
                    T1_dict[loss] = T1.nodes[nid].inten
                t1sum = 0
                for t1loss, t1inten in T1_dict.items():
                    t1sum = t1sum + t1inten*t1inten
                for j in range(n_train):
                    T2 = train_trees[j]
                    T2_rloss = T2.root_loss
                    T2_dict = {}
                    for nid, loss in T2_rloss:
                        T2_dict[loss] = T2.nodes[nid].inten
                    tsum = 0
                    for t1loss, t1inten in T1_dict.items():
                        if t1loss in T2_dict:
                            tsum = tsum + t1inten*T2_dict[t1loss]
                    t2sum = 0
                    for t2loss, t2inten in T2_dict.items():
                        t2sum = t2sum + t2inten*t2inten
                    if len(T1_rloss) == 0 or len(T2_rloss) == 0:
                        km[i,j] = 0
                    else:
                        km[i,j] = tsum/numpy.sqrt(t1sum*t2sum)
        else:
            raise Exception("Error in tree test kernel type")
        return km


    def _gaussproduct(self,m1,m2,i1,i2,sm,si):
        #return 1/numpy.sqrt((4*numpy.pi*si*sm))*numpy.exp(-0.25*numpy.square(i1-i2)/numpy.sqrt(si)-0.25*numpy.square(m1-m2)/numpy.sqrt(sm))
        return 0.25/(numpy.pi*numpy.sqrt(sm*si))*numpy.exp(-0.25*numpy.square(i1-i2)/si-0.25*numpy.square(m1-m2)/sm)

    def _normalize_km(self, km):
        n = len(km)
        for i in range(n):
            if km[i,i] == 0:
                km[i,i] = 1.0/100000
        return km / numpy.array(numpy.sqrt(numpy.mat(numpy.diag(km)).T * numpy.mat(numpy.diag(km))))

        
