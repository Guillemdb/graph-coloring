import os
import numpy as np
import networkx as nx
import pandas as pd
from time import perf_counter
from IPython.core.display import clear_output


class AntSolver(object):
    """This class implements a graph coloring algorithm described in the paper
    A Multiagent System for Frequency Assignment in Cellular Radio Networks
    (EEE TRANSACTIONS ON VEHICULAR TECHNOLOGY, VOL. 49, NO. 5, SEPTEMBER 2000),
    which is an stochastic algorithm that allows to solve the graph coloring problem
    
    >>> ant = AntSolver(G,n_ants=100,max_iters=15)
    >>> ant.solve()
    
    
    Parameters
    ----------
    g: networkx Graph; Graph to be colored. Color will be added as a node attribute.
    
    n_colors: int or None; Number of colors used to color the graph.
    
    pre_colored: bool, If set to true, the graph is assigned a random coloring during initialization.
    n_ants: int 15; Number of ants to use in each iteration.
    
    pc: float, 0<pc<1; Probability of moving to worst node. 
    
    pn: float, 0<pc<1; Probability of chosing best color. 
    
    max_iters: int; Number of max iterations for each ant.
    
    order_2: False; Use also the conflict level of direct neighbors when choosing
             the best color and worst node. Do not set to True, still experimental.
             
    verbose: True; Display information about the solver status when running
    
    Atributes
    ---------
    g: networkx Graph; Target graph where the optimization takes place.
    
    best: network Graph; Target graph with best colors found so far.
    
    cost_best: int; Total disagreement in the best graph.
    
    Methods
    -------
    solve: return None; Executes the ANT algorithm to color the graph
    
    random_node_coloring: see docs
    
    conflict_level: see docs
    
    """
    def __init__(self,
                 g,
                 n_colors=None,
                 n_ants=5,
                 pre_colored=True,
                 pc=0.8,
                 pn=0.8,
                 max_iters=1000000,
                 max_time=100000,
                 order_2=False,
                 verbose=True):
        
        self.order_2 = order_2
        self.verbose = verbose
        self.max_time = max_time
        #max number of colors
        max_c = max(nx.degree(g).values())+1
        self.n_colors =  max_c if n_colors is None else n_colors
        
        self.max_iters = max_iters
        self.n_ants = min(n_ants,g.number_of_nodes())
        self.pc = pc
        self.pn = pn
        self.g = g.copy()
        self.i = 0
        if not pre_colored:
            self.g = self.random_node_coloring(self.g,self.n_colors).copy()
        self.init_solver()
        

    @staticmethod
    def random_node_coloring(g,n=None):
        """Color a graph using a random scheme.
        n: int or None; Number of colors to use.
           If not n will be max degree(g)+1
        """
        max_c = max(nx.degree(g).values())+1
        n_colors =  max_c if n is None or n>max_c else n
        colors = np.arange(n_colors)
        for n in g.nodes_iter():
            g.node[n]['color'] = np.random.choice(colors)
        return g
    
    @staticmethod
    def conflict_level(g,n):
        """Conflict level in a node. This is the number
        of times two different colors are repeated on a node"""
        #get values of the colors conected to a given node
        colors =  [g.node[nei]for nei in nx.neighbors(g,n)]
        #how many of them are the same as the target node color
        return sum([g.node[n]['color'] == c['color'] for c in colors])
        
   
    def worst_adjacent_node(self,n,order_2=False):
        "node with greater number of conflicts"
        if order_2:
            adj_nodes = np.array([(nei,self.local_potential(self.g,nei)) for nei in nx.neighbors(self.g,n)])
        else:
            adj_nodes = np.array([(nei,self.g.node[nei]['cost']) for nei in nx.neighbors(self.g,n)])
        #neighbor with max cost
        return adj_nodes[adj_nodes[:,1].argmax(),0]
    
    
    def total_cost(self):
        """total number of conflicts in current graph"""
        return sum([self.g.node[n]['cost'] for n in self.g.nodes_iter()])
       
            
    def init_solver(self):
        """Initialization phase of the ANT algorithm. Called at instantiation"""       
        #init ants on random location
        self.ants = np.random.choice(np.arange(self.n_ants),self.n_ants,replace=False)
        #initialize local cost function
        for n in self.g.nodes_iter():
            self.g.node[n]['cost'] = self.conflict_level(self.g,n)
        self.best = self.g.copy()
        self.cost_best = self.total_cost()
        self.curr_cost = self.cost_best
        self._first_best = int(self.curr_cost)
    
  
    def local_potential(self,g,node):
        """Defined as the total numbre of discrepancies in the targed and adjacent nodes"""
        neigs = sum([self.conflict_level(g,nei) for nei in nx.neighbors(g,node)])
        return  self.conflict_level(g,node)+ neigs 
    
    def change_to_best(self,ant,order_2=True):
        """Change color to the best possible color for the current graph """
        g_test = self.g.subgraph([ant]+nx.neighbors(self.g,ant))
        cost = []
        #cost for changing to every color
        for c in range(self.n_colors):
            g_test.node[ant]['color'] = c
            cost.append(self.local_potential(g_test,ant) if order_2 else self.conflict_level(g_test,ant))
        #set the best color
        best = np.array(cost).argmin()
        self.g.node[ant]['color'] = best
    
    def change_to_random(self,ant):
        """Change one node to a random color"""
        self.g.node[ant]['color'] = np.random.choice(np.arange(self.n_colors))
    
    def print_status(self):
        """Display information about the solving process"""
        clear_output(True)
        imp = self.cost_best/self._first_best-1
        text = "--"*20+"""\n Solving G with {} nodes and {} edges: 
              \n {} colors needed and initial cost {} 
              \n- Epoch: {} Ant: {}
              \n- Best found: {} 
              \n- Current cost: {} 
              \n- Gain {}
              """.format(self.g.number_of_nodes(),
                                    self.g.number_of_edges(),
                                    self.n_colors,
                                    self._first_best,
                                    self.i,self._ant_i,
                                    self.cost_best,
                                    self.curr_cost,
                                    imp,)
        print(text)  
    
    def solve(self):
        """Execute ANT algorithm"""
        
        tic = perf_counter()
        self.i = 0
        print("--"*20+"""\n Solving g with {} nodes and {} edges: 
              \n {} colors needed and initial cost {}""".format(self.g.number_of_nodes(),
                                                                self.g.number_of_edges(),
                                                                self.n_colors,
                                                                self.cost_best))
        #iters and time limitations also
        while self.cost_best>0 and self.i<self.max_iters and perf_counter()-tic < self.max_time:
            self._ant_i = 0
            for ant in range(self.n_ants):
                if self.verbose:
                    self.print_status()
                n0 = self.ants[ant]
                if np.random.random() < self.pc:
                    #move to the worst adjacent node
                    n1 = self.worst_adjacent_node(n0,order_2=self.order_2)
                else:
                    #move randomly to an adjacent node
                    n1 = np.random.choice(nx.neighbors(self.g,n0))
                self.ants[ant] = n1 
                if np.random.random() < self.pn:
                    #Change color to best color
                    self.change_to_best(n1,order_2=self.order_2)
                else:
                    #change color to random color
                    self.change_to_random(n1)
                self.g.node[n1]['cost'] = self.conflict_level(self.g,n1)
                for nei in nx.neighbors(self.g,n1):
                    self.g.node[nei]['cost'] = self.conflict_level(self.g,nei)
                
                self._ant_i += 1
            self.curr_cost = self.total_cost()
            if self.curr_cost < self.cost_best:
                self.best = self.g.copy()
                self.cost_best = float(self.curr_cost)
            if self.cost_best == 0 or perf_counter()-tic > self.max_time:
                print("Finished in {0:.2f} seconds. Achieved cost {1:.0f}".format(perf_counter()-tic,self.cost_best))
                break
                self.i += 1
        if self.cost_best != 0 :
	   print("Ended after {0} iterations in {1:.2f} seconds and final conflict of {2:.0f}".format(self.i, perf_counter()-tic,self.cost_best))
