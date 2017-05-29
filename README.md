# ANT graph coloring algorithm implementation for networkx Graphs


This class implements a graph coloring algorithm described in the paper
[A Multiagent System for Frequency Assignment in Cellular Radio Networks](https://github.com/Guillem-db/graph-coloring/blob/master/IEEE.Freq.Assign.pdf)
(EEE TRANSACTIONS ON VEHICULAR TECHNOLOGY, VOL. 49, NO. 5, SEPTEMBER 2000),
which is an stochastic algorithm that allows to solve the graph coloring problem

	>>> ant = AntSolver(G,n_ants=100,max_iters=15)

	>>> ant.solve()
	
	----------------------------------------
   Solving G with 1000 nodes and 449735 edges: 
                
   930 colors needed and initial cost 978 
                
  - Epoch: 0 Ant: 49
                
  - Best found: 0 
                
  - Current cost: 2 
                
  - Gain -0.9979550102249489
                
  Finished in 1185.79 seconds. Achieved cost 0


## Parameters

**g**: networkx Graph; Graph to be colored. Color will be added as a node attribute.

**n_colors**: int or None; Number of colors used to color the graph.

**pre_colored**: bool, If set to true, the graph is assigned a random coloring during initialization.

**n_ants**: int 15; Number of ants to use in each iteration.

**pc**: float, 0<pc<1; Probability of moving to worst node. 

**pn**: float, 0<pc<1; Probability of chosing best color. 

**max_iters**: int; Number of max iterations for each ant.

**order_2**: False; Use also the conflict level of direct neighbors when choosing
     the best color and worst node. Do not set to True, still experimental.
     
**verbose**: True; Display information about the solver status when running

## Atributes

**g**: networkx Graph; Target graph where the optimization takes place.

**best**: network Graph; Target graph with best colors found so far.

**cost_best**: int; Total disagreement in the best graph.

## Methods

**solve**: return None; Executes the ANT algorithm to color the graph

**random_node_coloring**: see docs

**conflict_level**: see docs
    
