# TIFS-OptimumFusionWithByzantines
 Simulation files for the paper "A Game-Theoretic Framework for Optimum Decision Fusion in the Presence of Byzantines" published at IEEE transactions on Information Forensics and Security

 Authors: Andrea Abrardo, Mauro Barni, Kassem Kallas, Benedetta Tondi

 Abstract:
 Optimum decision fusion in the presence of malicious nodes - often referred to as Byzantines - is hindered by the necessity of exactly knowing the statistical behavior of Byzantines. In this paper, we focus on a simple, yet widely adopted, setup in which a fusion center (FC) is asked to make a binary decision about a sequence of system states by relying on the possibly corrupted decisions provided by local nodes. We propose a game-theoretic framework, which permits to exploit the superior performance provided by optimum decision fusion, while limiting the amount of a priori knowledge required. We use numerical simulations to derive the optimum behavior of the FC and the Byzantines in a game-theoretic sense, and to evaluate the achievable performance at the equilibrium point of the game. We analyze several different setups, showing that in all cases, the proposed solution permits to improve the accuracy of data fusion. We also show that, in some cases, it is preferable for the Byzantines to minimize the mutual information between the status of the observed system and the reports submitted to the FC, rather than always flipping the decision made by the local nodes.

 The repository contains three different folders for three different node states as explained in the paper which are: fixed node states, independent node states and maximum entropy node states. Inside each of the folders there is a main .m file to setup the parameters and run the simulations.

 The zerosum.m file is used to solve using linear programming the two-players zero sum game obtained by the simulations.
 The fnk.m and FNKmatrix.m files implement the dynamic programming algorithm explained in the paper in order to reduce the comutational cost of the algorithm.

 The paper can be found at: https://ieeexplore.ieee.org/abstract/document/7401067

 For any inquiries please contact Dr. Kassem Kallas at k_kallas@hotmail.com

