package pgr.year201516.p5;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;

/** This algorithm computes the bidirectional ring which minimizes the total ring length, using an ACO (Ant Colony Optimization) heuristic, described as
 * Ant System in the literature: M. Dorigo, V. Maniezzo, A. Colorni, "Ant system: optimization by a colony of cooperating agents", IEEE T on Cybernetics, 1996.
 * The cost of a link equals the euclidean distance between link end nodes. The algorithm executes ACO iterations until the maxExecTime is reached.
 * In each ACO iteration, a loop for each ant is executed. Each ant, creates a greedy-randomized solution using the pheromones information of each potential link, 
 * and each link length, as follows. Each ant starts in a node chosen randomly. At each iteration, an ant in node n decides the next node to visit 
 * randomly among the non-visited nodes, being the probability of choosing node n' proportional to ph_nn'^alpha b_nn'^beta. ph_nn' is the amount of pheromones 
 * associated to link nn', b_nn' is the inverse of the distance between both nodes. alpha and beta are parameters tuning the importance of pheromones 
 * and link distances respectively. After all ants have finished, an evaporation strategy is executed, where each link nn' looses pheromones multiplicatively 
 * (pheromones are multiplied by 1-r, where r is a 0...1 evaporation factor). After evaporation phase, a reinforcement step is completed, where each ant a adds 
 * a quantity 1/La to the pheromones of all links traversed, being La the total distance of its ring. 
 * @author Pablo Pavon-Marino
 * @version 1.0, April 2014
 * @since 1.0 */
public class Offline_tca_acoTSP_template implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Initialize some variables */
		final int N = netPlan.getNumberOfNodes();
		final double maxExecTimeSecs = Double.parseDouble(algorithmParameters.get("maxExecTimeSecs"));
		final int numAnts = Integer.parseInt(algorithmParameters.get("numAnts"));
		final double alpha = Double.parseDouble(algorithmParameters.get("alpha"));
		final double beta = Double.parseDouble(algorithmParameters.get("beta"));
		final double evaporationFactor = Double.parseDouble(algorithmParameters.get("evaporationFactor"));
		final double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));
		final long randomSeed = Long.parseLong(algorithmParameters.get("randomSeed"));
		final long algorithmEndTime = System.nanoTime() + (long) (1E9 * maxExecTimeSecs);
		Random r = new Random (randomSeed);

		/* Check input parameters */
		if (numAnts <= 0) throw new Net2PlanException ("A positive number of ants is needed");
		if ((evaporationFactor < 0) || (evaporationFactor > 1))  throw new Net2PlanException ("The evaporationFactor must be between 0 and 1");
		if (linkCapacities < 0) throw new Net2PlanException ("Link capacities must be a non-negative number");
		if (maxExecTimeSecs <= 0) throw new Net2PlanException ("Algorithm running time must be a non-negative number");

		/* Compute the distances between each node pair */
		double [][] c_ij = netPlan.getMatrixNode2NodeEuclideanDistance().toArray();
		
		/* The best solution found so far (incumbent solution) is stored in these variables */
		ArrayList <Node> best_nodeSequence = new ArrayList<Node> ();
		double best_cost = Double.MAX_VALUE;

		/* Initialize some ACO control variables: the pheromones */
		double [][] pheromones_ij = new double [N][N];
		for (int n1 = 0 ; n1 < N ; n1 ++) for (int n2 = 0 ; n2 < N ; n2 ++) if (n1 != n2) pheromones_ij [n1][n2] = 1;

		/* Main loop. Stop when maximum execution time is reached */
		while (System.nanoTime () < algorithmEndTime)
		{
			ArrayList<Pair<ArrayList<Node>,Double>> antSolutions = new ArrayList<Pair<ArrayList<Node>,Double>> (numAnts); 
			for (int a = 0 ; a < numAnts ; a ++)
			{
				/* Build a greedy-random solution using pheromones info */
				Pair<ArrayList<Node>,Double> sol = computeAntSolution (netPlan , r , alpha , beta , pheromones_ij , c_ij);
				//checkRing (sol ,c_ij , "SOLUTION FOUND");

				/* Update incumbent solution */
				if (sol.getSecond() < best_cost)
				{
					best_cost = sol.getSecond();
					best_nodeSequence = (ArrayList<Node>) sol.getFirst().clone();
					checkRing (sol ,c_ij , "-- IMPROVING SOLUTION");
				}
				antSolutions.add(sol);
			}
			
			/* Apply evaporation strategy */
			for (int n1 = 0; n1 < N ; n1 ++)
				for (int n2 = 0; n2 < N ; n2 ++)
					if (n1 != n2)
						pheromones_ij [n1][n2] *= (1 - evaporationFactor);

			/* Apply reinforcement strategy */
			for (int a = 0 ; a < numAnts ; a ++)
			{
				ArrayList<Node> sol = antSolutions.get (a).getFirst();
				double benefit = 1 / antSolutions.get (a).getSecond();
				for (int cont = 0; cont < N-1 ; cont ++) 
					pheromones_ij [sol.get(cont).getIndex ()][sol.get(cont+1).getIndex ()] += benefit;
				pheromones_ij [sol.get(N-1).getIndex ()][sol.get(0).getIndex ()] += benefit;
			}
		}	

		/* Save the best solution found into the netPlan object */
		netPlan.removeAllLinks();
		for (int cont = 0 ; cont < N-1 ; cont ++)
		{
			final Node n1 = best_nodeSequence.get(cont); final Node n2 = best_nodeSequence.get(cont+1);
			netPlan.addLink(n1, n2, linkCapacities, c_ij[n1.getIndex ()][n2.getIndex ()] , 200000 , null); 
			netPlan.addLink(n2, n1, linkCapacities, c_ij[n2.getIndex ()][n1.getIndex ()] , 200000 , null); 
		}
		final Node firstRingNode =  best_nodeSequence.get(0); 
		final Node lastRingNode = best_nodeSequence.get(best_nodeSequence.size() - 1);
		netPlan.addLink(firstRingNode, lastRingNode, linkCapacities, c_ij[firstRingNode.getIndex ()][lastRingNode.getIndex ()] , 200000 , null); 
		netPlan.addLink(lastRingNode, firstRingNode, linkCapacities, c_ij[lastRingNode.getIndex ()][firstRingNode.getIndex ()] , 200000 , null); 
		
		/* Return printing the total distance of the solution */
		return "Ok! Cost : " + best_cost;
	}

	@Override
	public String getDescription()
	{
		return "This algorithm computes the bidirectional ring which minimizes the total ring length, using an ACO (Ant Colony Optimization) heuristic, described as" + 
		" Ant System in the literature: M. Dorigo, V. Maniezzo, A. Colorni, \"Ant system: optimization by a coloyn of cooperating agents\", IEEE T on Cybernetics, 1996." + 
		" The cost of a link equals the euclidean distance between link end nodes. The algorithm executes ACO iterations until the maxExecTime is reached." + 
		 " In each ACO iteration, a loop for each ant is executed. Each ant, creates a greedy-randomized solution using the pheromones information of each potential link," +  
		 " and each link length, as follows. Each ant starts in a node chosen randomly. At each iteration, an ant in node n decides the next node to visit " + 
		 " randomly among the non-visited nodes, being the probability of choosing node n' proportional to ph_nn'^alpha b_nn'^beta. ph_nn' is the amount of pheromones" +  
		 " associated to link nn', b_nn' is the inverse of the distance between both nodes. alpha and beta are parameters tuning the importance of pheromones " + 
		 " and link distances respectively. After all ants have finished, an evaporation strategy is executed, where each link nn' looses pheromones multiplicatively" +  
		 " (pheromones are multiplied by 1-r, where r is a 0...1 evaporation factor). After evaporation phase, a reinforcement step is completed, where each ant a adds " + 
		 " a quantity 1/La to the pheromones of all links traversed, being La the total distance of its ring.";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new ArrayList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "10", "Execution time of the algorithm. The algorithm will stop after this time, returning the best solution found"));	
		algorithmParameters.add(Triple.of("numAnts", "10", "Number of ants in the colony"));	
		algorithmParameters.add(Triple.of("alpha", "1", "Alpha factor tuning the pheromone influence in the ant movement"));	
		algorithmParameters.add(Triple.of("beta", "1", "Beta factor tuning the link distance influence in the ant movement"));	
		algorithmParameters.add(Triple.of("evaporationFactor", "0.5", "Factor controlling the evaporation of pheromones"));	
		algorithmParameters.add(Triple.of("linkCapacities", "100", "Capacities to set to the links"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "Seed for the random number genertor"));
		
		return algorithmParameters;
	}

	/* This function implements the greedy-randomized computation of a ring by an ant, as described in the class documentation */
	private Pair<ArrayList<Node>,Double> computeAntSolution (NetPlan netPlan , Random r , double alpha , double beta , double [][] pheromones_ij , double [][] c_ij)
	{
		/* TO BE IMPLEMENTED BY THE STUDENTS */
	}

	/* This function receives a solution (node sequence and cost), and checks it is a valid ring, and that the cost is well calculated */
	private void checkRing (Pair<ArrayList<Node>,Double> sol , double [][] c_ij , String message)
	{
		System.out.println(message + ": Sequence of nodes: " + sol.getFirst());
		final int N = c_ij.length;
		boolean [] inRing = new boolean [N];
		final ArrayList<Node> sol_ns = sol.getFirst();
		double checkCost = 0;
		if (sol_ns.size () != N) throw new RuntimeException ("The solution is not a ring");
		for (int cont = 0; cont < N ; cont ++) 
			{ int a_e = sol_ns.get(cont).getIndex (); int b_e = (cont == N-1)? sol_ns.get(0).getIndex () : sol_ns.get(cont+1).getIndex (); if (inRing [a_e] == true) throw new RuntimeException ("The solution is not a ring"); else { checkCost += c_ij [a_e][b_e]; inRing [a_e] = true; } }		if (Math.abs(checkCost-sol.getSecond ()) > 1e-3)  throw new RuntimeException ("The solution cost is not well calculated. Solution true cost: " + checkCost + ", solution advertised cost: " + sol.getSecond ());
		System.out.println(message + " -- Cost: " + checkCost + ": Sequence of nodes: " + sol.getFirst());
	}

	
	/* Receives a vector with values proportional to the probabilities p[s] of each option s. Returns the index of the sample chosen. 
	 * Each sample s has a probability to be chosen proportional to p[s] */
	private int sampleUniformDistribution (double [] p , double totalSumProbabilities , Random r)
	{
		final double comp = r.nextDouble() * totalSumProbabilities;
		double accumValue = 0;
		
		double checkTotalSum = 0; for (int cont = 0 ; cont < p.length ; cont ++) checkTotalSum += p [cont];
		if (Math.abs(checkTotalSum - totalSumProbabilities) > 1E-3) throw new RuntimeException ("Bad"); 
		for (int cont = 0 ; cont < p.length ; cont ++)
			{ accumValue += p [cont]; if (accumValue >= comp) return cont; }
		return p.length-1;
	}

}
	
	
	
