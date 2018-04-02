/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mariño.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * 
 * Contributors:
 *     Pablo Pavon Mariño - initial API and implementation
 ******************************************************************************/



 




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
import com.net2plan.utils.Constants.OrderingType;
import com.net2plan.utils.DoubleUtils;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;

/**
 * This algorithm computes the bidirectional ring which minimizes the total 
 * ring length, using a GRASP heuristic. The cost of a link equals the euclidean 
 * distance between link end nodes. The algorithm executes GRASP iterations until 
 * the maxExecTime is reached. In each GRASP iteration, a solution is first created 
 * using a greedy-randomized approach. Then, this solution is the starting point 
 * of a local search (first-fit) heuristic, where two rings are considered neighbors 
 * when they have all but 2 bidirectional links in common. The greedy randomized 
 * approach is based on the concept of restricted candidate list (RCL). The ring 
 * starts with one single node chosen randomly. In each greedy iteration one node 
 * is added to the ring, randomly chosen from a RCL created in this iteration. 
 * The RCL contains the non-visited nodes which are closer to current node.
 * In particular, the nodes in the RCL are those which are at a distance from 
 * c_min to c_min + alpha * (c_max - c_min). c_min and c_max are the distances 
 * from last added node, to closest and furthest non-visited nodes. alpha is an 
 * algroithm parameter between 0 and 1, tuning the size of the RCL. When alpha=0, 
 * the algorithm is equal to pure greedy nearest neighbor heuristic. If alpha = 1, 
 * randomness is maximum: a non-visited node is chosen randomly.
 * 
 * @author Pablo Pavon-Marino, Jose-Luis Izquierdo-Zaragoza
 * @version 1.1, May 2015
 */
public class Offline_tca_graspTSP_template implements IAlgorithm
{
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		/* Initialize some variables */
		final double maxExecTimeSecs = Double.parseDouble(algorithmParameters.get("maxExecTimeSecs"));
		final double alpha = Double.parseDouble(algorithmParameters.get("alpha"));
		final double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));
		final long randomSeed = Long.parseLong(algorithmParameters.get("randomSeed"));
		final long algorithmEndTime = System.nanoTime() + (long) (1E9 * maxExecTimeSecs);
		Random r = new Random(randomSeed);

		/* Check input parameters */
		if ((alpha < 0) || (alpha > 1)) throw new Net2PlanException("The alpha factor must be between 0 and 1");
		if (linkCapacities < 0) throw new Net2PlanException("Link capacities must be a non-negative number");
		if (maxExecTimeSecs <= 0) throw new Net2PlanException("Algorithm running time must be a non-negative number");

		/* Compute the distances between each node pair */
		double[][] c_ij = netPlan.getMatrixNode2NodeEuclideanDistance().toArray();

		/* The best solution found so far (incumbent solution) is stored in these variables */
		ArrayList<Node> best_nodeSequence = new ArrayList<Node>();
		double best_cost = Double.MAX_VALUE;

		/* Main loop. Stop when maximum execution time is reached */
		while (System.nanoTime() < algorithmEndTime)
		{
			/******************* GREEDY-RANDOM STEP ***************************/
			Pair<ArrayList<Node>, Double> greedySolution = computeGreedySolution(netPlan , r, alpha, c_ij);

			/*** COMPUTE COST, AND CHECK IT IS A RING ***/
			checkRing(greedySolution, c_ij, "GREEDY SOLUTION");

			/******************* LOCAL SEARCH STEP ***************************/
			Pair<ArrayList<Node>, Double> lsSolution = computeLocalSearchStep(netPlan , greedySolution, c_ij);
			checkRing(lsSolution, c_ij, "LOCAL SEARCH SOLUTION");

			/* Check is the solution found improves best solution so far */
			if (lsSolution.getSecond() < best_cost)
			{
				checkRing(lsSolution, c_ij, "IMPROVEMENT!!");
				best_cost = lsSolution.getSecond();
				best_nodeSequence = (ArrayList<Node>) lsSolution.getFirst().clone();
			}
		}

		/* Save the best solution found into the netPlan object */
		netPlan.removeAllLinks();
		for (int cont = 0; cont < N - 1; cont++)
		{
			final Node n1 = best_nodeSequence.get(cont);
			final Node n2 = best_nodeSequence.get(cont + 1);
			netPlan.addLinkBidirectional(n1, n2, linkCapacities, c_ij[n1.getIndex ()][n2.getIndex ()], 200000 , null);
		}
		
		final Node firstRingN = best_nodeSequence.get(0);
		final Node lastRingN = best_nodeSequence.get(best_nodeSequence.size() - 1);
		netPlan.addLinkBidirectional(firstRingN, lastRingN, linkCapacities, c_ij[firstRingN.getIndex ()][lastRingN.getIndex ()], 200000 , null);

		/* Return printing the total distance of the solution */
		return "Ok! Cost : " + best_cost;
	}

	@Override
	public String getDescription()
	{
		return "This algorithm computes the bidirectional ring which minimizes the total ring length, using a GRASP heuristic. "
				+ "The cost of a link equals the euclidean distance between link end nodes. The algorithm executes GRASP iterations until the maxExecTime is reached."
				+ "In each GRASP iteration, a solution is first created using a greedy-randomized approach. Then, this solution is the starting point of a local search "
				+ "(first-fit) heuristic, where two rings are considered neighbors when they have all but 2 bidirectional links in common. "
				+ "The greedy randomized approach is based on the concept of restricted candidate list (RCL). The ring starts with one single node chosen randomly. In each greedy iteration"
				+ "one node is added to the ring, randomly chosen from a RCL created in this iteration. The RCL contains the non-visited nodes which are closer to current node."
				+ "In particular, the nodes in the RCL are those which are at a distance from c_min to c_min + alpha * (c_max - c_min). c_min and c_max are the distances "
				+ "from last added node, to closest and furthest non-visited nodes. alpha is an algroithm parameter between 0 and 1, tuning the size of the RCL. When alpha=0, the algorithm is equal to pure greedy"
				+ "nearest neighbor heuristic. If alpha = 1, randomness is maximum: a non-visited node is chosen randomly.";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new ArrayList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("alpha", "0.5", "Alpha factor for the greedy randomized step"));
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "10", "Execution time of the algorithm. The algorithm will stop after this time, returning the best solution found"));
		algorithmParameters.add(Triple.of("linkCapacities", "100", "Capacities to set to the links"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "Seed for the random number generator"));
		return algorithmParameters;
	}

	/**
	 * This function computes the randomized-greedy step in the GRASP. The ring 
	 * starts with a randomly chosen node. In each iteration, a restricted candidate 
	 * list is created with potential nodes to be added to the ring. These are 
	 * those non-visited nodes n, whose links (currentNode , n) are between c_min 
	 * and c_min + alpha (c_max - c_min), where c_min and c_max are the shortest 
	 * and longest distances to the non-visited nodes. Alpha is an algorithm 
	 * parameter 0...1. Then, next node is chosen randomly from the restricted 
	 * candidate list.
	 * 
	 * @since 1.0
	 */
	private Pair<ArrayList<Node>, Double> computeGreedySolution(NetPlan netPlan , Random r, double alpha, double[][] c_ij)
	{
		/* TO BE COMPLETED BY THE STUDENT */
	}

	/**
	 * This function receives a solution (node sequence and cost), and checks it 
	 * is a valid ring, and that the cost is well calculated.
	 * 
	 * @since 1.0
	 */
	private void checkRing(Pair<ArrayList<Node>, Double> sol, double[][] c_ij, String message)
	{
		final int N = c_ij.length;
		boolean[] inRing = new boolean[N];
		final ArrayList<Node> sol_ns = sol.getFirst();
		double checkCost = 0;
		if (sol_ns.size() != N) throw new RuntimeException("The solution created by the greed step is not a ring");

		for (int cont = 0; cont < N; cont++)
		{
			int a_e = sol_ns.get(cont).getIndex ();
			int b_e = (cont == N - 1) ? sol_ns.get(0).getIndex () : sol_ns.get(cont + 1).getIndex ();
			if (inRing[a_e] == true)
			{
				throw new Net2PlanException("The solution created by the greed step is not a ring");
			}
			else
			{
				checkCost += c_ij[a_e][b_e];
				inRing[a_e] = true;
			}
		}
		if (Math.abs(checkCost - sol.getSecond()) > 1e-3)
		{
			throw new RuntimeException("The solution cost is not well calculated. Solution true cost: " + checkCost + ", solution advertised cost: " + sol.getSecond());
		}
		System.out.println(message + " -- Cost: " + checkCost + ": Sequence of nodes: " + sol.getFirst());
	}


	/**
	 * This function starts from an initial solution. Applies a local search with 
	 * first-fit. Two neighbor solutions are those which have all the links in 
	 * common but two.
	 * 
	 * @since 1.0
	 */
	private Pair<ArrayList<Node>, Double> computeLocalSearchStep(NetPlan netPlan , Pair<ArrayList<Node>, Double> initialSolution, double[][] c_ij)
	{
		ArrayList<Node> nodeSequence = initialSolution.getFirst();
		double nodeSequenceCost = initialSolution.getSecond();
		final int N = c_ij.length;
		boolean costImprovementThisIteration = true;
		while (costImprovementThisIteration)
		{
			/* Generate the neighbors of current solution */
			costImprovementThisIteration = false;
			for (int n1 = 0; n1 < N - 2; n1++)
			{
				if (costImprovementThisIteration) break;

				final int e1_a = nodeSequence.get(n1).getIndex ();
				final int e1_b = nodeSequence.get(n1 + 1).getIndex ();
				for (int n2 = n1 + 2; n2 < N; n2++)
				{
					final int e2_a = nodeSequence.get(n2).getIndex ();
					final int e2_b = (n2 == N - 1) ? nodeSequence.get(0).getIndex () : nodeSequence.get(n2 + 1).getIndex ();
					if ((e1_b == e2_a) || (e2_b == e1_a)) continue;
					if ((e1_b == e2_b) || (e1_a == e2_a)) throw new RuntimeException("Bad");

					final double candidateCost = nodeSequenceCost + (c_ij[e1_a][e2_a] + c_ij[e1_b][e2_b] - c_ij[e1_a][e1_b] - c_ij[e2_a][e2_b]);
					if (candidateCost < nodeSequenceCost)
					{
						/* improvement */
						nodeSequenceCost = candidateCost;
						costImprovementThisIteration = true;
						ArrayList<Node> oldNodeSequence = (ArrayList<Node>) nodeSequence.clone();
						for (int cont = 0; cont < n2 - n1; cont++)
							nodeSequence.set(n1 + 1 + cont, oldNodeSequence.get(n2 - cont));
						nodeSequence.set((n2 + 1 == N) ? 0 : n2 + 1, netPlan.getNode (e2_b));
						break;
					}
				}
			}

		}
		return Pair.of(nodeSequence, nodeSequenceCost);
	}
}
