/*******************************************************************************
 * Copyright (c) 2016 Pablo Pavon Mariño.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 ******************************************************************************/
package pgr.year201516.p2;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.net2plan.interfaces.networkDesign.Demand;
import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Route;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.DoubleUtils;
import com.net2plan.utils.Triple;

/**
 * Given a set of nodes and links, with their capacities, an offered traffic vector, and a set of admissible paths for each demand (given by the k-shortest paths
 * in km for each demand, being k a user-defined parameter), this algorithm uses a simulated annealing metaheuristic to find the non-bifurcated routing solution
 * that minimizes the network congestion. Two solutions are considered neighbors if they differ the routing ofone demand. The number of iterations
 * in the inner (all with the same temperature) and outer loop (decreasing the temperature) are user-defined parameters.
 * The initial temperature is computed such that we accept a solution with a 5% of more congestion with 99% of probability. The last
 * temperature is such that with accept a solution with a 5% of more congestion with a 0.01% of probability. Temperature decreases geometrically, the alpha
 * factor of this progression is computed to match previous parameters.
 */
public class Offline_fa_sanMinCongestion_template implements IAlgorithm
{
	private int N , E , D;
	private double [] u_e;
	private Map<Demand,List<List<Link>>> cpl;
	private NetPlan netPlan;
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		this.N = netPlan.getNumberOfNodes();
		this.D = netPlan.getNumberOfDemands();
		this.E = netPlan.getNumberOfLinks();
		this.netPlan = netPlan;
		if (N == 0 || E == 0 || D == 0) throw new Net2PlanException("This algorithm requires a topology with links and a demand set");

		/* Initialize some variables */
		u_e = netPlan.getVectorLinkCapacity().toArray ();
	
		final int k = Integer.parseInt(algorithmParameters.get("k")); /* number of loopless candidate paths per demand */
		final int san_numOuterIterations = Integer.parseInt(algorithmParameters.get("san_numOuterIterations"));
		final int san_numInnerIterations = Integer.parseInt(algorithmParameters.get("san_numInnerIterations"));
		final long randomSeed = Long.parseLong (algorithmParameters.get("randomSeed"));
		Random rnd = new Random(randomSeed); /* random number generator */
		
		/* Initial temperature so that we accept bad solution with extra congestion 0.05 with probability 0.99 */
		final double san_initialTemperature = -1; // MODIFY BY THE STUDENT

		/* Final temperature so that we accept bad solution with extra congestion 0.05 with probability 0.001 */
		final double san_finalTemperature = -1; // MODIFY BY THE STUDENT

		/* Alpha factor so that starting in the initial temp, after exactly san_numOuterIterations we have a temperature equal to the final temp  */
		/* In each outer iteration, the new temperature t(i+1) = alpha*t(i) */
		final double san_alpha = -1; // MODIFY BY THE STUDENT

		System.out.println("Initial temp: " + san_initialTemperature + ", end temp: " + san_finalTemperature + ", alpha: " + san_alpha);

		/* Create one route per admissible path. The routes carry zero traffic */
		netPlan.removeAllUnicastRoutingInformation();
		netPlan.setRoutingType (RoutingType.SOURCE_ROUTING);
		/* map with the length of each link. Needed to compute the k admissible paths */
		this.cpl = netPlan.computeUnicastCandidatePathList(netPlan.getVectorLinkLengthInKm().toArray() , "K" , "" + k);

		/* SOLUTION CODING: A solution is of integers, with one element per demand. 
		 * Position d in the array contains the order of the path in the cpl object carrying the traffic of demand of index d */
		
		/* Construct an initial solution: each demand is carried through the shortest path in km (the first route in the list for that demand) */
		int [] currentSolution_d = new int [D]; 
		double currentCongestion = computeNetCongestion(currentSolution_d);

		/* Update the best solution found so far to the initial solution */
		int [] bestSolution_d = Arrays.copyOf (currentSolution_d , D);
		double bestCongestion = currentCongestion;

		/* Main algorithm loop --> to do by the students */


		/* Save the best solution found */
		netPlan.removeAllUnicastRoutingInformation();
		netPlan.setRoutingType(RoutingType.SOURCE_ROUTING);
		for (Demand d : netPlan.getDemands ())
		{
			final double offeredTraffic = d.getOfferedTraffic();
			final int selectedPathIndex = bestSolution_d [d.getIndex ()];
			final List<Link> selectedPath = cpl.get(d).get(selectedPathIndex);
			netPlan.addRoute (d , offeredTraffic , offeredTraffic , selectedPath , null);
		}

		return "Ok! Congestion : " + bestCongestion;
	}

	@Override
	public String getDescription()
	{
		return "Given a set of nodes and links, with their capacities, an offered traffic vector, and a set of admissible paths for each demand (given by "
				+ "the k-shortest paths in km for each demand, being k a user-defined parameter), this algorithm uses a simulated annealing metaheuristic"
				+ " to find the non-bifurcated routing solution that minimizes the network congestion."
				+ "Two solutions are considered neighbors if they differ the routing ofone demand. The number of iterations in the inner (all with the"
				+ " same temperature) and outer loop (decreasing the temperature) are user-defined parameters. The initial temperature is computed such "
				+ "that we accept a solution with a 5% of more congestion with 99% of probability. The last temperature is such that with accept a "
				+ "solution with a 5% of more congestion with a 0.01% of probability. Temperature decreases geometrically, the alpha factor of this"
				+ " progression is computed to match previous parameters.";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new LinkedList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("k", "100", "Number of candidate paths per demand"));
		algorithmParameters.add(Triple.of("san_numOuterIterations", "50", "Number of iterations in the outer loop (changing the temperature)"));
		algorithmParameters.add(Triple.of("san_numInnerIterations", "1000", "Number of iterations in the inner loop (all with the same temperature)"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "The seed for the random number generator"));
		return algorithmParameters;
	}

	/*
	 * Computes the congestion in the network (the objective function of the problem)
	 */
	private double computeNetCongestion(int [] solution_d)
	{
		double [] y_e = new double [E];
		for (Demand d : netPlan.getDemands ())
		{
			final double offeredTrafficThisDemand = d.getOfferedTraffic();
			final int pathIndexThisDemand = solution_d [d.getIndex ()];
			final List<Link> seqLinksThisDemand = cpl.get (d).get (pathIndexThisDemand);
			for (Link e : seqLinksThisDemand) 
				y_e [e.getIndex ()] += offeredTrafficThisDemand;
		}
		return DoubleUtils.maxValue(DoubleUtils.divide(y_e, u_e));
	}
}
