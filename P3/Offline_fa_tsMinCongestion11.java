package pgr.year201516.p3;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.net2plan.interfaces.networkDesign.Demand;
import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.ProtectionSegment;
import com.net2plan.interfaces.networkDesign.Route;
import com.net2plan.utils.DoubleUtils;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;

/**  
 * @author Pablo Pavon-Marino
 * @version 1.0, May 2013
 * @since 1.0 */
public class Offline_fa_tsMinCongestion11 implements IAlgorithm
{
	private int N, E, D;
	private Map<Demand,List<Pair<List<Link>,List<Link>>>> cpl11;
	private double [] u_e;
	private double [] h_d;
	private NetPlan netPlan;
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Initialize some variables */
		this.N = netPlan.getNumberOfNodes();
		this.D = netPlan.getNumberOfDemands();
		this.E = netPlan.getNumberOfLinks();
		this.h_d = netPlan.getVectorDemandOfferedTraffic().toArray();
		this.u_e = netPlan.getVectorLinkCapacity().toArray();
		this.netPlan = netPlan;

		final int k = Integer.parseInt(algorithmParameters.get("k")); // number of loopless candidate paths per demand  
		final int ts_maxNumIt = Integer.parseInt(algorithmParameters.get("ts_maxNumIt"));
		final int ts_tabuListTenure = (int) Math.floor(D * Double.parseDouble(algorithmParameters.get("ts_tabuListTenureFraction")));

		/* Construct the candidate path list */
		Map<Demand,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(netPlan.getVectorLinkLengthInKm().toArray() , "K" , "" + k);
		this.cpl11 = NetPlan.computeUnicastCandidate11PathList(cpl , netPlan.getVectorLinkLengthInKm() , 2);

		/* Solution coding: an int[D] with the index of the 1+1 pair associated to this demand.   */
		int [] currentSolution_d = new int [D]; // initialize it so all demands use the first 1+1 pair in its 1+1 pair list  
		double currentCost = computeNetCongestion (currentSolution_d);
		
		/* Initialize tabu list: list of last demands changed */
		LinkedList<Integer> tabuList = new LinkedList<Integer> ();
		
		/* Main loop */
		double bestCost = currentCost;
		int [] bestSolution_d = Arrays.copyOf(currentSolution_d,D);
		for (int numIt = 0 ; numIt < ts_maxNumIt ; numIt ++)
		{
			/* For debug purposes: print information of the iteration */
			System.out.println("Starting iteration : " + numIt + ". Best cost: " + bestCost + ". Current cost: " + currentCost ); // . Current solution: " + Arrays.toString(current_solution));

			/* compute the neighborhood: change one demand */
			int [] bestAcceptableNeighbor_solution = null; // best neighbor solution found in this iteration
			double bestAcceptableNeighbor_cost = Double.MAX_VALUE; // cost of best neighbor solution found in this iteration
			int bestAcceptableNeighbor_changingDemand = -1; // the demand that changes in the neighbor respect to current solution: needed for log-term memory update
			
			/* Iterate for each neighbor: a neighbor solution differs from current in the 1+1 pair used in ONE demand  */
			for (int neighborChangingDemand = 0 ; neighborChangingDemand < D ; neighborChangingDemand ++)
			{
				for (int neighborNewPathPairId = 0 ; neighborNewPathPairId < cpl11.get(netPlan.getDemand(neighborChangingDemand)).size () ; neighborNewPathPairId ++) 
				{
					/* If the neighbor is the same as current solution, skip this iteration */
					if (neighborNewPathPairId == currentSolution_d [neighborChangingDemand]) continue; 

					/* create the neighbor solution, copying it in a new array */
					int [] neighborSolution_d = Arrays.copyOf(currentSolution_d, D);
					neighborSolution_d [neighborChangingDemand] = neighborNewPathPairId;
					
					/* compute the cost of the neighbor solution */
					double neighborCost = computeNetCongestion (neighborSolution_d);
					
					/* Ambition criteria: if the cost is better than best solution so far: accept it as next potential solution 
					 * (best solution, and best neighbor are updated)  */
					if (neighborCost < bestCost) 
					{ 
						bestSolution_d = Arrays.copyOf(neighborSolution_d,D);
						bestCost = neighborCost; 
						bestAcceptableNeighbor_solution = Arrays.copyOf(neighborSolution_d,D);
						bestAcceptableNeighbor_cost = neighborCost; 
						continue;
					}
					
					/* check if the neighbor solution is tabu: if so, skip it */
					if (tabuList.contains(neighborChangingDemand)) continue;
						
					/* the neighbor solution is not tabu: check if it is the best neighbor found so far */
					if (neighborCost < bestAcceptableNeighbor_cost)
					{
						bestAcceptableNeighbor_cost = neighborCost;
						bestAcceptableNeighbor_solution = Arrays.copyOf(neighborSolution_d, D);
						bestAcceptableNeighbor_changingDemand = neighborChangingDemand;
					}
				}
			}
			
			/* Move to the best neighbor found: current solution is the best neighbor */
			currentSolution_d = Arrays.copyOf(bestAcceptableNeighbor_solution, D);
			currentCost = bestAcceptableNeighbor_cost;
			
			/* Update the tabu list: add the new move to the list */
			tabuList.addLast(bestAcceptableNeighbor_changingDemand);
			/* Update the tabu list: if the tabu list is oversized, remove the first (oldest) element  */
			if (tabuList.size() > ts_tabuListTenure) tabuList.removeFirst();  
		}
		
		/* Save the best solution found into the netPlan object */
		netPlan.removeAllRoutes();
		netPlan.removeAllProtectionSegments();
		for (Demand d : netPlan.getDemands ())
		{
			final int bestPair = bestSolution_d [d.getIndex ()];
			final Pair<List<Link>,List<Link>> primaryAndBackupPath = cpl11.get(d).get(bestPair);
			final Route primary = netPlan.addRoute(d, d.getOfferedTraffic() , d.getOfferedTraffic() , primaryAndBackupPath.getFirst() , null);
			final ProtectionSegment backup = netPlan.addProtectionSegment(primaryAndBackupPath.getSecond(), d.getOfferedTraffic() , null);
			primary.addProtectionSegment(backup);
		}
		return "Ok! Cost : " + bestCost;
	}

	@Override
	public String getDescription()
	{
		return "";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new ArrayList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("k", "10", "Number of candidate paths per demand"));	
		algorithmParameters.add(Triple.of("ts_maxNumIt", "1000", "Number of iterations in the outer loop"));	
		algorithmParameters.add(Triple.of("ts_tabuListTenureFraction", "0.5", "Size of the tabu list as a fraction of the number of demands"));
		return algorithmParameters;
	}

	/* computes the congestion in the network (the objective function of the problem) */
	private double computeNetCongestion (int [] solution_d)
	{
		double [] y_e = new double [E];
		for (Demand d : netPlan.getDemands ())
		{
			final double offeredTrafficThisDemand = d.getOfferedTraffic();
			final int pathIndexThisDemand = solution_d [d.getIndex ()];
			final List<Link> seqLinksThisDemand_primary = cpl11.get (d).get (pathIndexThisDemand).getFirst();
			final List<Link> seqLinksThisDemand_backup = cpl11.get (d).get (pathIndexThisDemand).getSecond();
			for (Link e : seqLinksThisDemand_primary) 
				y_e [e.getIndex ()] += offeredTrafficThisDemand;
			for (Link e : seqLinksThisDemand_backup) 
				y_e [e.getIndex ()] += offeredTrafficThisDemand;
		}
		return DoubleUtils.maxValue(DoubleUtils.divide(y_e, u_e));
	}

}

	
	
	
