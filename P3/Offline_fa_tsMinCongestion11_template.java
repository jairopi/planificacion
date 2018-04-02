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
public class Offline_fa_tsMinCongestion11_template implements IAlgorithm
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

			/* STUDENTS' CODE GOES HERE */
		
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

	
	
	
