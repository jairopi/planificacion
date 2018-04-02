import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cern.jet.math.tdouble.DoubleFunctions;
import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.utils.Triple;


public class TCA_LS_nodeLocation implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		/* Initialize some input variables */
		
		

		
		
		/* Compute the cost of each possible link */
		double[][] c_ij = netPlan.getMatrixNode2NodeEuclideanDistance().assign (DoubleFunctions.mult (linkCostPerKm)).toArray();
		

		/* Compute the cost of the initial solution */












		/*loop*/
		while (true)
		{
			/* Compute cost of neighbours */










			/* No improvement: end */
			if (bestNeighbor_newCost >= current_cost) break;

			/* Otherwise, take the best improving solution */
			
		




		}

		/* Remove all previous demands, links, protection segments, routes */
		netPlan.removeAllLinks();
		
		/* Save the new network links */








		return "Ok! Number of backbone nodes: " + (N - current_inactiveCoreNodes.size()) + ", cost: " + current_cost;
	}

	@Override
	public String getDescription()
	{
		return "Given a set of access nodes, this algorithm computes the subset of access nodes which have a core node located next to it (in the same place), "
				+ " and the links access-to-core nodes, so that the network cost is minimized. This cost is given by a cost per core node (always 1) plus a cost "
				+ "per link, given by the product of link distance and the user-defined parameter linkCostPerKm. Access-core link capacities are fixed "
				+ "to the user-defined parameter linkCapacities. This problem is solved using a local-search heuristic. Initially, every access node has a core "
				+ "node located next to it. In each iteration, for every core node, I check the cost if we eliminated it and reconnected the access nodes. "
				+ "The best solution among these that also reduces the cost of the current solution becomes the next solution. If no node can be eliminated "
				+ "without worsening the solution, a local optimal has been reached and the algorithm stops.";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new LinkedList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("linkCapacities", "100", "The capacities to set in the links"));
		algorithmParameters.add(Triple.of("linkCostPerKm", "0.001", "The cost per km of the links between access and core nodes"));
		algorithmParameters.add(Triple.of("initialNodeIndex", "0", "The index of the initial node in the algorithm"));
		return algorithmParameters;
	}

	private static double newCostAfterAddingCoreNode(int newCoreNode, int currentNumberOfCoreNodes, int[] currentConnectivity, double[][] c_ij)
	{
		final double N = currentConnectivity.length;
		




		
		return newCost;
	}
}
