import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.Triple;

public class TCA_nearestNeighborTSP implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		/* Initialize some input variables */
		



		
		
		
		/* Remove all the links in the current design */
		netPlan.removeAllLinks();






		/* loop */

		





















		/* Add the link to the solution */
		








		return "Ok! Ring total length: " + totalRingLength;
	}

	@Override
	public String getDescription()
	{
		return "Given a set of nodes, this heuristic tries to find a (possibly sub-optimal) minimum cost bidirectional ring (where the cost of a link "
				+ "is given by its length in km) using the nearest-neighbor greedy heuristic. The algorithm starts in a user-defined initial node, and in "
				+ "each iteration sets the next node to visit as the closest one to current node, that has not been visited yet. Link capacities "
				+ "are fixed to a user-defined constant value.";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new LinkedList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("linkCapacities", "100", "The capacities to set in the links"));
		algorithmParameters.add(Triple.of("initialNodeIndex", "0", "Index of the initial node in the algorithm"));
		return algorithmParameters;
	}
}
