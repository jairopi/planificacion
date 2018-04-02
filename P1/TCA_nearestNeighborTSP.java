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

/**
 * Given a set of nodes, this heuristic tries to find a (possibly sub-optimal) minimum cost bidirectional ring (where the cost of a link is given by its
 * length in km) using the nearest-neighbor greedy heuristic. The algorithm starts in a user-defined initial node, and in each iteration sets
 * the next node to visit as the closest one to current node, that has not been visited yet. Link capacities are fixed to a user-defined constant value.
 *
 * @author Pablo Pavon-Marino, Jose-Luis Izquierdo-Zaragoza
 * @version 1.1, May 2015
 */
public class TCA_nearestNeighborTSP implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		/* Initialize some variables */
		final double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));
		final int initialNodeIndex = Integer.parseInt(algorithmParameters.get("initialNodeIndex"));
		
		if ((initialNodeIndex < 0) || (initialNodeIndex >= netPlan.getNumberOfNodes ())) throw new Net2PlanException ("Wrong initial node index");
		final Node initialNode = netPlan.getNode (initialNodeIndex);
		
		/* Remove all the links in the current design */
		netPlan.removeAllLinks();

		double totalRingLength = 0;
		Node currentNode = initialNode;
		Set<Node> notVisitedNode = new LinkedHashSet<Node>();
		for (Node n : netPlan.getNodes ())
			if (n != initialNode)
				notVisitedNode.add(n);

		while (notVisitedNode.size() > 0)
		{
			double bestDistance = Double.MAX_VALUE;
			Node bestNeighbor = null;
			for (Node node : notVisitedNode)
			{
				if (netPlan.getNodePairEuclideanDistance(currentNode, node) < bestDistance)
				{
					bestDistance = netPlan.getNodePairEuclideanDistance(currentNode, node);
					bestNeighbor = node;
				}
			}
			
			/* Add the link to the solution */
			netPlan.addLinkBidirectional(currentNode, bestNeighbor, linkCapacities, bestDistance, 200000 , null);
			totalRingLength += bestDistance;

			/* Update the current node */
			currentNode = bestNeighbor;
			notVisitedNode.remove(currentNode);
		}

		/* Add the link to the solution */
		final double lastLinkDistance = netPlan.getNodePairEuclideanDistance(initialNode, currentNode);
		netPlan.addLinkBidirectional(currentNode, initialNode, linkCapacities, lastLinkDistance, 200000 , null);
		totalRingLength += lastLinkDistance;

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
