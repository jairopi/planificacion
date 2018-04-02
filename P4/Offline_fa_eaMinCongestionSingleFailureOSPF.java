package pgr.year201516.p4;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import cern.colt.matrix.tdouble.DoubleFactory1D;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.libraries.IPUtils;
import com.net2plan.utils.Triple;

/** This algorithm computes the OSPF routing that minimizes the worst congestion between (i) the congestion when no failure occurs in the network, (ii) 
				the worse congestion ranging the cases when one single SRG is failed. An evolutionary algorithm is 
				implemented. Each solution is coded as an array of as many elements as links, storing the OSPF weight in each link. 
				Details on the particular form in which the evolutionary operators are implemented can be seen in the source code. 
				The returning design includes the routes defined by the OSPF weights, and an attribute per link with its OSPF weight  
 * @author Pablo Pavon-Marino */
public class Offline_fa_eaMinCongestionSingleFailureOSPF implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		int N = netPlan.getNumberOfNodes();
		int E = netPlan.getNumberOfLinks();
		int D = netPlan.getNumberOfDemands();
		if (N == 0 || D == 0 || E == 0) throw new Net2PlanException("This algorithm requires a topology with links, and a demand set");

		/* Initialize some variables */
		final int maxLinkWeight = Integer.parseInt(algorithmParameters.get("maxLinkWeight")); /* number of loopless candidate paths per demand */
		if (maxLinkWeight < 1) throw new Net2PlanException("The maximum link weight must be greater or equal than one" + maxLinkWeight);
		final double maxExecTimeSecs = Double.parseDouble(algorithmParameters.get("maxExecTimeSecs"));
		final int ea_populationSize = Integer.parseInt(algorithmParameters.get("ea_populationSize"));
		final int ea_offSpringSize = Integer.parseInt(algorithmParameters.get("ea_offSpringSize"));
		final double ea_fractionParentsChosenRandomly = Double.parseDouble(algorithmParameters.get("ea_fractionParentsChosenRandomly"));
		if (ea_offSpringSize > ea_populationSize) throw new Net2PlanException("The offspring size cannot exceed the population size divided by two");

		EvolutionaryAlgorithmCore gac = new EvolutionaryAlgorithmCore(netPlan, maxExecTimeSecs, maxLinkWeight, ea_populationSize, ea_fractionParentsChosenRandomly, ea_offSpringSize);
		gac.evolve();

		/* Save the best solution found into the netPlan object */
		IPUtils.setLinkWeights(netPlan, DoubleFactory1D.dense.make (gac.best_solution));
		IPUtils.setECMPForwardingRulesFromLinkWeights(netPlan, DoubleFactory1D.dense.make (gac.best_solution));

		return "Ok! Cost: " + gac.best_cost;
	}

	@Override
	public String getDescription()
	{
		return "This algorithm computes the OSPF routing that minimizes the worst congestion between (i) the congestion when no failure occurs in the network, (ii) " + "the worse congestion ranging the cases when one single SRG is failed. An evolutionary algorithm is " + "implemented. Each solution is coded as an array of as many elements as links, storing the OSPF weight in each link. "
			+ "Details on the particular form in which the evolutionary operators are implemented can be seen in the source code. " + "The returning design includes the routes defined by the OSPF weights, and an attribute per link with its OSPF weight";
	}


	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new ArrayList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("maxLinkWeight", "10", "Maximum OSPF weight associated to a link"));	
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "10", "Execution time of the algorithm. The algorithm will stop after this time, returning the best solution found"));	
		algorithmParameters.add(Triple.of("ea_populationSize", "1000", "Number of elements in the population"));
		algorithmParameters.add(Triple.of("ea_offSpringSize", "500", "Number of childs in the offspring every generation"));
		algorithmParameters.add(Triple.of("ea_fractionParentsChosenRandomly", "0.5", "Fraction of the parents that are selected randomly for creating the offspring"));
		return algorithmParameters;
	}

}
	
	
	
