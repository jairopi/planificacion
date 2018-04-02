package pgr.year201516.p4;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.interfaces.networkDesign.SharedRiskGroup;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.libraries.IPUtils;
import com.net2plan.utils.Constants.OrderingType;
import com.net2plan.utils.Constants.SearchType;
import com.net2plan.utils.DoubleUtils;

/** This class implements the core operations of the genetic algorithm to solve the network design problem of the FA_EA_minCongestionSingleFailureOSPF algorithm.
 * @author Pablo Pavon-Marino
 * @version 1.0, April 2014
 * @since 1.0 */
public class EvolutionaryAlgorithmCore_template
{
	private NetPlan netPlan;
	private final int E;
	private final long algorithmEndTime;
	private final int maxLinkWeight;
	private final int ea_populationSize;
	private final double ea_fractionChosenRandomly;
	private final int ea_offspringSize;
	private final Random rnd;
	private List<double []> population;
	private double[] costs;
	public double best_cost;
	public double [] best_solution;

	/**
	 * Initialize the object with the problem data and the genetic algorithm specific parameters.
	 * 
	 * @since 1.0
	 */
	public EvolutionaryAlgorithmCore_template(NetPlan netPlan, double maxExecTimeSecs, int maxLinkWeight, int ea_populationSize, double ea_fractionParentsChosenRandomly, int ea_offSpringSize)
	{
		this.netPlan = netPlan;
		E = netPlan.getNumberOfLinks();

		this.maxLinkWeight = maxLinkWeight;
		this.ea_populationSize = ea_populationSize;
		this.ea_fractionChosenRandomly = ea_fractionParentsChosenRandomly;
		this.ea_offspringSize = ea_offSpringSize;
		rnd = new Random();
		algorithmEndTime = System.nanoTime() + (long) (maxExecTimeSecs * 1E9);
	}

	/**
	 * Computes the cost of a given solution.
	 * 
	 * @since 1.0
	 */
	private double computeCostSolution(double [] solution_e)
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * This method is called after creating the object, to launch the genetic 
	 * algorithm. The method will first create an initial population. Then, 
	 * will complete a number of iterations given by ea_numberIterations. In each 
	 * iteration, a population is transformed using the evolutionary operators: 
	 * parent selection, crossover, mutation, selection. The method keeps track 
	 * of the best solution found, which is kept in the internal variables 
	 * best_solution, best_cost. Then, they are publicly accessible when the 
	 * method returns.
	 * 
	 * @since 1.0
	 */
	public void evolve()
	{
		/* Generate the initial population */
		generateInitialSolutions();

		/* Update the best solution found so far */
		final int bestSolutionId = DoubleUtils.maxIndexes(costs, SearchType.FIRST)[0];
		best_cost = costs[bestSolutionId];
		best_solution = Arrays.copyOf (population.get(bestSolutionId) , E);
		System.out.println("Initial population. Best solution cost: " + best_cost);

		/* Evolve: one iteration per generation */
		int it = 0;
		while (System.nanoTime() < algorithmEndTime)
		{
			LinkedList<Integer> parents = operator_parentSelection();
			List<double []> offspring = operator_crossover(parents);
			this.operator_mutate(offspring);
			this.operator_select(offspring);
			if (costs[0] < best_cost)
			{
				best_cost = costs[0];
				best_solution = Arrays.copyOf (this.population.get(0) , E);
			}
			
			System.out.println("Iteration: " + it + ". Best solution cost: " + best_cost + ", best solution: " + best_solution);
			it++;
		}
	}

	/**
	 * Generates an initial population of solutions. Each solution is generated 
	 * choosing randomly for each link, a weight within the range 1...maxLinkWeight.
	 * 
	 * @since 1.0
	 */
	private void generateInitialSolutions()
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * Receive a set of parents of size ea_offspringSize * 2. From them, 
	 * generates an offspring of size ea_offspringSize. Couples are generated 
	 * randomly, each parent appears in only one couple. Given two parents, 
	 * the descendant is created by, for each link, picking from any parent 
	 * (chosen randomly), the link weight.
	 * 
	 * @since 1.0
	 */
	private List<double []> operator_crossover(LinkedList<Integer> parents)
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * Receives a population (typically an offspring). For each element in the 
	 * population, mutates the solution. A mutation consists of choosing a link 
	 * randomly, and changing randomly the link weight in the range 1...maxLinkWeight.
	 * 
	 * @since 1.0
	 */
	private void operator_mutate(List<double []> offspring)
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * Selects from the population ea_offspringSize*2 elements that will 
	 * become parents, generating the offspring. A fraction of the parents 
	 * (1-ea_fractionChosenRandomly) is chosen from the best in the population. 
	 * A fraction ea_fractionChosenRandomly is chosen randomly from all the 
	 * population, so that repetitions can occur.
	 *
	 * @since 1.0
	 */
	private LinkedList<Integer> operator_parentSelection()
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * Given the population from previous generation, and a new offspring, 
	 * selects the elements that will pass to the next generation. These will 
	 * be the best ea_populationSize among all the aggregated population.
	 * 
	 * @since 1.0
	 */
	private void operator_select(List<double []> offspring)
	{
		/* TO BE IMPLEMENTED */
	}

	/**
	 * Auxiliary function: Sorts a population in ascending order according to 
	 * its cost (best cost first).
	 * 
	 * @since 1.0
	 */
	private void sortAscending(List<double []> population, double[] costs)
	{
		double[] copyCosts = Arrays.copyOf(costs, costs.length);
		List<double []> copyPopulation = new ArrayList<double []>(population);

		int[] orderedIndexes = DoubleUtils.sortIndexes(costs, OrderingType.ASCENDING);
		for (int newPosition = 0; newPosition < orderedIndexes.length; newPosition++)
		{
			final int oldIndex = orderedIndexes[newPosition];
			costs[newPosition] = copyCosts[oldIndex];
			population.set(newPosition, copyPopulation.get(oldIndex));
		}
	}

	/**
	 * Auxiliary function: Computes the congestion in a network for a given 
	 * set of OSPF weights. Those links with a weight equal to Double.MAX_VALUE 
	 * are supposed to not exist => the routing behaves as if those links did 
	 * not exist.
	 *
	 * @since 1.0
	 */
	private double computeCongestion(double [] solution_e)
	{
		/* Compute the congestion when there are no failures */
		DoubleMatrix2D f_te = IPUtils.computeECMPRoutingTableMatrix_fte (netPlan.getNodes () , netPlan.getLinks () , DoubleFactory1D.dense.make (solution_e));
		DoubleMatrix1D y_e = GraphUtils.convert_fte2xde(netPlan.getNodes () , netPlan.getLinks () , netPlan.getDemands () , f_te).getThird();
		double congestion = 0;
		for(Link e : netPlan.getLinks ())
		{
			double rho_e = y_e.get(e.getIndex ()) == 0 ? 0 : y_e.get(e.getIndex ()) / e.getCapacity();
			if (rho_e > congestion) congestion = rho_e;
		}
		return congestion;
	}
	

}
