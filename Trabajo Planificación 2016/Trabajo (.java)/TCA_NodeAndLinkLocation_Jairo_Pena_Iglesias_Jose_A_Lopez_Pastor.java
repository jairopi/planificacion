package Trabajo;

import java.util.ArrayList;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import cern.jet.math.tdouble.DoubleFunctions;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.utils.Triple;


public class TCA_NodeAndLinkLocation_Jairo_Pena_Iglesias_Jose_A_Lopez_Pastor implements IAlgorithm
{
	private double linkCostPerKm;
	private double linkCapacities;
	private int maxNumAccessNodePerCoreNode;
	private double coreNodeCost;
	private double maxExecTimeSecs;
	private long randomSeed;
	
	private double[][] c_ij;
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		//Initial check
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		//Initial setup
		linkCostPerKm = Double.parseDouble(algorithmParameters.get("linkCostPerKm"));
		linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));
		maxNumAccessNodePerCoreNode = Integer.parseInt(algorithmParameters.get("maxNumAccessNodePerCoreNode"));
		coreNodeCost = Double.parseDouble(algorithmParameters.get("coreNodeCost"));
		maxExecTimeSecs = Double.parseDouble(algorithmParameters.get("maxExecTimeSecs"));
		randomSeed = Long.parseLong(algorithmParameters.get("randomSeed"));
		Random rnd = new Random(randomSeed);
		
		//Lists of nodes used in the algorithm
		List<CoreNode> coreNodes = new ArrayList<CoreNode>();
		List<Integer> accessNodes = new ArrayList<Integer>();
		List<CoreNode> bestSolutionNodeTrocal = new ArrayList<CoreNode>();
		List<Integer> bestSolutionAccesNode = new ArrayList<Integer>();
		
		double totalCost = 0;
		double bestCost = 0;
		
		//Used in the 60 secs loop 
		int h = 0;
		
		//Constant max time execution
		final long algorithmEndTime = System.nanoTime() + (long) (1E9 * maxExecTimeSecs);
		
		//Clearing links and routes
		netPlan.removeAllLinks();
		netPlan.removeAllRoutes();
		
		//Compute the cost of one link
		//The cost of a link ij is given by c_ij
		c_ij = netPlan.getMatrixNode2NodeEuclideanDistance().assign(DoubleFunctions.mult(linkCostPerKm)).toArray();
		int[] current_coreNodeConnectedToAccessNode = new int[N];
		
		//Min number of cores depending of maxNumAccessNodePerCoreNode
		int minCoreNodes = (N + maxNumAccessNodePerCoreNode - 1) / maxNumAccessNodePerCoreNode; 
		
		//First node choosen randomly
		int random = rnd.nextInt(N-1);
		
		//Look for the 40 nearest node
		//We use a List of objects with the distances sorted
		List<Distance> distances = new ArrayList<Distance>();
		for(int i=0; i<N; i++){
			Distance d = new Distance();
			d.setPosition(i);
			d.setDistance(c_ij[random][i]);
			distances.add(d);
		}
		
		Collections.sort(distances, new Comparator<Distance>() {
		    @Override
		    public int compare(Distance a1, Distance a2) {
		    	Double result = a1.getDistance() - a2.getDistance();
		        return result.intValue();
		    }
		});
				
		// The nearest core nodes are choosen
		for(int i = 0; i<minCoreNodes; i++){
			CoreNode nT = new CoreNode();
			nT.setCoreNode(distances.get(i).getPosition());
			coreNodes.add(nT);
			current_coreNodeConnectedToAccessNode[distances.get(i).getPosition()] = 1;
		}

		//The other ones are access nodes
		for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
			if(current_coreNodeConnectedToAccessNode[i]==0){
				accessNodes.add(i);
			}
		}
		
		//For each core node the nearest access cores are choosen. 
		//One access core is choosen in each iteration to minimize cost
		//This is the initial solution
		for(int j = 0; j< maxNumAccessNodePerCoreNode; j++){
			for(CoreNode currentNodoTroncal : coreNodes){
				double bestDistance = Double.MAX_VALUE;
				int bestNeighbor = -1;
				int bestNeighborId = -1;
				for(int i = 0; i<accessNodes.size(); i++){
					double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getCoreNode()), netPlan.getNode(accessNodes.get(i)));
					//The best neighbor is stored
					if(distanceNodes<bestDistance){
						bestDistance = distanceNodes;
						bestNeighbor = accessNodes.get(i);
						bestNeighborId = i;
					}
				}
				if(bestNeighbor==-1) break;
				currentNodoTroncal.addConnectedNode(bestNeighbor);
				accessNodes.remove(bestNeighborId);
			}
		}
		
		//Calculate the initial cost
		totalCost = calculateCost(coreNodes);
		
		//Debug. Comment in production
		System.out.println("Initial cost of core nodes are: " +  coreNodes.size()*coreNodeCost);
		System.out.println("Total cost is: " + totalCost );
		
		//The initial solution is stored as bestSolution
		bestSolutionNodeTrocal.addAll(coreNodes);
		bestSolutionAccesNode.addAll(accessNodes);
		bestCost = totalCost;
		
		//The iteration for 60 secs
		while(System.nanoTime()<algorithmEndTime){
			
			int[] newCurrent_coreNodeConnectedToAccessNode = new int[N];
			
			List<CoreNode> currentCoreNodes = new ArrayList<CoreNode>();
			List<Integer> currentAccessNodes = new ArrayList<Integer>();
			
			double currentTotalCost = 0;
			
			//New random initial core Node
			Random newRnd = new Random(randomSeed+h+1);
			int newRandom = newRnd.nextInt(N-1);
			h++;
			
			List<Distance> currenDistances = new ArrayList<Distance>();
			for(int i=0; i<N; i++){
				Distance d = new Distance();
				d.setPosition(i);
				d.setDistance(c_ij[newRandom][i]);
				currenDistances.add(d);
			}
			
			Collections.sort(currenDistances, new Comparator<Distance>() {
			    @Override
			    public int compare(Distance a1, Distance a2) {
			    	Double result = a1.getDistance() - a2.getDistance();
			        return result.intValue();
			    }
			});
			
			// The nearest core nodes are choosen
			for(int i = 0; i<minCoreNodes; i++){
				CoreNode nT = new CoreNode();
				nT.setCoreNode(currenDistances.get(i).getPosition());
				currentCoreNodes.add(nT);
				newCurrent_coreNodeConnectedToAccessNode[currenDistances.get(i).getPosition()] = 1;
			}
			//The other ones are access nodes
			for(int i = 0; i<newCurrent_coreNodeConnectedToAccessNode.length; i++){
				if(newCurrent_coreNodeConnectedToAccessNode[i]==0){
					currentAccessNodes.add(i);
				}
			}
			//For each core node the nearest access cores are choosen. 
			//One access core is choosen in each iteration to minimize cost		
			for(int j = 0; j< maxNumAccessNodePerCoreNode; j++){
				for(CoreNode currentNodoTroncal : currentCoreNodes){
					double bestDistance = Double.MAX_VALUE;
					int bestNeighbor = -1;
					int bestNeighborId = -1;
					for(int i = 0; i<currentAccessNodes.size(); i++){
						double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getCoreNode()), netPlan.getNode(currentAccessNodes.get(i)));
						if(distanceNodes<bestDistance){
							bestDistance = distanceNodes;
							bestNeighbor = currentAccessNodes.get(i);
							bestNeighborId = i;
						}
					}
					if(bestNeighbor==-1) break;
					currentNodoTroncal.addConnectedNode(bestNeighbor);
					currentAccessNodes.remove(bestNeighborId);
				}
			}
			
			//Calculate the current solution cost
			currentTotalCost = calculateCost(currentCoreNodes);
			System.out.println("en la iteracion " + h + " el coste es " + currentTotalCost);
			
			//If currentSolution is better than best solution, then best solution is updated
			if(currentTotalCost<bestCost){
				bestCost = currentTotalCost;
				bestSolutionNodeTrocal = new ArrayList<CoreNode>();
				bestSolutionNodeTrocal.addAll(currentCoreNodes);
			}
		}
		 
		//Links are added
		for(CoreNode nodoTroncal : bestSolutionNodeTrocal){
			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
				netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getCoreNode()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getCoreNode()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
			}
		}
		
		return "The total cost is: " + bestCost ;
	}

	//Calculate the cost of solution
	public double calculateCost(List<CoreNode> coreNodes){
		double totalCost = 0;
		
		//Cost of links between access and cores
		for(CoreNode nodoTroncal : coreNodes){
			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
				totalCost = totalCost + c_ij[nodoTroncal.getCoreNode()][nodoTroncal.getConnectedNodes().get(i)] * linkCostPerKm;
			}
		}
		
		//Cost of core nodes
		totalCost = totalCost + coreNodes.size()*coreNodeCost;
		
		//Cost of cores nodes links between them
		for(int i = 0; i< coreNodes.size(); i++){
			for(int j = i+1; j < coreNodes.size(); j++){
				totalCost = totalCost + c_ij[coreNodes.get(i).getCoreNode()][coreNodes.get(j).getCoreNode()] * linkCostPerKm;	
			}
		}
		return totalCost;
	}
	
	@Override
	public String getDescription()
	{
		return "The JAIRO - JALP ALGORITHM [JJA]. \r\n"
				+ "Based in greedy algorithm, we look for the nearest core nodes with the aim of minimize the cost of generate core-core links." 
				+ "Due to high number of cores, the distance between them is more important than the distance between access to core nodes."
				+ "Once the core nodes are choosen, a greedy algorithm is used to choose the nearest acccess nodes. ";
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new LinkedList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("linkCapacities", "1", "The capacities to set in the links"));
		algorithmParameters.add(Triple.of("linkCostPerKm", "1", "The cost per km of the links between access and core nodes"));
		algorithmParameters.add(Triple.of("maxNumAccessNodePerCoreNode", "5", "Max number of links to core node"));
		algorithmParameters.add(Triple.of("coreNodeCost", "100", "Cost of a core node"));
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "60", "Max execution time"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "The seed for the random number generator"));
		
		return algorithmParameters;
	}
	
	public class Distance{
		
		int position;
		double distance;
		
		public Distance(){}

		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public double getDistance() {
			return distance;
		}

		public void setDistance(double distance) {
			this.distance = distance;
		}
		
		
	}
	
	public class CoreNode{
		
		int coreNode;
		List<Integer> connectedNodes = new ArrayList<Integer>();
		
		public CoreNode(){}

		public int getCoreNode() {
			return coreNode;
		}

		public void setCoreNode(int coreNode) {
			this.coreNode = coreNode;
		}

		public List<Integer> getConnectedNodes() {
			return connectedNodes;
		}

		public void setConnectedNoder(List<Integer> connectedNodes) {
			this.connectedNodes = connectedNodes;
		}
		
		public void addConnectedNode(Integer node){
			this.connectedNodes.add(node);
		}
		
		public int getNumberOfConnectedNodes(){
			return this.connectedNodes.size();
		}
	} 
}


