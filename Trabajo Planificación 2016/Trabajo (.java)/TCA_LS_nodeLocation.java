package Trabajo;
import java.util.ArrayList;
/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * 
 * Contributors:
 *     Pablo Pavon Mari�o - initial API and implementation
 ******************************************************************************/
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import cern.jet.math.tdouble.DoubleFunctions;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.utils.Triple;

/**
 * Given a set of access nodes, this algorithm computes the subset of access nodes which have a core node located next to it (in the same place),
 * and the links access-to-core nodes, so that the network cost is minimized. This cost is given by a cost per core node (always 1) plus a cost per link,
 * given by the product of link distance and the user-defined parameter linkCostPerKm. Access-core link capacities are fixed to the user-defined parameter
 * linkCapacities. This problem is solved using a local-search heuristic. Initially, every access node has a core node located next to it. In each iteration,
 * for every core node, I check the cost if we eliminated it and reconnected the access nodes. The best solution among these that also reduces the cost of
 * the current solution becomes the next solution. If no node can be eliminated without worsening the solution, a local optimal has been reached and the
 * algorithm stops.
 *
 * @author Pablo Pavon-Marino, Jose-Luis Izquierdo-Zaragoza
 * @version 1.1, May 2015
 */
public class TCA_LS_nodeLocation implements IAlgorithm
{
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();
		if (N == 0) throw new Net2PlanException("The input design must have nodes");

		/* Initialize some variables */
		final double linkCostPerKm = Double.parseDouble(algorithmParameters.get("linkCostPerKm"));
		final double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));
		final int initialNodeIndex = Integer.parseInt(algorithmParameters.get("initialNodeIndex"));
		final int maxNumAccessNodePerCoreNode = Integer.parseInt(algorithmParameters.get("maxNumAccessNodePerCoreNode"));
		final double coreNodeCost = Double.parseDouble(algorithmParameters.get("coreNodeCost"));
		final double maxExecTimeSecs = Double.parseDouble(algorithmParameters.get("maxExecTimeSecs"));
		final long randomSeed = Long.parseLong(algorithmParameters.get("randomSeed"));
		Random rnd = new Random(randomSeed);
		
		final long algorithmEndTime = System.nanoTime() + (long) (1E9 * maxExecTimeSecs);
		/* UN NODO CORE SOLO SE CONECTA A CINCO NODOS
		 * EL COSTE DE UN NODO TRONCAL YA NO ES 1 SI NO QUE ES CORENODECOST
		 *  
		 * 
		 */
		netPlan.removeAllLinks();
		netPlan.removeAllRoutes();
		
		if ((initialNodeIndex < 0) || (initialNodeIndex >= netPlan.getNumberOfNodes ())) throw new Net2PlanException ("Wrong initial node index");
		
		/* Compute the cost of one link */
		double[][] c_ij = netPlan.getMatrixNode2NodeEuclideanDistance().assign(DoubleFunctions.mult(linkCostPerKm)).toArray();
		
		/* The cost of a link ij is given by c_ij. 
		 * The cost of a core node in location j is coreNodeCost */
		int[] current_coreNodeConnectedToAccessNode = new int[N];
		
		//Calculamos el número de cores que tiene que tener nuestra solución 
		//Debido a la restricción del número máximo de enlaces a un nodo tenemos que cumplir la siguiente restricción
		int numCoreNodes = (N + maxNumAccessNodePerCoreNode - 1) / maxNumAccessNodePerCoreNode; 
	
		//Debugging
		System.out.println("Topology nodes: " + netPlan.getNumberOfNodes());
		System.out.println("Core nodes calculated: " + numCoreNodes);
		
		//Debugging
		//System.out.println("\r\n Los números aleatorios generados son: ");
		//Tenemos que generar numCoreNodes posiciones de current_coreNodeConnectedToAccessNode a 1
		List<Integer> initialCoreNodes = new ArrayList<Integer>();
		int VC = 0;
		while(initialCoreNodes.size()<numCoreNodes){
			int random = rnd.nextInt(N-1);
			if(!initialCoreNodes.contains(random)) {
				initialCoreNodes.add(random);	
				System.out.println(VC+ ": "+ random);
			}
			VC++;
		}
		
		System.out.println("initialCoreNodes.size() " + initialCoreNodes.size());
		//Me creo una lista de objetos NodoTroncal con los nodos troncales y me creo una lista de int con Nodos de Acceso
		List<NodoTroncal> nodosTroncales = new ArrayList<NodoTroncal>();
		List<Integer> nodosAcceso = new ArrayList<Integer>();
		
		//Tengo que rellenar con un 1 el array current_coreNodeConnectedToAccessNode en las posiciones de initialNodePositions
		//Para cada elemento del current_coreNodeConnectedToAccessNode tengo que comprobar si es troncal o enlace
				//Y añadirlo así a una lista u otra
		for(int i = 0; i< initialCoreNodes.size(); i++ ){
			current_coreNodeConnectedToAccessNode[initialCoreNodes.get(i)] = 1;
			NodoTroncal nodoTroncal = new NodoTroncal();
			nodoTroncal.setNodoTroncal(initialCoreNodes.get(i));
			nodosTroncales.add(nodoTroncal);
		}
		//Aquí relleno los nodos de acceso
		for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
			if(current_coreNodeConnectedToAccessNode[i]==0){
				nodosAcceso.add(i);
			}
		}
		
		//Debugging
		System.out.println("\r\n El array  current_coreNodeConnectedToAccessNode es: ");
		for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
			System.out.println(i+ ": "+ current_coreNodeConnectedToAccessNode[i]);	
		}
		
		//DEBUGGING
		System.out.println("El tamaño de nodosTroncales es" + nodosTroncales.size());
		System.out.println("El tamaño de nodosAcceso es" + nodosAcceso.size());
		System.out.println("\r\n El array de nodos trocales es: ");
		for(int i = 0; i<nodosTroncales.size(); i++){
			System.out.println(i+": " + nodosTroncales.get(i).getNodoTroncal());
		}
		System.out.println("\r\n El array de nodos de acceso es: ");
		for(int i = 0; i<nodosAcceso.size(); i++){
			System.out.println(i+ ": " + nodosAcceso.get(i));
		}
		int k = 0;
//		//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso más cercanos
		for(NodoTroncal currentNodoTroncal : nodosTroncales){
			//Debugging
			System.out.println(k +": Estamos trabajando con el nodo: " + currentNodoTroncal.getNodoTroncal());
			while(currentNodoTroncal.getNumberOfConnectedNodes()<maxNumAccessNodePerCoreNode-1){
				double bestDistance = Double.MAX_VALUE;
				int bestNeighborg = -1;
				int betNeighborgId = -1;
				for(int i = 0; i<nodosAcceso.size(); i++){
					double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getNodoTroncal()), netPlan.getNode(nodosAcceso.get(i)));
					if(distanceNodes<bestDistance){
						bestDistance = distanceNodes;
						bestNeighborg = nodosAcceso.get(i);
						betNeighborgId = i;
					}
				}
				if(bestNeighborg==-1) break;
				currentNodoTroncal.addConnectedNode(bestNeighborg);
				nodosAcceso.remove(betNeighborgId);
				System.out.println("El nodo más cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
				
			}
			k++;
		}
		
		//Debugging
		//Vamos a mostrar los nodosCore y sus enlaces
//		for(int i = 0; i<nodosTroncales.size(); i++){
//			System.out.println(i + ": " + nodosTroncales.get(i).getNodoTroncal());
//			System.out.println("Los nodos asociados al trocal son: ");
//			for(int j = 0; j < nodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
//				System.out.println(j + ": " +nodosTroncales.get(i).getConnectedNodes().get(j));
//			}
//		}
		double totalCost = 0;
		double totalKm = 0;
		double totalKmTroncales = 0;
		//Me voy a hacer un método que me genera todos los enlaces de los Core a los Accesos
		for(NodoTroncal nodoTroncal : nodosTroncales){
			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
				//netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
				totalCost = totalCost + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] * linkCostPerKm;
				totalKm = totalKm + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)];
				System.out.println("La distancia entre el nodo " + nodoTroncal.getNodoTroncal() + " y el nodo " + nodoTroncal.getConnectedNodes().get(i) + " es " + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] );
			}
		}

		//Coste de los nodos trocales con nodos de acceso
		totalCost = totalCost + nodosTroncales.size()*coreNodeCost;
	
		
		//Me falta la suma de los costes de los nodos trocales entre ellos
		for(int i = 0; i< nodosTroncales.size(); i++){
			for(int j = i+1; j < nodosTroncales.size(); j++){
				totalKmTroncales = totalKmTroncales + c_ij[nodosTroncales.get(i).getNodoTroncal()][nodosTroncales.get(j).getNodoTroncal()]; 
				totalCost = totalCost + c_ij[nodosTroncales.get(i).getNodoTroncal()][nodosTroncales.get(j).getNodoTroncal()] * linkCostPerKm;	
				System.out.println(i +" " + j + "Estoy sumando el coste de los enlaces entre los nodos: " + nodosTroncales.get(i).getNodoTroncal() + " y " + nodosTroncales.get(j).getNodoTroncal()   );	
			}
		}
		
		System.out.println("El coste de los nodo es: " +  nodosTroncales.size()*coreNodeCost);
		System.out.println("Los kilómetros totales de Core con enlaces son: " + totalKm);
		System.out.println("Los kilómetro de unir todos los troncales es: " + totalKmTroncales );
		System.out.println("El coste total de todo es : " + totalCost );
		
		

		double bestCost = totalCost;
		/*
		 * 
		 * 
		 * 
		 *Aquí tengo el bucle de mejora 
		 * */
		int h = 0;
		while(System.nanoTime()<algorithmEndTime){
			//System.out.println("Topology nodes: " + netPlan.getNumberOfNodes());
			//System.out.println("Core nodes calculated: " + numCoreNodes);
			
			Random newRnd = new Random(randomSeed+h+1);
			//System.out.println("Genero nuevos números aleatorios");
			List<Integer> newInitialCoreNodes = new ArrayList<Integer>();
			
			int VC2 = 0;
			while(newInitialCoreNodes.size()<numCoreNodes){
				int random = newRnd.nextInt(N-1);
				if(!newInitialCoreNodes.contains(random)) {
					newInitialCoreNodes.add(random);	
					//Debugging
				//	System.out.println(VC2 + ": "+ random);
				}
				VC2++;
			}
			
			//System.out.println("initialCoreNodes.size() " + newInitialCoreNodes.size());
			//Me creo una lista de objetos NodoTroncal con los nodos troncales y me creo una lista de int con Nodos de Acceso
			List<NodoTroncal> newNodosTroncales = new ArrayList<NodoTroncal>();
			List<Integer> newNodosAcceso = new ArrayList<Integer>();
			
			//Tengo que rellenar con un 1 el array current_coreNodeConnectedToAccessNode en las posiciones de initialNodePositions
			//Para cada elemento del current_coreNodeConnectedToAccessNode tengo que comprobar si es troncal o enlace
					//Y añadirlo así a una lista u otra
			for(int i = 0; i< newInitialCoreNodes.size(); i++ ){
				current_coreNodeConnectedToAccessNode[initialCoreNodes.get(i)] = 1;
				NodoTroncal nodoTroncal = new NodoTroncal();
				nodoTroncal.setNodoTroncal(newInitialCoreNodes.get(i));
				newNodosTroncales.add(nodoTroncal);
			}
			//Aquí relleno los nodos de acceso
			for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
				if(current_coreNodeConnectedToAccessNode[i]==0){
					newNodosAcceso.add(i);
				}
			}
			
			//Debugging
			//System.out.println("\r\n El array  current_coreNodeConnectedToAccessNode es: ");
			for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
				//System.out.println(i+ ": "+ current_coreNodeConnectedToAccessNode[i]);	
			}
			
			//DEBUGGING
			//System.out.println("El tamaño de nodosTroncales es" + newNodosTroncales.size());
			//System.out.println("El tamaño de nodosAcceso es" + newNodosAcceso.size());
			//System.out.println("\r\n El array de nodos trocales es: ");
			for(int i = 0; i<newNodosTroncales.size(); i++){
				//System.out.println(i+": " + newNodosTroncales.get(i).getNodoTroncal());
			}
			//System.out.println("\r\n El array de nodos de acceso es: ");
			for(int i = 0; i<newNodosAcceso.size(); i++){
				//System.out.println(i+ ": " + newNodosAcceso.get(i));
			}
			int newk = 0;
//			//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso más cercanos
			for(NodoTroncal currentNodoTroncal : newNodosTroncales){
				//Debugging
				//System.out.println(newk +": Estamos trabajando con el nodo: " + currentNodoTroncal.getNodoTroncal());
				while(currentNodoTroncal.getNumberOfConnectedNodes()<maxNumAccessNodePerCoreNode-1){
					double bestDistance = Double.MAX_VALUE;
					int bestNeighborg = -1;
					int betNeighborgId = -1;
					for(int i = 0; i<newNodosAcceso.size(); i++){
						double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getNodoTroncal()), netPlan.getNode(newNodosAcceso.get(i)));
						if(distanceNodes<bestDistance){
							bestDistance = distanceNodes;
							bestNeighborg = newNodosAcceso.get(i);
							betNeighborgId = i;
						}
					}
					if(bestNeighborg==-1) break;
					currentNodoTroncal.addConnectedNode(bestNeighborg);
					newNodosAcceso.remove(betNeighborgId);
					//System.out.println("El nodo más cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
					
				}
				newk++;
			}
			
			//Debugging
			//Vamos a mostrar los nodosCore y sus enlaces
//			for(int i = 0; i<nodosTroncales.size(); i++){
//				System.out.println(i + ": " + nodosTroncales.get(i).getNodoTroncal());
//				System.out.println("Los nodos asociados al trocal son: ");
//				for(int j = 0; j < nodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
//					System.out.println(j + ": " +nodosTroncales.get(i).getConnectedNodes().get(j));
//				}
//			}
			double newTotalCost = 0;
			double newTotalKm = 0;
			double newTotalKmTroncales = 0;
			
			//Me voy a hacer un método que me genera todos los enlaces de los Core a los Accesos
			for(NodoTroncal nodoTroncal : newNodosTroncales){
				for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
					//netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
					newTotalCost = newTotalCost + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] * linkCostPerKm;
					newTotalKm = newTotalKm + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)];
					//System.out.println("La distancia entre el nodo " + nodoTroncal.getNodoTroncal() + " y el nodo " + nodoTroncal.getConnectedNodes().get(i) + " es " + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] );
				}
			}

			//Coste de los nodos trocales con nodos de acceso
			newTotalCost = newTotalCost + newNodosTroncales.size()*coreNodeCost;
		
			
			//Me falta la suma de los costes de los nodos trocales entre ellos
			for(int i = 0; i< newNodosTroncales.size(); i++){
				for(int j = i+1; j < newNodosTroncales.size(); j++){
					newTotalKmTroncales = newTotalKmTroncales + c_ij[newNodosTroncales.get(i).getNodoTroncal()][newNodosTroncales.get(j).getNodoTroncal()]; 
					newTotalCost = newTotalCost + c_ij[newNodosTroncales.get(i).getNodoTroncal()][newNodosTroncales.get(j).getNodoTroncal()] * linkCostPerKm;	
					//System.out.println(i +" " + j + "Estoy sumando el coste de los enlaces entre los nodos: " + newNodosTroncales.get(i).getNodoTroncal() + " y " + newNodosTroncales.get(j).getNodoTroncal()   );	
				}
			}
			System.out.println("Para la iteración " + h + " los valores son:");
			System.out.println("El coste de los nodo es: " +  newNodosTroncales.size()*coreNodeCost);
			System.out.println("Los kilómetros totales de Core con enlaces son: " + newTotalKm);
			System.out.println("Los kilómetro de unir todos los troncales es: " + newTotalKmTroncales );
			System.out.println("El coste total de todo es : " + newTotalCost );
			
			if(newTotalCost<bestCost){
				bestCost = newTotalCost;
				System.out.println("La iteración " + h  + " produce un valor mejor " + bestCost);
			}
			h++;
		}
		return "Hellom my best cost is: " + bestCost ;
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
		algorithmParameters.add(Triple.of("linkCapacities", "1", "The capacities to set in the links"));
		algorithmParameters.add(Triple.of("linkCostPerKm", "1", "The cost per km of the links between access and core nodes"));
		algorithmParameters.add(Triple.of("initialNodeIndex", "0", "The index of the initial node in the algorithm"));
		algorithmParameters.add(Triple.of("maxNumAccessNodePerCoreNode", "5", "Max number of links to core node"));
		algorithmParameters.add(Triple.of("coreNodeCost", "100", "Cost of a core node"));
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "60", "Max execution time"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "The seed for the random number generator"));
		
		return algorithmParameters;
	}

	private static double newCostAfterAddingCoreNode(int newCoreNode, int currentNumberOfCoreNodes, int[] currentConnectivity, double[][] c_ij)
	{
		final double N = currentConnectivity.length;
		double newCost = currentNumberOfCoreNodes + 1; /* cost of adding one core node */
		for (int i = 0; i < N; i++)
		{
			final double costIfOldLink = c_ij[i][currentConnectivity[i]];
			final double costIfNewLink = c_ij[i][newCoreNode];
			newCost += Math.min(costIfNewLink, costIfOldLink);
		}
		
		return newCost;
	}
	
	
	public class NodoTroncal{
		
		int nodoTroncal;
		List<Integer> connectedNodes = new ArrayList<Integer>();
		
		public NodoTroncal(){}

		public int getNodoTroncal() {
			return nodoTroncal;
		}

		public void setNodoTroncal(int nodoTroncal) {
			this.nodoTroncal = nodoTroncal;
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


