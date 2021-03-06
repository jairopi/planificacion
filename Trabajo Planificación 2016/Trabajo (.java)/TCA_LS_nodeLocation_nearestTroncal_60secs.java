package Trabajo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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
public class TCA_LS_nodeLocation_nearestTroncal_60secs implements IAlgorithm
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
		
		//Calculamos el n�mero de cores que tiene que tener nuestra soluci�n 
		//Debido a la restricci�n del n�mero m�ximo de enlaces a un nodo tenemos que cumplir la siguiente restricci�n
		int numCoreNodes = (N + maxNumAccessNodePerCoreNode - 1) / maxNumAccessNodePerCoreNode; 
	
		//Debugging
		System.out.println("Topology nodes: " + netPlan.getNumberOfNodes());
		System.out.println("Core nodes calculated: " + numCoreNodes);
		
		// Cogemos un nodo aleatorio
		int random = rnd.nextInt(N-1);
		System.out.println("Cominenzo en el nodo: " + random);
		
		//Tengo que buscar los 40 nodos m�s cercanos
		List<Distancia> arrayDeDistancias = new ArrayList<Distancia>();
		for(int i=0; i<N; i++){
			Distancia d = new Distancia();
			d.setPosicion(i);
			d.setDistancia(c_ij[random][i]);
			arrayDeDistancias.add(d);
		}
		
		Collections.sort(arrayDeDistancias, new Comparator<Distancia>() {
		    @Override
		    public int compare(Distancia a1, Distancia a2) {
		    	Double result = a1.getDistancia() - a2.getDistancia();
		        return result.intValue();
		    }
		});
		
		//Debugging
		for(int i = 0; i<arrayDeDistancias.size(); i++){
			System.out.println(arrayDeDistancias.get(i).getPosicion() + ": " + arrayDeDistancias.get(i).getDistancia());
		}
		
		
		List<NodoTroncal> nodosTroncales = new ArrayList<NodoTroncal>();
		List<Integer> nodosAcceso = new ArrayList<Integer>();
		// Cojo los  nodostrocales m�s cercanos
		for(int i = 0; i<numCoreNodes; i++){
			NodoTroncal nT = new NodoTroncal();
			nT.setNodoTroncal(arrayDeDistancias.get(i).getPosicion());
			nodosTroncales.add(nT);
			//Esto lo uso para hacerme la lista de nodos de acceso
			current_coreNodeConnectedToAccessNode[arrayDeDistancias.get(i).getPosicion()] = 1;
		}
		
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
		System.out.println("El tama�o de nodosTroncales es" + nodosTroncales.size());
		System.out.println("El tama�o de nodosAcceso es" + nodosAcceso.size());
		System.out.println("\r\n El array de nodos trocales es: ");
		for(int i = 0; i<nodosTroncales.size(); i++){
			System.out.println(i+": " + nodosTroncales.get(i).getNodoTroncal());
		}
		System.out.println("\r\n El array de nodos de acceso es: ");
		for(int i = 0; i<nodosAcceso.size(); i++){
			System.out.println(i+ ": " + nodosAcceso.get(i));
		}
		
		//Para cada nodo troncal cojo los nodos de acceso m�s cercanos
		int k = 0;
		//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso m�s cercanos
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
				System.out.println("El nodo m�s cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
				
			}
			k++;
		}
		//Tengo la soluci�n inicial
		
		//Debugging
		//Vamos a mostrar los nodosCore y sus enlaces
		for(int i = 0; i<nodosTroncales.size(); i++){
			System.out.println(i + ": " + nodosTroncales.get(i).getNodoTroncal());
			System.out.println("Los nodos asociados al trocal son: ");
			for(int j = 0; j < nodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
				System.out.println(j + ": " +nodosTroncales.get(i).getConnectedNodes().get(j));
			}
		}
		
		double totalCost = 0;
		double totalKm = 0;
		double totalKmTroncales = 0;
		
		//Me voy a hacer un m�todo que me genera todos los enlaces de los Core a los Accesos
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
		System.out.println("Los kil�metros totales de Core con enlaces son: " + totalKm);
		System.out.println("Los kil�metro de unir todos los troncales es: " + totalKmTroncales );
		System.out.println("El coste total de todo es : " + totalCost );
		
		
		List<NodoTroncal> bestSolutionNodeTrocal = new ArrayList<NodoTroncal>();
		List<Integer> bestSolutionAccesNode = new ArrayList<Integer>();
		
		bestSolutionNodeTrocal.addAll(nodosTroncales);
		bestSolutionAccesNode.addAll(nodosAcceso);
		double bestCost = totalCost;
		

		/*
		 *	AHORA COMENZAMOS LA ITERACI�N HASTA 6O SECS 
		 * 
		 */
		
		int h = 0;
		while(System.nanoTime()<algorithmEndTime){
			
			Random newRnd = new Random(randomSeed+h+1);
			int newRandom = rnd.nextInt(N-1);
			
			//System.out.println(" En la iteraci�n " + h + " cominenzo en el nodo: " + newRandom);
			
			//Tengo que buscar los 40 nodos m�s cercanos
			List<Distancia> newArrayDeDistancias = new ArrayList<Distancia>();
			for(int i=0; i<N; i++){
				Distancia d = new Distancia();
				d.setPosicion(i);
				d.setDistancia(c_ij[newRandom][i]);
				newArrayDeDistancias.add(d);
			}
			
			Collections.sort(newArrayDeDistancias, new Comparator<Distancia>() {
			    @Override
			    public int compare(Distancia a1, Distancia a2) {
			    	Double result = a1.getDistancia() - a2.getDistancia();
			        return result.intValue();
			    }
			});
			
			//Debugging
			for(int i = 0; i<arrayDeDistancias.size(); i++){
				//System.out.println(arrayDeDistancias.get(i).getPosicion() + ": " + arrayDeDistancias.get(i).getDistancia());
			}
			
			int[] newCurrent_coreNodeConnectedToAccessNode = new int[N];
			List<NodoTroncal> newNodosTroncales = new ArrayList<NodoTroncal>();
			List<Integer> newNodosAcceso = new ArrayList<Integer>();
			// Cojo los  nodostrocales m�s cercanos
			for(int i = 0; i<numCoreNodes; i++){
				NodoTroncal nT = new NodoTroncal();
				nT.setNodoTroncal(newArrayDeDistancias.get(i).getPosicion());
				newNodosTroncales.add(nT);
				//Esto lo uso para hacerme la lista de nodos de acceso
				newCurrent_coreNodeConnectedToAccessNode[newArrayDeDistancias.get(i).getPosicion()] = 1;
			}
			
			
			for(int i = 0; i<newCurrent_coreNodeConnectedToAccessNode.length; i++){
				if(newCurrent_coreNodeConnectedToAccessNode[i]==0){
					newNodosAcceso.add(i);
				}
			}

			//Debugging
			//System.out.println("\r\n El array  current_coreNodeConnectedToAccessNode es: ");
			for(int i = 0; i<newCurrent_coreNodeConnectedToAccessNode.length; i++){
				//System.out.println(i+ ": "+ newCurrent_coreNodeConnectedToAccessNode[i]);	
			}
			
			//DEBUGGING
			//System.out.println("El tama�o de nodosTroncales es" + newNodosTroncales.size());
			//System.out.println("El tama�o de nodosAcceso es" + newNodosAcceso.size());
			//System.out.println("\r\n El array de nodos trocales es: ");
			for(int i = 0; i<newNodosTroncales.size(); i++){
				//System.out.println(i+": " + newNodosTroncales.get(i).getNodoTroncal());
			}
			
			//System.out.println("\r\n El array de nodos de acceso es: ");
			for(int i = 0; i<newNodosAcceso.size(); i++){
				//System.out.println(i+ ": " + newNodosAcceso.get(i));
			}
			
			//Para cada nodo troncal cojo los nodos de acceso m�s cercanos
			int newk = 0;
			//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso m�s cercanos
			for(NodoTroncal currentNodoTroncal : newNodosTroncales){
				//Debugging
				//System.out.println(k +": Estamos trabajando con el nodo: " + currentNodoTroncal.getNodoTroncal());
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
					//System.out.println("El nodo m�s cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
				}
			newk++;
			}

			//Debugging
			//Vamos a mostrar los nodosCore y sus enlaces
			for(int i = 0; i<newNodosTroncales.size(); i++){
				//System.out.println(i + ": " + newNodosTroncales.get(i).getNodoTroncal());
				//System.out.println("Los nodos asociados al trocal son: ");
				for(int j = 0; j < newNodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
					//System.out.println(j + ": " + newNodosTroncales.get(i).getConnectedNodes().get(j));
				}
			}
			
			double newTotalCost = 0;
			double newTotalKm = 0;
			double newTotalKmTroncales = 0;
			

			//Me voy a hacer un m�todo que me genera todos los enlaces de los Core a los Accesos
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
			
			
			//System.out.println("Para la iteraci�n " + h + " El newcoste de los nodo es: " +  newNodosTroncales.size()*coreNodeCost);
			//System.out.println("Para la iteraci�n " + h + " Los newkil�metros totales de Core con enlaces son: " + newTotalKm);
			//System.out.println("Para la iteraci�n " + h + " Los newkil�metro de unir todos los troncales es: " + newTotalKmTroncales );
			//System.out.println("Para la iteraci�n " + h + " El newcoste total de todo es : " + newTotalCost );
		
		
			if(newTotalCost<bestCost){
				bestCost = newTotalCost;
				bestSolutionNodeTrocal = new ArrayList<NodoTroncal>();
				bestSolutionNodeTrocal.addAll(newNodosTroncales);
				
				System.out.println("Para la iteraci�n " + h + " los valores son:");
				System.out.println("El coste de los nodo es: " +  newNodosTroncales.size()*coreNodeCost);
				System.out.println("Los kil�metros totales de Core con enlaces son: " + newTotalKm);
				System.out.println("Los kil�metro de unir todos los troncales es: " + newTotalKmTroncales );
				System.out.println("El coste total de todo es : " + newTotalCost );
				System.out.println("La iteraci�n " + h  + " produce un valor mejor " + bestCost);
			}
			h++;
		}
		 
		//A�ado los enlaces adem�s de devolver el coste
		for(NodoTroncal nodoTroncal : bestSolutionNodeTrocal){
			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
				//System.out.println("Para el troncal: " + nodoTroncal.getNodoTroncal() + " esta conectando a " + nodoTroncal.getConnectedNodes().get(i));
				netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
			}
		}
		
		return "Hellom my best cost is: " + bestCost ;
			
			
			
			

			
		

			
//		
//		//Debugging
//		//System.out.println("\r\n Los n�meros aleatorios generados son: ");
//		//Tenemos que generar numCoreNodes posiciones de current_coreNodeConnectedToAccessNode a 1
//		List<Integer> initialCoreNodes = new ArrayList<Integer>();
//		int VC = 0;
//		while(initialCoreNodes.size()<numCoreNodes){
//			int random = rnd.nextInt(N-1);
//			if(!initialCoreNodes.contains(random)) {
//				initialCoreNodes.add(random);	
//				System.out.println(VC+ ": "+ random);
//			}
//			VC++;
//		}
//		
//		System.out.println("initialCoreNodes.size() " + initialCoreNodes.size());
//		//Me creo una lista de objetos NodoTroncal con los nodos troncales y me creo una lista de int con Nodos de Acceso
//		List<NodoTroncal> nodosTroncales = new ArrayList<NodoTroncal>();
//		List<Integer> nodosAcceso = new ArrayList<Integer>();
//		
//		//Tengo que rellenar con un 1 el array current_coreNodeConnectedToAccessNode en las posiciones de initialNodePositions
//		//Para cada elemento del current_coreNodeConnectedToAccessNode tengo que comprobar si es troncal o enlace
//				//Y a�adirlo as� a una lista u otra
//		for(int i = 0; i< initialCoreNodes.size(); i++ ){
//			current_coreNodeConnectedToAccessNode[initialCoreNodes.get(i)] = 1;
//			NodoTroncal nodoTroncal = new NodoTroncal();
//			nodoTroncal.setNodoTroncal(initialCoreNodes.get(i));
//			nodosTroncales.add(nodoTroncal);
//		}
//		//Aqu� relleno los nodos de acceso
//		for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
//			if(current_coreNodeConnectedToAccessNode[i]==0){
//				nodosAcceso.add(i);
//			}
//		}
//		
//		//Debugging
//		System.out.println("\r\n El array  current_coreNodeConnectedToAccessNode es: ");
//		for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
//			System.out.println(i+ ": "+ current_coreNodeConnectedToAccessNode[i]);	
//		}
//		
//		//DEBUGGING
//		System.out.println("El tama�o de nodosTroncales es" + nodosTroncales.size());
//		System.out.println("El tama�o de nodosAcceso es" + nodosAcceso.size());
//		System.out.println("\r\n El array de nodos trocales es: ");
//		for(int i = 0; i<nodosTroncales.size(); i++){
//			System.out.println(i+": " + nodosTroncales.get(i).getNodoTroncal());
//		}
//		System.out.println("\r\n El array de nodos de acceso es: ");
//		for(int i = 0; i<nodosAcceso.size(); i++){
//			System.out.println(i+ ": " + nodosAcceso.get(i));
//		}
//		int k = 0;
////		//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso m�s cercanos
//		for(NodoTroncal currentNodoTroncal : nodosTroncales){
//			//Debugging
//			System.out.println(k +": Estamos trabajando con el nodo: " + currentNodoTroncal.getNodoTroncal());
//			while(currentNodoTroncal.getNumberOfConnectedNodes()<maxNumAccessNodePerCoreNode-1){
//				double bestDistance = Double.MAX_VALUE;
//				int bestNeighborg = -1;
//				int betNeighborgId = -1;
//				for(int i = 0; i<nodosAcceso.size(); i++){
//					double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getNodoTroncal()), netPlan.getNode(nodosAcceso.get(i)));
//					if(distanceNodes<bestDistance){
//						bestDistance = distanceNodes;
//						bestNeighborg = nodosAcceso.get(i);
//						betNeighborgId = i;
//					}
//				}
//				if(bestNeighborg==-1) break;
//				currentNodoTroncal.addConnectedNode(bestNeighborg);
//				nodosAcceso.remove(betNeighborgId);
//				System.out.println("El nodo m�s cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
//				
//			}
//			k++;
//		}
//		
//		//Debugging
//		//Vamos a mostrar los nodosCore y sus enlaces
////		for(int i = 0; i<nodosTroncales.size(); i++){
////			System.out.println(i + ": " + nodosTroncales.get(i).getNodoTroncal());
////			System.out.println("Los nodos asociados al trocal son: ");
////			for(int j = 0; j < nodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
////				System.out.println(j + ": " +nodosTroncales.get(i).getConnectedNodes().get(j));
////			}
////		}
//		double totalCost = 0;
//		double totalKm = 0;
//		double totalKmTroncales = 0;
//		//Me voy a hacer un m�todo que me genera todos los enlaces de los Core a los Accesos
//		for(NodoTroncal nodoTroncal : nodosTroncales){
//			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
//				//netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
//				totalCost = totalCost + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] * linkCostPerKm;
//				totalKm = totalKm + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)];
//				System.out.println("La distancia entre el nodo " + nodoTroncal.getNodoTroncal() + " y el nodo " + nodoTroncal.getConnectedNodes().get(i) + " es " + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] );
//			}
//		}
//
//		//Coste de los nodos trocales con nodos de acceso
//		totalCost = totalCost + nodosTroncales.size()*coreNodeCost;
//	
//		
//		//Me falta la suma de los costes de los nodos trocales entre ellos
//		for(int i = 0; i< nodosTroncales.size(); i++){
//			for(int j = i+1; j < nodosTroncales.size(); j++){
//				totalKmTroncales = totalKmTroncales + c_ij[nodosTroncales.get(i).getNodoTroncal()][nodosTroncales.get(j).getNodoTroncal()]; 
//				totalCost = totalCost + c_ij[nodosTroncales.get(i).getNodoTroncal()][nodosTroncales.get(j).getNodoTroncal()] * linkCostPerKm;	
//				System.out.println(i +" " + j + "Estoy sumando el coste de los enlaces entre los nodos: " + nodosTroncales.get(i).getNodoTroncal() + " y " + nodosTroncales.get(j).getNodoTroncal()   );	
//			}
//		}
//		
//		System.out.println("El coste de los nodo es: " +  nodosTroncales.size()*coreNodeCost);
//		System.out.println("Los kil�metros totales de Core con enlaces son: " + totalKm);
//		System.out.println("Los kil�metro de unir todos los troncales es: " + totalKmTroncales );
//		System.out.println("El coste total de todo es : " + totalCost );
//		
//		
//		
//		List<NodoTroncal> bestSolutionNodeTrocal = new ArrayList<NodoTroncal>();
//		List<Integer> bestSolutionAccesNode = new ArrayList<Integer>();
//		
//		bestSolutionNodeTrocal.addAll(nodosTroncales);
//		bestSolutionAccesNode.addAll(nodosAcceso);
//		double bestCost = totalCost;
//		/*
//		 * 
//		 * 
//		 * 
//		 *Aqu� tengo el bucle de mejora 
//		 * */
//		int h = 0;
//		while(System.nanoTime()<algorithmEndTime){
//			//System.out.println("Topology nodes: " + netPlan.getNumberOfNodes());
//			//System.out.println("Core nodes calculated: " + numCoreNodes);
//			
//			Random newRnd = new Random(randomSeed+h+1);
//			//System.out.println("Genero nuevos n�meros aleatorios");
//			List<Integer> newInitialCoreNodes = new ArrayList<Integer>();
//			
//			int VC2 = 0;
//			while(newInitialCoreNodes.size()<numCoreNodes){
//				int random = newRnd.nextInt(N-1);
//				if(!newInitialCoreNodes.contains(random)) {
//					newInitialCoreNodes.add(random);	
//					//Debugging
//				//	System.out.println(VC2 + ": "+ random);
//				}
//				VC2++;
//			}
//			
//			//System.out.println("initialCoreNodes.size() " + newInitialCoreNodes.size());
//			//Me creo una lista de objetos NodoTroncal con los nodos troncales y me creo una lista de int con Nodos de Acceso
//			List<NodoTroncal> newNodosTroncales = new ArrayList<NodoTroncal>();
//			List<Integer> newNodosAcceso = new ArrayList<Integer>();
//			int[] newCurrent_coreNodeConnectedToAccessNode = new int[N];
//			//Tengo que rellenar con un 1 el array current_coreNodeConnectedToAccessNode en las posiciones de initialNodePositions
//			//Para cada elemento del current_coreNodeConnectedToAccessNode tengo que comprobar si es troncal o enlace
//					//Y a�adirlo as� a una lista u otra
//			for(int i = 0; i< newInitialCoreNodes.size(); i++ ){
//				newCurrent_coreNodeConnectedToAccessNode[newInitialCoreNodes.get(i)] = 1;
//				NodoTroncal nodoTroncal = new NodoTroncal();
//				nodoTroncal.setNodoTroncal(newInitialCoreNodes.get(i));
//				newNodosTroncales.add(nodoTroncal);
//			}
//			//Aqu� relleno los nodos de acceso
//			for(int i = 0; i<newCurrent_coreNodeConnectedToAccessNode.length; i++){
//				if(newCurrent_coreNodeConnectedToAccessNode[i]==0){
//					newNodosAcceso.add(i);
//				}
//			}
//			
//			//Debugging
//			//System.out.println("\r\n El array  current_coreNodeConnectedToAccessNode es: ");
//			for(int i = 0; i<current_coreNodeConnectedToAccessNode.length; i++){
//				//System.out.println(i+ ": "+ current_coreNodeConnectedToAccessNode[i]);	
//			}
//			
//			//DEBUGGING
//			//System.out.println("El tama�o de nodosTroncales es" + newNodosTroncales.size());
//			//System.out.println("El tama�o de nodosAcceso es" + newNodosAcceso.size());
//			//System.out.println("\r\n El array de nodos trocales es: ");
//			for(int i = 0; i<newNodosTroncales.size(); i++){
//				//System.out.println(i+": " + newNodosTroncales.get(i).getNodoTroncal());
//			}
//			//System.out.println("\r\n El array de nodos de acceso es: ");
//			for(int i = 0; i<newNodosAcceso.size(); i++){
//				//System.out.println(i+ ": " + newNodosAcceso.get(i));
//			}
//			int newk = 0; 
////			//Ahora para cada nodo troncal voy a buscar los 4 nodos de acceso m�s cercanos
//			for(NodoTroncal currentNodoTroncal : newNodosTroncales){
//				//Debugging
//				//System.out.println(newk +": Estamos trabajando con el nodo: " + currentNodoTroncal.getNodoTroncal());
//				while(currentNodoTroncal.getNumberOfConnectedNodes()<maxNumAccessNodePerCoreNode-1){
//					double bestDistance = Double.MAX_VALUE;
//					int bestNeighborg = -1;
//					int betNeighborgId = -1;
//					for(int i = 0; i<newNodosAcceso.size(); i++){
//						double distanceNodes = netPlan.getNodePairEuclideanDistance(netPlan.getNode(currentNodoTroncal.getNodoTroncal()), netPlan.getNode(newNodosAcceso.get(i)));
//						if(distanceNodes<bestDistance){
//							bestDistance = distanceNodes;
//							bestNeighborg = newNodosAcceso.get(i);
//							betNeighborgId = i;
//						}
//					}
//					if(bestNeighborg==-1) break;
//					currentNodoTroncal.addConnectedNode(bestNeighborg);
//					newNodosAcceso.remove(betNeighborgId);
//					//System.out.println("El nodo m�s cercano del nodo " + currentNodoTroncal.getNodoTroncal() + " es " + bestNeighborg);
//					
//				}
//				newk++;
//			}
//			
//			//Debugging
//			//Vamos a mostrar los nodosCore y sus enlaces
////			for(int i = 0; i<nodosTroncales.size(); i++){
////				System.out.println(i + ": " + nodosTroncales.get(i).getNodoTroncal());
////				System.out.println("Los nodos asociados al trocal son: ");
////				for(int j = 0; j < nodosTroncales.get(i).getNumberOfConnectedNodes(); j++){
////					System.out.println(j + ": " +nodosTroncales.get(i).getConnectedNodes().get(j));
////				}
////			}
//			double newTotalCost = 0;
//			double newTotalKm = 0;
//			double newTotalKmTroncales = 0;
//			
//			//Me voy a hacer un m�todo que me genera todos los enlaces de los Core a los Accesos
//			for(NodoTroncal nodoTroncal : newNodosTroncales){
//				for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
//					//netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
//					newTotalCost = newTotalCost + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] * linkCostPerKm;
//					newTotalKm = newTotalKm + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)];
//					//System.out.println("La distancia entre el nodo " + nodoTroncal.getNodoTroncal() + " y el nodo " + nodoTroncal.getConnectedNodes().get(i) + " es " + c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] );
//				}
//			}
//
//			//Coste de los nodos trocales con nodos de acceso
//			newTotalCost = newTotalCost + newNodosTroncales.size()*coreNodeCost;
//		
//			
//			//Me falta la suma de los costes de los nodos trocales entre ellos
//			for(int i = 0; i< newNodosTroncales.size(); i++){
//				for(int j = i+1; j < newNodosTroncales.size(); j++){
//					newTotalKmTroncales = newTotalKmTroncales + c_ij[newNodosTroncales.get(i).getNodoTroncal()][newNodosTroncales.get(j).getNodoTroncal()]; 
//					newTotalCost = newTotalCost + c_ij[newNodosTroncales.get(i).getNodoTroncal()][newNodosTroncales.get(j).getNodoTroncal()] * linkCostPerKm;	
//					//System.out.println(i +" " + j + "Estoy sumando el coste de los enlaces entre los nodos: " + newNodosTroncales.get(i).getNodoTroncal() + " y " + newNodosTroncales.get(j).getNodoTroncal()   );	
//				}
//			}
//			
//			
//			if(newTotalCost<bestCost){
//				bestSolutionNodeTrocal = new ArrayList<NodoTroncal>();
//				bestSolutionNodeTrocal.addAll(newNodosTroncales);
//				
//				bestCost = newTotalCost;
//				
//				System.out.println("Para la iteraci�n " + h + " los valores son:");
//				System.out.println("El coste de los nodo es: " +  newNodosTroncales.size()*coreNodeCost);
//				System.out.println("Los kil�metros totales de Core con enlaces son: " + newTotalKm);
//				System.out.println("Los kil�metro de unir todos los troncales es: " + newTotalKmTroncales );
//				System.out.println("El coste total de todo es : " + newTotalCost );
//				System.out.println("La iteraci�n " + h  + " produce un valor mejor " + bestCost);
//			}
//			h++;
//		}
//		 
//		//A�ado los enlaces adem�s de devolver el coste
//		for(NodoTroncal nodoTroncal : bestSolutionNodeTrocal){
//			for(int i = 0; i<nodoTroncal.getNumberOfConnectedNodes(); i++){
//				System.out.println("Para el troncal: " + nodoTroncal.getNodoTroncal() + " esta conectando a " + nodoTroncal.getConnectedNodes().get(i));
//				netPlan.addLinkBidirectional(netPlan.getNode(nodoTroncal.getNodoTroncal()) , netPlan.getNode(nodoTroncal.getConnectedNodes().get(i)), linkCapacities, c_ij[nodoTroncal.getNodoTroncal()][nodoTroncal.getConnectedNodes().get(i)] / linkCostPerKm, 200000 , null);
//			}
//		}
//		
//		return "Hellom my best cost is: " + bestCost ;
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
		algorithmParameters.add(Triple.of("maxExecTimeSecs", "20", "Max execution time"));
		algorithmParameters.add(Triple.of("randomSeed", "1", "The seed for the random number generator"));
		
		return algorithmParameters;
	}
	
	public class Distancia{
		int posicion;
		double distancia;
		
		public Distancia(){}
		
		public int getPosicion() {
			return posicion;
		}
		public void setPosicion(int posicion) {
			this.posicion = posicion;
		}
		public double getDistancia() {
			return distancia;
		}
		public void setDistancia(double distancia) {
			this.distancia = distancia;
		}
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


