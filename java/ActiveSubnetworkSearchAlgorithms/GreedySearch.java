package ActiveSubnetworkSearchAlgorithms;

import ActiveSubnetworkSearchMisc.Subnetwork;
import Application.Parameters;
import Network.Network;
import Network.Node;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

/**
 *
 * @author Ozan Ozisik adapted from https://github.com/idekerlab/jActiveModules
 *
 */
public class GreedySearch {

    int max_depth, search_depth;
    HashMap<Node, Subnetwork> node2BestComponent;
    /**
     * Track the best score generated from the current starting point
     */
    double bestScore;

    /**
     * Map from a node to the number of nodes which are dependent on this node
     * for connectivity into the graph
     */
    HashMap<Node, Integer> node2DependentCount;

    /**
     * Map from a node to it's predecessor in the search tree When we remove
     * this node, that predecessor may be optionally added to the list of
     * removable nodes, dependending if it has any other predecessors
     */
    HashMap<Node, Node> node2Predecessor;

    /**
     * Lets us know if we need to repeat the greedy search from a new starting
     * point
     */
    boolean greedyDone;
    /**
     * Determines which nodes are within max depth of the starting point
     */
    HashSet<Node> withinMaxDepth;
    ArrayList<Node> nodeList;
    Network graph;

    public ArrayList<Subnetwork> greedySearch() {
        this.max_depth = Parameters.gr_maxDepth;
        this.search_depth = Parameters.gr_searchDepth;

        node2BestComponent = new HashMap<Node, Subnetwork>();
        nodeList = ActiveSubnetworkSearch.networkNodeList;
        graph = ActiveSubnetworkSearch.network;

        int percent=0;
        for(int nodeNo=0;nodeNo<nodeList.size();nodeNo++){
            Node seed=nodeList.get(nodeNo);
            int newPercent=(100*nodeNo)/nodeList.size();
            if(newPercent>percent){
                percent=newPercent;
                System.out.println(percent+"% of seeds checked");
            }
            
            
            withinMaxDepth = new HashSet();
            
            /**
             * determine which nodes are within max-depth
             * of this starting node and add them to a hash set
             * so we can easily identify them
             * if the user doesn't wish to limit the maximum
             * depth, just add every node into the max depth
             * hash, thus all nodes are accepted as possible
             * additions
             */
            if (max_depth==0) {
                for(Node node:nodeList){
                    withinMaxDepth.add(node);
                }
            } else {
                /**
                 * recursively find the nodes within a max depth
                 */ 
                initializeMaxDepth(seed, max_depth);
            }
            
            // set the neighborhood of nodes to initially be only
            // the single node we are starting the search from
            ArrayList<Node> nodeListForSubnetwork=new ArrayList();
            nodeListForSubnetwork.add(seed);
            Subnetwork component = new Subnetwork(nodeListForSubnetwork);
            
            // make sure that the seed is never added to the list of removables
            node2DependentCount = new HashMap();
            node2Predecessor = new HashMap();
            node2DependentCount.put(seed, 1);
            // we don't need to make a predecessor entry for the seed,
            // since it should never be added to the list of removable nodes
            HashSet<Node> removableNodes = new HashSet();
            bestScore = Double.NEGATIVE_INFINITY;
            runGreedySearchRecursive(search_depth, component, seed, removableNodes);            
            runGreedyRemovalSearch(component, removableNodes);
            
            for(Node node:component.getNodeList()){
                Subnetwork oldBest = node2BestComponent.get(node);
                if (oldBest == null || oldBest.getScore() < component.getScore()) {
                    node2BestComponent.put(node, component);
                }
            }
            
        }
        System.out.println("100%");
        
        ArrayList<Subnetwork> subnetworkList=new ArrayList<Subnetwork>(node2BestComponent.values());
        Collections.sort(subnetworkList,Collections.reverseOrder());
        
        System.out.println("Filtering");
        
        System.out.println("Subnetwork number"+subnetworkList.size());
        for(int i=subnetworkList.size()-1;i>=0;i--){
            if(subnetworkList.get(i).getScore()<=0){
                subnetworkList.remove(i);
            }
        }
        System.out.println("Subnetwork number"+subnetworkList.size());
        for(int i=subnetworkList.size()-1;i>=0;i--){
            if(subnetworkList.get(i).numberOfNodes()<2){
                subnetworkList.remove(i);
            }
        }
        System.out.println("Subnetwork number"+subnetworkList.size());
        
        return filterSubnetworkList(subnetworkList);
        
    }
    
    /**
     * Takes sorted subnetworkList and filters subnetworks using overlap threshold
     * @param subnetworkList 
     */
    private ArrayList<Subnetwork> filterSubnetworkList(ArrayList<Subnetwork> subnetworkList){
        ArrayList<Subnetwork> filteredSubnetworkList=new ArrayList<>();
        ArrayList<Subnetwork> subnetworkListToBeDeleted=new ArrayList<>();
        
        int percent=0;
        int i=0;
        while(i<subnetworkList.size()-1 && filteredSubnetworkList.size()<Parameters.gr_subnetworkNum){
            int newPercent=(100*i)/subnetworkList.size();
            if(newPercent>percent){
                percent=newPercent;
                System.out.println(percent+"% of subnetworks checked, "+(filteredSubnetworkList.size()+1) + " filtered subnetworks in the list");
            }
            Subnetwork subnetwork1=subnetworkList.get(i);
            if (!subnetworkListToBeDeleted.contains(subnetwork1)) {
                filteredSubnetworkList.add(subnetwork1);
                for (int j = i + 1; j < subnetworkList.size(); j++) {
                    Subnetwork subnetwork2 = subnetworkList.get(j);
                    if (!subnetworkListToBeDeleted.contains(subnetwork2)) {
                        int common = 0;
                        for (Node node1 : subnetwork1.getNodeList()) {
                            if(subnetwork2.contains(node1)){
                                common++;
                            }
                        }
                        int size;
                        if(subnetwork1.numberOfNodes()<subnetwork2.numberOfNodes()){
                            size=subnetwork1.numberOfNodes();
                        }else{
                            size=subnetwork2.numberOfNodes();
                        }
                        double overlap=common/((double)(size));
                        if(overlap>Parameters.gr_overlapThreshold){
                            //subnetwork2 is added because it has lower score
                            //subnetworkList is sorted
                            subnetworkListToBeDeleted.add(subnetwork2);
                        }
                    }
                }
            }
            i++;
        }
        System.out.println("100%");
        //subnetworkList.removeAll(subnetworkListToBeDeleted);
        return filteredSubnetworkList;
    }
    
    /**
     * Recursively find the nodes within a max depth
     */
    private void initializeMaxDepth(Node current, int depth) {
        withinMaxDepth.add(current);
        if (depth > 0) {
            for (Node neighbor : graph.getNeighborSet(current)) {
                if (!withinMaxDepth.contains(neighbor)) {
                    initializeMaxDepth(neighbor, depth - 1);
                }
            }
        }
    }


    /**
     * Recursive greedy search function. Called from greedySearch() to a
     * recursive set of calls to greedily identify high scoring networks. The
     * idea for this search is that we make a recursive call for each addition
     * of a node from the neighborhood. At each stage we check to see if we have
     * found a higher scoring network, and if so, store it in one of the global
     * variables. You know how in the Wonder Twins, one of them turned into an
     * elephant and the other turned into a bucket of water? This function is
     * like the elephant.
     *
     * @param depth The remaining depth allowed for this greed search.
     * @param component The current component we are branching from.
     * @param lastAdded The last node added.
     * @param removableNodes Nodes that can be removed.
     */
    private boolean runGreedySearchRecursive(int depth, Subnetwork component,
            Node lastAdded, HashSet<Node> removableNodes) {
        boolean improved = false;
        // score this component, check and see if the global top scores should
        // be updated, if we have found a better score, then return true
        if (component.getScore() > bestScore) {
            depth = search_depth;
            improved = true;
            bestScore = component.getScore();
        }

        if (depth > 0) {
            // if depth > 0, otherwise we are out of depth and the recursive
            // calls will end
            // Get an iterator of nodes which are next to the
            boolean anyCallImproved = false;
            removableNodes.remove(lastAdded);
            int dependentCount = 0;
            for(Node newNeighbor:graph.getNeighborSet(lastAdded)){
                //this node is only a new neighbor if it is not currently
                // in the component.
                if (withinMaxDepth.contains(newNeighbor)
                        && !component.contains(newNeighbor)) {
                    component.addNode(newNeighbor);
                    removableNodes.add(newNeighbor);
                    boolean thisCallImproved = runGreedySearchRecursive(
                            depth - 1, component, newNeighbor, removableNodes);
                    if (!thisCallImproved) {
                        component.removeNode(newNeighbor);
                        removableNodes.remove(newNeighbor);
                    }else {
                        dependentCount += 1;
                        anyCallImproved = true;
                        node2Predecessor.put(newNeighbor, lastAdded);
                    }
                }
            }
            improved = improved | anyCallImproved;
            if (dependentCount > 0) {
                removableNodes.remove(lastAdded);
                node2DependentCount.put(lastAdded, dependentCount);
            }

        }
        return improved;
    }

    private void runGreedyRemovalSearch(Subnetwork component, HashSet removableNodes) {
        LinkedList list = new LinkedList(removableNodes);
        while (!list.isEmpty()) {
            Node current = (Node) list.removeFirst();
            component.removeNode(current);
            double score = component.getScore();
            if (score > bestScore) {
                bestScore = score;
                Node predecessor = (Node) node2Predecessor.get(current);
                int dependentCount = node2DependentCount.get(predecessor);
                dependentCount -= 1;
                if (dependentCount == 0) {
                    removableNodes.add(predecessor);
                }else {
                    node2DependentCount.put(predecessor, dependentCount);
                }
            }else {
                component.addNode(current);
            }
        }
    }
}
