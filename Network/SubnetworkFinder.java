package Network;

import ActiveSubnetworkSearchAlgorithms.ActiveSubnetworkSearch;
import ActiveSubnetworkSearchMisc.Subnetwork;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;


/**
 *
 * @author Ozan Ozisik
 */
public class SubnetworkFinder {
    Network network;
    HashSet<Node> nodesOnSet;
    HashSet<Node> reached;
    
    public SubnetworkFinder(){
        network=ActiveSubnetworkSearch.network;
    }
    
    
    /**
     * Finds the connected subnetworks of the given nodes using depth first
     * search. This method may return empty ArrayList, this should be handled in the 
     * calling methods.
     * @param nodesOnSet
     * @return 
     */ 
    public ArrayList<Subnetwork> findSubnetworksDFS(HashSet<Node> nodesOnSet){
        this.nodesOnSet=nodesOnSet;
        ArrayList<Subnetwork> subnetworkList=new ArrayList<Subnetwork>();
        reached=new HashSet<>(2 * nodesOnSet.size());
        
        for(Node node:nodesOnSet){
            if(!reached.contains(node)){
                ArrayList<Node> subnetworkNodeList = new ArrayList<>();
                search(node, subnetworkNodeList);
                subnetworkList.add(new Subnetwork(subnetworkNodeList));
            }
        }
        return subnetworkList;
    }
    private void search(Node node, ArrayList<Node> subnetworkNodeList){
        reached.add(node);
        subnetworkNodeList.add(node);
        HashSet<Node> neighborNodesSet=network.getNeighborSet(node);
        for(Node neighborNode:neighborNodesSet){
            if((nodesOnSet.contains(neighborNode))&&(!reached.contains(neighborNode))){
                search(neighborNode, subnetworkNodeList);
            }
        }
    }
    
    
    public ArrayList<Subnetwork> findSubnetworksDFSNonRecursive(HashSet<Node> nodesOnSet){
        ArrayList<Subnetwork> subnetworkList=new ArrayList<Subnetwork>();
        HashSet<Node> reached=new HashSet<>(2 * nodesOnSet.size());
        
        for(Node node:nodesOnSet){
            if(!reached.contains(node)){
                ArrayList<Node> subnetworkNodeList = new ArrayList<>();
                LinkedList<Node> nodesToBeChecked=new LinkedList<>();
                nodesToBeChecked.add(node);
                while(!nodesToBeChecked.isEmpty()){
                    Node curNode=nodesToBeChecked.pop();
                    if(!reached.contains(curNode)){
                        reached.add(curNode);
                        subnetworkNodeList.add(curNode);
                        for(Node neighborNode:network.getNeighborSet(curNode)){
                            if((nodesOnSet.contains(neighborNode))&&(!reached.contains(neighborNode))){
                                nodesToBeChecked.push(neighborNode);
                            }
                        }
                    }
                }
                subnetworkList.add(new Subnetwork(subnetworkNodeList));
            }
        }
        return subnetworkList;
    }
    
    
    
    
    public ArrayList<Subnetwork> findSubnetworksBFS(HashSet<Node> nodesOnSet){
        ArrayList<Subnetwork> subnetworkList=new ArrayList<Subnetwork>();
        HashSet<Node> reached=new HashSet<>(2 * nodesOnSet.size());
        
        for(Node node:nodesOnSet){
            if(!reached.contains(node)){
                ArrayList<Node> subnetworkNodeList = new ArrayList<>();
                LinkedList<Node> nodesToBeChecked=new LinkedList<>();
                nodesToBeChecked.add(node);
                while(!nodesToBeChecked.isEmpty()){
                    Node curNode=nodesToBeChecked.pop();
                    if(!reached.contains(curNode)){
                        reached.add(curNode);
                        subnetworkNodeList.add(curNode);
                        for(Node neighborNode:network.getNeighborSet(curNode)){
                            if((nodesOnSet.contains(neighborNode))&&(!reached.contains(neighborNode))){
                                nodesToBeChecked.add(neighborNode);
                            }
                        }
                    }
                }
                subnetworkList.add(new Subnetwork(subnetworkNodeList));
            }
        }
        return subnetworkList;
    }
}
