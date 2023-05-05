package ActiveSubnetworkSearchMisc;

import Network.Network;
import ActiveSubnetworkSearchAlgorithms.ActiveSubnetworkSearch;
import Network.Node;
import java.util.ArrayList;
import java.util.HashSet;

/**
 *
 * @author Ozan Ozisik
 */
public class Subnetwork implements Comparable<Object> {
    
    private Network network;
    private ArrayList<Node> nodeList;
    private ScoreCalculations scoreCalculations;
    private double score;
    private double zSum;
    private HashSet<Node> neighborSet;
    private ArrayList<Node> neighborList;
    
    public Subnetwork(ArrayList<Node> nodeList){
        this.nodeList=nodeList;
        this.scoreCalculations=ActiveSubnetworkSearch.scoreCalculations;
        neighborSet=new HashSet<>();
        neighborList=new ArrayList();
        network = ActiveSubnetworkSearch.network;
        
        zSum=0;
        for(Node node:nodeList){
            zSum=zSum+scoreCalculations.getZScore(node);
        }
        this.score=scoreCalculations.calculateScoreOfSubnetwork(nodeList.size(), zSum, true);
    }
    
    //TODO It may be better to return a copy of the list or Collections.unmodifiableList(nodeList), here and in other private collection returning areas
    public ArrayList<Node> getNodeList(){
        return nodeList;
    }
    
    public HashSet<Node> getNeighborSet(){
        if(neighborSet.isEmpty()){
            extractNeighborSet();
        }
        return neighborSet;
    }
    public ArrayList<Node> getNeighborList(){
        if(neighborSet.isEmpty()){
            extractNeighborSet();
        }
        return neighborList;
    }
    
    
    public int numberOfNodes(){
        return nodeList.size();
    }
    
    public double getScore(){
        return score;
    }

    @Override
    public int compareTo(Object o) {
        return (int)Math.signum(this.getScore()-((Subnetwork)o).getScore());
    }
    
    public boolean contains(Node node){
        return nodeList.contains(node);
    }
    
    public void addNode(Node node){
        nodeList.add(node);
        zSum=zSum+scoreCalculations.getZScore(node);
        this.score=scoreCalculations.calculateScoreOfSubnetwork(nodeList.size(), zSum, true);
        
        //neighborSet is cleared for reextraction in case of need.
        //It could also be updated here.
        neighborSet.clear();
    }
    
    public void removeNode(Node node){
        if(nodeList.contains(node)){
            nodeList.remove(node);
            zSum=zSum-scoreCalculations.getZScore(node);
            this.score=scoreCalculations.calculateScoreOfSubnetwork(nodeList.size(), zSum, true);
            
            //neighborSet is cleared for reextraction in case of need.
            //It could also be updated here.
            neighborSet.clear();
        }
    }
    
    private void extractNeighborSet(){
        neighborSet.clear();
        for(Node node:nodeList){
            neighborSet.addAll(network.getNeighborSet(node));
        }
        neighborSet.removeAll(nodeList);
        neighborList.addAll(neighborSet);
    }
}
