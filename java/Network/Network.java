package Network;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Ozan Ozisik
 */
public class Network {
    
    private HashMap<Node, HashSet<Node>> adjacency;

    public Network() {
        adjacency=new HashMap<Node, HashSet<Node>>();
    }
    
    public void addInteraction(String strNode1, String strNode2){
        Node node1=new Node(strNode1);
        Node node2=new Node(strNode2);
        addInteraction(node1, node2);
    }
    
    public void addInteraction(Node node1, Node node2){
        if(node1.equals(node2)){
            Logger.getLogger(Network.class.getName()).log(Level.WARNING, "Self interaction discarded");
        }else{
            if(adjacency.get(node1)==null){
            adjacency.put(node1, new HashSet<Node>());
            }
            if(adjacency.get(node2)==null){
                adjacency.put(node2, new HashSet<Node>());
            }
            adjacency.get(node1).add(node2);
            adjacency.get(node2).add(node1);
        }
        
    }
    
    public HashSet<Node> getNeighborSet(Node node){
        return adjacency.get(node);
    }
    
    public ArrayList<Node> getNodeList(){
        return new ArrayList<>(adjacency.keySet());
    }
    
    public boolean areAdjacent(Node node1, Node node2){
        return adjacency.get(node1).contains(node2);
    }
    
    public int getNumberOfNodes(){
        return adjacency.keySet().size();
    }
    
    public int getNumberOfInteractions(){
        int interactionNumber=0;
        ArrayList<Node> nodeList=getNodeList();
        for(int i=0;i<nodeList.size()-1;i++){
            for(int j=i+1;j<nodeList.size();j++){
                if(areAdjacent(nodeList.get(i), nodeList.get(j))){
                    interactionNumber++;
                }
            }
        }
        return interactionNumber;
    }
}
