package ActiveSubnetworkSearchAlgorithms;

import ActiveSubnetworkSearchMisc.Subnetwork;
import Network.Node;
import Network.SubnetworkFinder;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author Ozan Ozisik
 */
public class GAIndividual implements Comparable<Object>{
    
    private ArrayList<Boolean> representationBoolean;
    private ArrayList<Node> networkNodeList;
    private HashSet<Node> nodesOnSet;
    private ArrayList<Subnetwork> subnetworkList;
    
    public GAIndividual(HashSet<Node> nodesOnSet){
        this.nodesOnSet=nodesOnSet;
        this.networkNodeList=ActiveSubnetworkSearch.networkNodeList;
        representationBoolean=new ArrayList<>();
        for(Node node:networkNodeList){
            if(nodesOnSet.contains(node)){
                representationBoolean.add(Boolean.TRUE);
            }else{
                representationBoolean.add(Boolean.FALSE);
            }
        }
        subnetworkList=(new SubnetworkFinder()).findSubnetworksDFS(nodesOnSet);
        Collections.sort(subnetworkList,Collections.reverseOrder());
    }
    
    public GAIndividual(ArrayList<Boolean> representationBoolean){
        this.representationBoolean=representationBoolean;
        this.networkNodeList=ActiveSubnetworkSearch.networkNodeList;
        nodesOnSet=new HashSet<>();
        for(int i=0;i<representationBoolean.size();i++){
            if(representationBoolean.get(i)){
                nodesOnSet.add(networkNodeList.get(i));
            }
        }
        subnetworkList=(new SubnetworkFinder()).findSubnetworksDFS(nodesOnSet);
        Collections.sort(subnetworkList,Collections.reverseOrder());
    }

//    @Override
//    public int compareTo(Object o) {
//        int result=0; 
//        result=(int)Math.signum(this.getScore()-((GAIndividual)o).getScore());
//        return result;
//    }

    @Override
    public int compareTo(Object o) {
        int result=0; 
        
        boolean decision = false;

        Iterator<Subnetwork> ownIt = this.subnetworkList.iterator();
        Iterator<Subnetwork> otherIt = ((GAIndividual)o).getSubnetworkList().iterator();

        while (!decision && (ownIt.hasNext() && otherIt.hasNext())) {
            Subnetwork subnetworkOwn=ownIt.next();
            Subnetwork subnetworkOther=otherIt.next();
            if (subnetworkOwn.getScore() > subnetworkOther.getScore()) {
                result=1;
                decision = true;
            }else if (subnetworkOwn.getScore() < subnetworkOther.getScore()) {
                result=-1;
                decision = true;
            }
        }
        
        return result;
    }

    public ArrayList<Boolean> getRepresentationBoolean() {
        return representationBoolean;
    }

    public ArrayList<Node> getNetworkNodeList() {
        return networkNodeList;
    }

    public HashSet<Node> getNodesOnSet() {
        return nodesOnSet;
    }

    public ArrayList<Subnetwork> getSubnetworkList() {
        return subnetworkList;
    }
    
    public Subnetwork getHighestScoringSubnetwork(){
        return subnetworkList.get(0);
    }
    
    
    /**
     * @return score of highest scoring subnetwork in the individual
     */
    public double getScore(){
        if(subnetworkList.isEmpty()){
            return 0;
        }else{
            return subnetworkList.get(0).getScore();
        }
    }
    
    public String toString(){
        String str="";
        for(Subnetwork subnetwork:subnetworkList){
            if(subnetwork.numberOfNodes()>1){
                str+=subnetwork.numberOfNodes()+" ";
                str+=(new DecimalFormat("###.##")).format(subnetwork.getScore())+", ";
            }
        }
        return str;
    }
}
