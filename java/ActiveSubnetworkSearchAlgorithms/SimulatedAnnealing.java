package ActiveSubnetworkSearchAlgorithms;

import ActiveSubnetworkSearchMisc.Subnetwork;
import ActiveSubnetworkSearchMisc.ScoreCalculations;
import Application.Parameters;
import Network.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

/**
 *
 * @author Ozan Ozisik Some code parts are
 * adapted from  https://github.com/idekerlab/jActiveModules
 *
 * Simulated Annealing for active subnetwork search.
 *
 * Notes: In code from idekerLab, all nodes are "on" in the beginning, which
 * cost lots of iterations to clean up. Here, if the related parameter is set, 
 * all the nodes with positive z-scores are on, others are off; else
 * they are set randomly. Keeping change was default behavior when
 * score was not improved in any subnetwork in the list and randomness did not
 * reject the change. I changed default to false, although this situation can
 * occur rarely.
 */
public class SimulatedAnnealing {

    /**
     * 
     * @param network
     * @param scoreCalculations 
     */
    public ArrayList<Subnetwork> simulatedAnnealing() {

        Network network=ActiveSubnetworkSearch.network;
        ScoreCalculations scoreCalculations=ActiveSubnetworkSearch.scoreCalculations;
                
        ArrayList<Node> nodeList = ActiveSubnetworkSearch.networkNodeList;

        HashSet<Node> nodesOnSet = new HashSet<>();
        HashSet<Node> nodesOffSet = new HashSet<>(nodeList);

        if(Parameters.startWithAllPositiveZScoreNodes){
            for (Node node : nodesOffSet) {
                if (scoreCalculations.getZScore(node) > 0) {
                    nodesOnSet.add(node);
                }
            }
        }else{
            for (Node node : nodesOffSet) {
                if(Math.random()<Parameters.geneInitialAdditionProbability){
                    nodesOnSet.add(node);
                }
            }
        }
        nodesOffSet.removeAll(nodesOnSet);

        SubnetworkFinder subnetworkFinder = new SubnetworkFinder();

        ArrayList<Subnetwork> subnetworkList = subnetworkFinder.findSubnetworksDFS(nodesOnSet);

        //subnetworkList.sort((subnetwork1, subnetwork2) -> - (int)Math.signum(subnetwork1.getScore()-subnetwork2.getScore()));
        Collections.sort(subnetworkList, Collections.reverseOrder());
        for (int i = 0; i < subnetworkList.size(); i++) {
            if (subnetworkList.get(i).numberOfNodes() > 1) {
                System.out.print(subnetworkList.get(i).numberOfNodes() + " " + (new DecimalFormat("###.##")).format(subnetworkList.get(i).getScore()) + ", ");
            }
        }
        System.out.println("");
        
        
        double initialTemperature = Parameters.sa_initialTemperature;
        double finalTemperature = Parameters.sa_finalTemperature;
        int totalIterations = Parameters.sa_totalIterations;

        double T = initialTemperature;
        double temp_step = 1 - Math.pow((finalTemperature / initialTemperature), (1.0 / totalIterations));
        Random rand = new Random();

        System.out.println("Percentage of finished job, node number and score of modules that have more than one node are as follows:");
        int percent=0;
        System.out.println("0%");
        //TODO: There should be another stop mechanism, not only iteration number
        for (int iteration = 0; iteration < totalIterations; iteration++) {
            int newPercent=(100*iteration)/totalIterations;
            if(newPercent>percent){
                percent=newPercent;
                System.out.println(percent+"% ");
                printSituation(subnetworkList);
            }
            
            

            Node node = nodeList.get(rand.nextInt(nodeList.size()));

            toggleNodeState(nodesOnSet, nodesOffSet, node);

            ArrayList<Subnetwork> newSubnetworkList = subnetworkFinder.findSubnetworksDFS(nodesOnSet);
            //newSubnetworkList.sort((subnetwork1, subnetwork2) -> -(int)Math.signum(subnetwork1.getScore()-subnetwork2.getScore()));
            Collections.sort(newSubnetworkList,Collections.reverseOrder());

            boolean decision = false;
            boolean keep = false;//was true in IdekerLab code

            Iterator<Subnetwork> oldIt = subnetworkList.iterator();
            Iterator<Subnetwork> newIt = newSubnetworkList.iterator();

            
            //Note: There is a higher chance of accepting a change
            while (!decision && (newIt.hasNext() && oldIt.hasNext())) {
                Subnetwork subnetworkOld=oldIt.next();
                Subnetwork subnetworkNew=newIt.next();
                double delta = subnetworkNew.getScore() - subnetworkOld.getScore();
                if (delta > .001) {
                    keep = true;
                    decision = true;
                }else if (rand.nextDouble() > Math.exp(delta / T)) {
                    keep = false;
                    decision = true;
                }
            }

            if (keep) {
                subnetworkList = newSubnetworkList;
            } else {
                toggleNodeState(nodesOnSet, nodesOffSet, node);
            }

            T = T * (1 - temp_step);
        }
        System.out.println("100%");
        printSituation(subnetworkList);
        
        return subnetworkList;
    }

    /**
     * Moves node from nodesOnSet to nodesOff set or vice versa.
     * @param nodesOnSet
     * @param nodesOffSet
     * @param node 
     */
    public void toggleNodeState(HashSet<Node> nodesOnSet, HashSet<Node> nodesOffSet, Node node) {
        if (nodesOnSet.contains(node)) {
            nodesOnSet.remove(node);
            nodesOffSet.add(node);
        } else {
            nodesOffSet.remove(node);
            nodesOnSet.add(node);
        }
    }
    
    public void printSituation(ArrayList<Subnetwork> subnetworkList){
        for (int i = 0; i < subnetworkList.size(); i++) {
            if (subnetworkList.get(i).numberOfNodes() > 1) {
                System.out.print(subnetworkList.get(i).numberOfNodes() + " " + (new DecimalFormat("###.##")).format(subnetworkList.get(i).getScore()) + ", ");
            }
        }
        System.out.println("");
    }
            
}
