package ActiveSubnetworkSearchAlgorithms;

import ActiveSubnetworkSearchMisc.ScoreCalculations;
import ActiveSubnetworkSearchMisc.Subnetwork;
import Application.AppActiveSubnetworkSearch;
import Application.Parameters;
import File.ExperimentFileReader;
import File.SIFReader;
import Network.Network;
import Network.Node;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Ozan Ozisik
 */
public class ActiveSubnetworkSearch {
    
    /**
     * scoreCalculations and network are used in other classes
     */
    public static ScoreCalculations scoreCalculations;
    public static Network network;
    public static ArrayList<Node> networkNodeList;
    
    public static void activeSubnetworkSearch(){
        
        network=SIFReader.readSIF(Parameters.sifPath);
        if(network==null){
            Logger.getLogger(AppActiveSubnetworkSearch.class.getName()).log(Level.SEVERE, "SIF file could not be loaded");
            System.exit(0);
        }
        networkNodeList=network.getNodeList();
        
        ArrayList<SimpleEntry<String, Double>> namePValuePairList=ExperimentFileReader.readExperimentFile(Parameters.experimentFilePath);
        if(namePValuePairList==null){
            Logger.getLogger(AppActiveSubnetworkSearch.class.getName()).log(Level.SEVERE, "Experiment file could not be loaded");
            System.exit(0);
        }
        
        scoreCalculations=new ScoreCalculations(namePValuePairList);
        
        ArrayList<Subnetwork> subnetworkList;
        
        if(Parameters.useSAorGAorGR==Parameters.SearchMethod.GA){
            GeneticAlgorithm geneticAlgorithm=new GeneticAlgorithm();
            subnetworkList=geneticAlgorithm.geneticAlgorithm();
        }else if(Parameters.useSAorGAorGR==Parameters.SearchMethod.SA){
            SimulatedAnnealing simulatedAnnealing=new SimulatedAnnealing();
            subnetworkList=simulatedAnnealing.simulatedAnnealing();
        }else{
            GreedySearch greedySearch=new GreedySearch();
            subnetworkList=greedySearch.greedySearch();
        }
        
        try {
            BufferedWriter bw=new BufferedWriter(new FileWriter(Parameters.resultFilePath));
            for(Subnetwork subnetwork:subnetworkList){
                if(subnetwork.getScore()>0){
                    bw.write(subnetwork.getScore()+" ");
//                    bw.write(subnetwork.numberOfNodes()+" ");
                    for(Node node:subnetwork.getNodeList()){
                        bw.write(node.getName()+" ");
                    }
                    bw.newLine();
                }
            }
            bw.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ActiveSubnetworkSearch.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ActiveSubnetworkSearch.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
}
