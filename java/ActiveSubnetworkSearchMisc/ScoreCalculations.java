package ActiveSubnetworkSearchMisc;

import Application.Parameters;
import ActiveSubnetworkSearchAlgorithms.ActiveSubnetworkSearch;
import Network.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.AbstractMap.SimpleEntry;
import java.util.Collections;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Ozan Ozisik
 * Some code parts from https://github.com/idekerlab/jActiveModules
 * 
 * Everything related to score calculation is in this class.
 * 
 * As Monte Carlo approach includes randomness, the score calibrated by this
 * approach will be different in each run. 
 * 
 */
public class ScoreCalculations {

    private HashMap<Node, Double> nodeToPValueMap;
    private HashMap<Node, Double> nodeToZScoreMap;
    private final ArrayList<Node> networkNodeList;
    private double[] samplingScoreMeans;
    private double[] samplingScoreStds;
    private double[] samplingScoreMins;
    private double[] samplingScoreMaxs;
    
    private double MIN_SIG = 0.0000000000001;
    private double MAX_SIG = 1 - MIN_SIG;

    public ScoreCalculations(ArrayList<SimpleEntry<String, Double>> namePValuePairList) {
        this.networkNodeList=ActiveSubnetworkSearch.networkNodeList;
        fillNodeToPValueMap(namePValuePairList);
        process();
    }

    private void fillNodeToPValueMap(ArrayList<SimpleEntry<String, Double>> namePValuePairList) {
        nodeToPValueMap = new HashMap<Node, Double>();
        int geneFromExperimentNotExisingInNetwork = 0;
        for (SimpleEntry<String, Double> entry : namePValuePairList) {
            Node node = new Node(entry.getKey());
            if (networkNodeList.contains(node)) {
                double pValue = entry.getValue();
                if(pValue<MIN_SIG){
                    pValue=MIN_SIG;
                }else if(pValue>MAX_SIG){
                    pValue=MAX_SIG;
                }
                double existingPValue=nodeToPValueMap.get(node) == null ? 1 : nodeToPValueMap.get(node);
                if (pValue < existingPValue) {
                    nodeToPValueMap.put(node, pValue);
                }
            } else {
                geneFromExperimentNotExisingInNetwork++;
            }
        }
        System.out.println(nodeToPValueMap);
        if(geneFromExperimentNotExisingInNetwork>0){
            Logger.getLogger(ScoreCalculations.class.getName()).log(Level.WARNING, "{0} genes in experiment file does not exist in the network", geneFromExperimentNotExisingInNetwork);
        }
        
        //Assign p-value to genes that do not exist in the experiment file.
        for (Node node : networkNodeList) {
            if (!nodeToPValueMap.containsKey(node)) {
                nodeToPValueMap.put(node, Parameters.pForNonSignificantNodes);
            }
        }
    }

    public void process() {
        boolean tmpPenaltyForSize=Parameters.penaltyForSize;
        Parameters.penaltyForSize=false;
        calculateZScores();
        calculateMeanAndStdForMonteCarlo();
        Parameters.penaltyForSize=tmpPenaltyForSize;
    }

    public Double getPValue(Node node) {
        return nodeToPValueMap.get(node);
    }
    
    public Double getZScore(Node node) {
        return nodeToZScoreMap.get(node);
    }

    private void calculateZScores() {
        nodeToZScoreMap=new HashMap<Node, Double>();
        for (Node node : networkNodeList) {
            double pValue = nodeToPValueMap.get(node);
            nodeToZScoreMap.put(node, ZStatistics.oneMinusNormalCDFInverse(pValue));
        }
    }
    
    
    private void calculateMeanAndStdForMonteCarlo() {
        int numberOfNodes = networkNodeList.size();
        samplingScoreMeans = new double[numberOfNodes+1];//0th position is not used
        samplingScoreStds = new double[numberOfNodes+1];//0th position is not used
        samplingScoreMins = new double[numberOfNodes+1];//0th position is not used
        samplingScoreMaxs = new double[numberOfNodes+1];//0th position is not used
        
        double[] samplingScoreSums=new double[numberOfNodes+1];//0th position is not used
        double[] samplingScoreSquareSums=new double[numberOfNodes+1];//0th position is not used
                
        for (int i = 0; i < numberOfNodes+1; i++) {
            samplingScoreSums[i] = 0;
            samplingScoreSquareSums[i] = 0;
            samplingScoreMins[i]=Double.MAX_VALUE;
            samplingScoreMaxs[i]=Double.MIN_VALUE;
        }
        int numberOfTrials=2000;
        
        ArrayList<Node> nodeListForSampling=new ArrayList<>(networkNodeList);
        ArrayList<Node> significantNodesList=new ArrayList<>();
        ArrayList<Node> nonsignificantNodesList=new ArrayList<>();
        for(Node node:networkNodeList){
            if(nodeToZScoreMap.get(node)>0){
                significantNodesList.add(node);
            }else{
                nonsignificantNodesList.add(node);
            }
        }
//        System.out.println(""+significantNodesList.size()+" "+nonsignificantNodesList.size());
        
        Random random=new Random(Parameters.seedForRandom);
        

        for (int trial = 0; trial < numberOfTrials; trial++) {
//            long start=System.nanoTime();

            Collections.shuffle(nodeListForSampling, random);

            //These code can be used to first add significant nodes and start 
            //sampling with positive scored nodes
//            Collections.shuffle(significantNodesList, random);
//            Collections.shuffle(nonsignificantNodesList, random);
//            nodeListForSampling.clear();
//            nodeListForSampling.addAll(significantNodesList);
//            nodeListForSampling.addAll(nonsignificantNodesList);
            
            
            double zSum=0;
            int numberOfNodesInSubnetwork=0;
            for(Node node:nodeListForSampling){
                zSum=zSum+nodeToZScoreMap.get(node);
                numberOfNodesInSubnetwork++;
                double score=ScoreCalculations.this.calculateScoreOfSubnetwork(numberOfNodesInSubnetwork,zSum,false);
                samplingScoreSums[numberOfNodesInSubnetwork]+=score;
                samplingScoreSquareSums[numberOfNodesInSubnetwork]+=score*score;
                
                if(score<samplingScoreMins[numberOfNodesInSubnetwork]){
                    samplingScoreMins[numberOfNodesInSubnetwork]=score;
                }
                if(score>samplingScoreMaxs[numberOfNodesInSubnetwork]){
                    samplingScoreMaxs[numberOfNodesInSubnetwork]=score;
                }
                
            }
//            long stop=System.nanoTime();
//            System.out.println((stop-start)/1000);//ms
        }
        
        for(int i=1;i<=numberOfNodes;i++){
            samplingScoreMeans[i]=samplingScoreSums[i]/numberOfTrials;
            
            /**
             * var = SUM((x-xmean)^2) / N 
             * var = SUM(x^2 - 2*xmean*x + xmean^2)/N
             * var = SUM(x^2)/N - (2*xmean*SUM(x))/N + (N*xmean^2)/N 
             * var = SUM(x^2 )/N - 2*xmean^2 + xmean^2
             * var = SUM(x^2 )/N - xmean^2
             */
            samplingScoreStds[i]=samplingScoreSquareSums[i]/numberOfTrials - samplingScoreMeans[i]*samplingScoreMeans[i];
            samplingScoreStds[i]=Math.sqrt(samplingScoreStds[i]+0.0000001);
        }
        
    }
    
    /**
     * Calculates score of subnetwork. Returns zero for one node subnetworks.
     * @param subnetwork
     * @param subnetworkScoreNormalization
     * @return 
     */
    public double calculateScoreOfSubnetwork(Subnetwork subnetwork, boolean subnetworkScoreNormalization) {
        return ScoreCalculations.this.calculateScoreOfSubnetwork(subnetwork.getNodeList(), subnetworkScoreNormalization);
    }
    
    /**
     * Calculates score using node list. Returns zero for one node subnetworks.
     * @param nodeList
     * @param subnetworkScoreNormalization
     * @return 
     */
    public double calculateScoreOfSubnetwork(ArrayList<Node> nodeList, boolean subnetworkScoreNormalization) {
        int numberOfNodes=nodeList.size();
        double zSum=0;
        for(Node node:nodeList){
            zSum=zSum+nodeToZScoreMap.get(node);
        }
        return ScoreCalculations.this.calculateScoreOfSubnetwork(numberOfNodes, zSum, subnetworkScoreNormalization);
    }
    
    /**
     * Calculates score using z score sum and number of nodes. 
     * Returns zero for one node subnetworks.
     * @param numberOfNodes
     * @param zSum
     * @param subnetworkScoreNormalization
     * @return 
     */
    public double calculateScoreOfSubnetwork(int numberOfNodes, double zSum, boolean subnetworkScoreNormalization) {
        if(numberOfNodes==1){
            return 0;
        }
        double score=zSum/Math.sqrt(numberOfNodes);
        if(subnetworkScoreNormalization){
            score=normalizeScore(score, numberOfNodes);
        }
        if(Parameters.penaltyForSize){
            score=penaltyForSize(score, numberOfNodes);
        }
        return score;
    }
    
    private double normalizeScore(double score, int numberOfNodes){
        return (score-samplingScoreMeans[numberOfNodes])/samplingScoreStds[numberOfNodes];
    }
    
    private double penaltyForSize(double score, int numberOfNodes){
        score=score*Gaussian.cdf(score, 100, 30)*1000;
        return score;
    }
    

}
