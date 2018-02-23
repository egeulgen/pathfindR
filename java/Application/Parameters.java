package Application;

/**
 *
 * @author Ozan Ozisik
 */


public class Parameters {
    public static String sifPath="BIOGRID-ORGANISM-Homo_sapiens-3.4.155.OzCleaned.sif";
    public static String experimentFilePath="Behcet_jp_GWASPvalue.txt";
    public static String resultFilePath="resultActiveSubnetworkSearch.txt";
    
    public enum SearchMethod{GR, SA, GA};
    public static SearchMethod useSAorGAorGR=SearchMethod.GR;//(default GR)
    public static boolean startWithAllPositiveZScoreNodes=false;//(default false)
    public static double geneInitialAdditionProbability=0.1;//(default 0.1)
    
    public static boolean penaltyForSize=false;
    public static double pForNonSignificantNodes=0.5;//0.9999999999999
    
    public static double sa_initialTemperature=1.0;//(default 1.0)
    public static double sa_finalTemperature=0.01;//(default 0.01)
    public static int sa_totalIterations=10000;//(default 10000)
    
    
    public static int ga_populationSize=400;//(default 400)
    public static int ga_totalIterations=200;//(default 200)
    public static int ga_threadNumber=5;//(default 5)
    public static boolean ga_mutation=false;//(default mutation off)
    public static double ga_mutationRate=0.05;
    public static boolean ga_Elitism=true;
    
    public static int gr_maxDepth=1;//(default 1)
    public static int gr_searchDepth=1;//(default 1)
    public static double gr_overlapThreshold=0.5;//(default 0.5)
    public static double gr_subnetworkNum=1000;//(default 1000)
    
}