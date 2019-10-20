package Application;

import ActiveSubnetworkSearchAlgorithms.ActiveSubnetworkSearch;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Ozan Ozisik
 */
public class AppActiveSubnetworkSearch {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try{
            processArguments(args);
        }catch(Exception e){
            Logger.getLogger(AppActiveSubnetworkSearch.class.getName()).log(Level.SEVERE, "Please check the arguments");
            System.exit(0);
        }        
        ActiveSubnetworkSearch.activeSubnetworkSearch();
        
    }
    
    public static void processArguments(String[] args) throws Exception {
        String helpText;
        helpText = "It is recommended to increase Stack size while running. You can use\n"
                + "java -Xss4m -jar ActiveSubnetworkSearch.jar\n"
                + "Options of the application are\n"
                + "-sif=<path>            \tuses the given interaction file\n"
                + "-sig=<path>            \tuses the given experiment file (gene p-value pairs)\n"
                + "-method=[GR|SA|GA]     \truns greedy search, simulated annealing or genetic algorithm for the search (default GR)\n"
                + "-useAllPositives       \tif used adds an individual with all positive nodes in GA, initializes candidate solution with all positive nodes in SA (default false)\n"
                + "-geneInitProb=<value>  \tprobability of adding a gene in inital solution for SA and GA (default 0.1)\n"
                + "-saTemp0=<value>       \tinitial temperature for SA (default 1.0)\n"
                + "-saTemp1=<value>       \tfinal temperature for SA (default 0.01)\n"
                + "-saIter=<value>        \titeration number for SA (default 10000)\n"
                + "-gaPop=<value>         \tpopulation size for GA (default 400)\n"
                + "-gaIter=<value>        \titeration number for GA (default 200)\n"
                + "-gaThread=<value>      \tnumber of threads to be used in GA (default 5)\n"
                + "-gaCrossover=<value>   \tapplies crossover with given probability (default 1)\n"
                + "-gaMut=<value>         \tapplies mutation with given rate (default 0)\n"
                + "-grMaxDepth=<value>    \tsets max depth in greedy search, 0 for no limit (default 1)\n"
                + "-grSearchDepth=<value> \tsets search depth in greedy search (default 1)\n"
                + "-grOverlap=<value>     \tsets overlap threshold for results of greedy search (default 0.5)\n"
                + "-grSubNum=<value>      \tsets number of subnetworks to be presented in the results (default 1000)\n";
        if(args.length==0 || args[0].equals("-h") || args[0].equals("help") || args[0].equals("-help")){
            System.out.println(helpText);
        }else{
            for(int i=0;i<args.length;i++){
                System.out.println(""+args[i]);
            }
            for(int i=0;i<args.length;i++){
                String[] str=args[i].split("=");
                String argType=str[0];
                String value="";
                if(str.length>1){
                    value=str[1];
                }
                
                switch(argType){
                    case "-sif":Parameters.sifPath=value;break;
                    case "-sig":Parameters.experimentFilePath=value;break;
                    case "-method":Parameters.useSAorGAorGR=Parameters.SearchMethod.valueOf(value);break;
                    case "-useAllPositives":Parameters.startWithAllPositiveZScoreNodes=true;break;
                    case "-geneInitProb":Parameters.geneInitialAdditionProbability=Double.parseDouble(value);break;
                    case "-saTemp0":Parameters.sa_initialTemperature=Double.parseDouble(value);break;
                    case "-saTemp1":Parameters.sa_finalTemperature=Double.parseDouble(value);break;
                    case "-saIter":Parameters.sa_totalIterations=Integer.parseInt(value);break;
                    case "-gaPop":Parameters.ga_populationSize=Integer.parseInt(value);break;
                    case "-gaIter":Parameters.ga_totalIterations=Integer.parseInt(value);break;
                    case "-gaThread":Parameters.ga_threadNumber=Integer.parseInt(value);break;
                    case "-gaMut":Parameters.ga_mutationRate=Double.parseDouble(value);break;
                    case "-grMaxDepth":Parameters.gr_maxDepth=Integer.parseInt(value);break;
                    case "-grSearchDepth":Parameters.gr_searchDepth=Integer.parseInt(value);break;
                    case "-grOverlap":Parameters.gr_overlapThreshold=Double.parseDouble(value);break;
                    case "-grSubNum":Parameters.gr_subnetworkNum=Integer.parseInt(value);break;
                    default:System.out.println("Unknown argument: "+argType);
                }
            }
        }
    }
}
