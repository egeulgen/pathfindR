package File;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;

/**
 *
 * @author Ozan Ozisik
 */
public class ExperimentFileReader {
    
    public static ArrayList<SimpleEntry<String, Double>> readExperimentFile(String path){
        
        try {
            //A list of pairs is used to allow multiple p-values for the same gene.
            //These multiple p-values situation is not business of this class
            ArrayList<SimpleEntry<String, Double>> namePValuePairList=new ArrayList<>();
            
            BufferedReader bufReader=new BufferedReader(new FileReader(path));
                        
            String line;
            String[] strArr;
            int lineNo=1;
            while ((line = bufReader.readLine()) != null) {
                strArr=line.split("[ \\t]");
                if(strArr.length==2){
                    try{
                        namePValuePairList.add(new SimpleEntry<>(strArr[0], Double.parseDouble(strArr[1])));
                    }catch(NumberFormatException nfe){
                        Logger.getLogger(ExperimentFileReader.class.getName()).log(Level.WARNING, "Unexpected number format in experiment file line {0}, discarded", lineNo);
                    }
                }else{
                    Logger.getLogger(ExperimentFileReader.class.getName()).log(Level.WARNING, "Unexpected column number in experiment file line {0}, discarded", lineNo);
                }
                lineNo++;
            }
            return namePValuePairList;
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ExperimentFileReader.class.getName()).log(Level.SEVERE, "Experiment file not found", ex);
            return null;
        } catch (IOException ex) {
            Logger.getLogger(ExperimentFileReader.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
        
    }
            
}
