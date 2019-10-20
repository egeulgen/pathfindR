package File;

import Network.Network;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Ozan Ozisik
 */
public class SIFReader {
    
    public static Network readSIF(String path){
    
        try {
            int columnNumber;
            
            BufferedReader bufReader=new BufferedReader(new FileReader(path));
            Network network=new Network();
            
            bufReader.mark(300);
            String line;
            line = bufReader.readLine();
            String[] strArr=line.split("[ \\t]");
            if(strArr.length==2){
                columnNumber=2;
            }else if(strArr.length==3){
                columnNumber=3;
            }else{
                Logger.getLogger(SIFReader.class.getName()).log(Level.SEVERE, "SIF file must have 2 or 3 columns");
                return null;
            }
            
            bufReader.reset();
            
            int lineNo=1;
            while ((line = bufReader.readLine()) != null) {
                strArr=line.split("[ \\t]");
                if(strArr.length==columnNumber){
                    String strNode1, strNode2;
                    if(columnNumber==2){
                        strNode1=strArr[0];
                        strNode2=strArr[1];
                    }else{
                        strNode1=strArr[0];
                        strNode2=strArr[2];
                    }
                    strNode1=strNode1.toUpperCase();
                    strNode2=strNode2.toUpperCase();
                    network.addInteraction(strNode1,strNode2);
                }else{
                    Logger.getLogger(SIFReader.class.getName()).log(Level.WARNING, "Unexpected column number in SIF line {0}, discarded", lineNo);
                }
                lineNo++;
            }
            return network;
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(SIFReader.class.getName()).log(Level.SEVERE, "SIF file not found", ex);
            return null;
        } catch (IOException ex) {
            Logger.getLogger(SIFReader.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
        
    }
}
