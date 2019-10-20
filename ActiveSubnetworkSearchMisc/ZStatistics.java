package ActiveSubnetworkSearchMisc;

/**
 *
 * @author Ozan Ozisik
 * adapted from  https://github.com/idekerlab/jActiveModules
 */
public class ZStatistics {
    
    public static double oneMinusNormalCDFInverse(double p) {
        if (p <= 0.5) {
            if (p > 0) {
                return oneMinusNormalCDFInversePLT5(p);
            } else {
                return Double.POSITIVE_INFINITY;
            }
        } else if (p < 1) {
            return -oneMinusNormalCDFInversePLT5(1 - p);
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }
    
    //from 26.2.23, page 933, Handbook of Mathematical Functions, NBS, 1964
    //Requires 0 < p <= 0.5
    private static double oneMinusNormalCDFInversePLT5(double p) {

        double t, temp;

        if (p < 0) {
            throw new IllegalArgumentException("oneMinusNormalCDFInversePLT5 called with negative p\n");
        } else if (p > 0.5) {
            throw new IllegalArgumentException("oneMinusNormalCDFInversePLT5 called with p > 0.5\n");
        } else {
            t = Math.sqrt(-2 * Math.log(p));
            temp = 2.515517 + 0.802853 * t + 0.010328 * t * t;
            temp = t - temp / (1 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);
            return temp;
        }
    }
    
    
}
