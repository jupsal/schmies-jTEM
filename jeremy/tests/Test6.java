// This test is simply to write the 3rd derivative of LogTheta. We first try to
// figure out exactly how they are even taking the second derivative.

package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class Test6 {
    
    public static void main(String args[]) throws IOException{

        KPData data = example(0);
        data.printKPData();
        double x = 0.054;
        double y = 1.000;
        double t = 0;
        Complex soln;
        soln = data.KPSolutionAt( x, y, t );
        System.out.println(soln);
    }

    static private KPData example(int exampleNum) {
        Complex[] A;
        Complex[] B;
        Complex[] mu;
        int genus;
        KPData kdata;

        if( exampleNum == 0 ) {
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 4, 0 );
            B[0] = new Complex( -4, 0 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

