/**
 I don't know what is happening
 AHHH
 Why is this red?
 **/

//package nothing;
//

import de.jtem.mfc.field.Complex;

public class SchottkyGroupExample{

    public static void main (String args[]) {
        Complex[] A1 = new Complex[1];
        Complex[] A2 = new Complex[1];
        A1[0] = new Complex(-4,0);
        A2[0] = new Complex(-4,1);
        Complex[] mu = new Complex[1];
        mu[0] = new Complex(0.01,0);
        GetSchottkyData groupOne = new GetSchottkyData( A1, mu );
        GetSchottkyData groupTwo = new GetSchottkyData( A2, mu );

        // Print out what A and mu are
        groupOne.printTest();
        groupTwo.printTest();
    }
}
