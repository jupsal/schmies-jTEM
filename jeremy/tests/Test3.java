package jeremy.tests;

import de.jtem.mfc.field.Complex;

public class Test3 {
    
    public static void main(String args[]) {

        KPData data = example(0);
        data.printKPData();

    }

    static private KPData example(int exampleNum) {
        int num = 1;
        Complex[] A = new Complex[num];
        Complex[] B = new Complex[num];
        Complex[] mu = new Complex[num];
        KPData kdata;

        if( exampleNum == 0 ) {
            A[0] = new Complex( 4, 0 );
            B[0] = new Complex( -4, 0 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 1 ) {
            A[0] = new Complex( 1, 1 );
            B[0] = new Complex( 1, -1 );
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

