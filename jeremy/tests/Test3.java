package jeremy.tests;

import de.jtem.mfc.field.Complex;

public class Test3 {
    
    public static void main(String args[]) {

        KPData data = example(0);
        data.printKPData();

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
        else if( exampleNum == 1 ) {
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 1, 1 );
            B[0] = new Complex( 1, -1 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 2 ) {
            // This is a genus 2 example of a KP2 solution as seen on p.163 of
            // Algebro-geometric approach. It is a UII uniformization.
            genus = 2;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            Complex a0 = new Complex( -0.024, 0.507 );
            Complex a1 = new Complex( 0.588, 0.556 );
            A[0] = a0; B[0] = a0.conjugate();
            A[1] = a1; B[1] = a1.conjugate();
            System.out.println( "A[0]="+A[0]);
            System.out.println( "A[1]="+A[1]);
            System.out.println( "B[0]="+B[0]);
            System.out.println( "B[1]="+B[1]);
            mu[0] = new Complex( 0.076, 0 );
            mu[1] = new Complex( 0.014, 0 );
            kdata = new KPData(A, B, mu);
        }
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

