/**
 I don't know what is happening
 AHHH
 Why is this red?
 **/

//package nothing;

import de.jtem.blas.*;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.*;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.theta.*;
//import junit.framework.TestCase;


class justASchottkyGroup {
    public static void main(String[] args) {
        new ComplexMatrix( new double [][]
                {{ -3.6 }},
                        new double [][]
                {{ 3.14}} );
        RealMatrix Bmat = new RealMatrix( 2 );
        Bmat.set(0,0,1);
        Bmat.set(0,1,0);
        Bmat.set(1,0,0);
        Bmat.set(1,1,1);
        System.out.println( Bmat );
        System.out.println("hello World!");

        Schottky schottky; // Gives 'schottky' \in class(Schottky)
        Complex[] A = new Complex[2];
        A[0] = new Complex( 4.0,0 );
        A[1] = new Complex( -4.0,0 );
        System.out.println(A[0]);
        System.out.println(A[1]);
        System.out.println( A.length );
        SchottkyData data;
        // SchottkyData data = new SchottkyData( A.length );
       // Complex A = new Complex(4.0,0);
       // Complex B = new Complex(-4.0,0);
       // Complex mu = new Complex(1./10000,0);
       // data.setA( 0, A );
       // data.setB( 0, B );
       // data.setMu( 0, mu);
       // if( !data.isClassical() )
       //     throw new IllegalArgumentException
       //           ( "schottky data is not classical");  

        //element.assign( new Complex( 3.0, -5.0),
        //                new Complex( -4.0, 1.0),
        //                new Complex( 0.0, 1.0) );
        //System.out.println(element);
    }
}
