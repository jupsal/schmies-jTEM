package jeremy.tests;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.Schottky;


public class Test1 {
    public static void main(String args[]) {

        Complex comp = new Complex( 4, 0 );
        Complex [] A = new Complex[1];
        A[0] = comp;
        Complex B = new Complex( -4, 0 );
        System.out.println("A[0]="+A[0]);
        System.out.println("B[0]="+B);
        Complex [] mu = new Complex[1];
        mu[0] = new Complex(1./10000, 0);
        SchottkyData data = new SchottkyData( A.length );
        data.setA( 0, A[0] );
        data.setB( 0, B );
        data.setMu( 0, mu[0] );
        double[] radius;
        radius = data.getRadius();
        System.out.println("Radius="+radius[0]);
        Complex[][] centers;
        centers = data.getCenters();
        System.out.println("centers=" + centers[0][0]);
        System.out.println("Classical or not? = "+data.isClassical());

        Schottky schottky = new Schottky( data );
        System.out.println( "num Elements=" + schottky.getNumElements() );
        System.out.println( "Is this series evaluable? = "+schottky.isSeriesEvaluable() );
        

    }
}



    

