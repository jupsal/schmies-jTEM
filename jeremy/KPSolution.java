package jeremy;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.Schottky;
import de.jtem.blas.*;
import de.jtem.riemann.theta.*;


public class KP {
    
    Schottky schottky;
    Theta theta;
    ComplexVector U, V, W, Z, T;
    Complex c; 
    SchottkyData data;
    Schottky schottky = new Schottky( data );

    public KP( Complex [] A, Complex [] B, Complex [] mu ) {
        data = new SchottkyData( A.length );
        for( int i=0; i<A.length; i++ ) {
            data.setA( i, A[i] );
            data.setB( i, B[i] );
            data.setMu( i, mu[i] );
        }

        ComplexMatrix periodMatrix;
        periodMatrix = schottky.getPeriodMatrix();
        System.out.println(periodMatrix);

        U = schottky.getV();
        V = schottky.getV(2);
        W = schottky.getV(3);
        c = schottky.gamma();

        
        
        Z = new ComplexVector( U.size() );
        T = new ComplexVector( U.size() );

        Theta theta = new Theta( schottky.getPeriodMatrix() );
        // Now we are in a place that we can write down a KP solution
    }

    public Complex valueAt( double x, double y, double t ) {
        Z.assignTimes( U, x );
        T.assignTimes( V, y ); Z.assignPlus( T );
        T.assignTimes( W, t ); Z.assignPlus( T );
        Complex result = theta.ddLogTheta( Z, U, U );
        result.assignPlus(c);
        result.assignTimes(2);
        return result;
    }

    public PrintDiagnostics() {
        double[] radius;
        radius = data.getRadius();
        System.out.println("Radius="+radius[0]);
        Complex[][] centers;
        centers = data.getCenters();
        System.out.println("centers=" + centers[0][0]);
        System.out.println("Classical or not? = "+data.isClassical());
        System.out.println( "num Elements=" + schottky.getNumElements() );
        System.out.println( "Is this series evaluable? = "+schottky.isSeriesEvaluable() );
        System.out.println( "U =" + U );
        System.out.println( "V =" + V );
        System.out.println( "W =" + W );


    }
}



    

