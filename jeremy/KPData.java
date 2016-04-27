package jeremy.tests;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.Schottky;
import de.jtem.blas.*;
import de.jtem.riemann.theta.*;

public class KPData {

    SchottkyData data; // Now we can mess with this.
    Schottky schottky;
    ComplexVector U, V, W, Z, T;
    Complex c;
    ComplexMatrix periodMatrix;

    public KPData(Complex[] A, Complex[] B, Complex[] mu) {
        data = new SchottkyData( A.length );
        for( int i=0; i<A.length; i++ ) {
            data.setA( i, A[i] );
            data.setB( i, B[i] );
            data.setMu( i, mu[i] );
        }
        System.out.println(data.getUniformizationData());
        schottky = new Schottky( data );
        System.out.println("Number of elements = " + schottky.getNumElements());
        printNumElements();
        System.out.println("Is differential series evaluable? " + schottky.isDifferentialSeriesEvaluable()); // THIS ALSO CAUSES THE BUG
        printNumElements();
        System.out.println("Is Integral series evaluable? " + schottky.isIntegralSeriesEvaluable()); // THIS CAUSES THE BUG. After this numelements increases.
        printNumElements();
        printseriesEvaluableOrNot();
        printNumElements();
        // periodMatrix = schottky.getPeriodMatrix(); // Don't need this now?
        U = schottky.getV();    // U
        V = schottky.getV(2);   // V
        W = schottky.getV(3);   // W, from \theta(Ux + Vy + tW) in KP Soln
        //printKPVectors();
        //THIS DOESN'T WORK!! //c = schottky.gamma();   // Constant at end of formula for KP soln
        //System.out.println( "WTF?" + schottky.isSeriesEvaluable() );
    }

    public void printKPData() {
        printRadii();
        printCenters();
        printClassicallness();
        printNumElements();
        //printseriesEvaluableOrNot();
        printNumElements();
        printKPVectors();
    }

    public void printRadii() {
        double[] radius;
        radius = data.getRadius();
        for( int i=0; i<radius.length; i++ ) {
            System.out.println( "Circle-"+i+" has radius = "+radius[i] );
        }
    }

    public void printCenters() {
        Complex[][] centers;
        centers = data.getCenters();
        for( int i=0; i<centers.length; i++ ) {
            System.out.println( "Circle-"+i+" has center = "+centers[i][0] );
            System.out.println( "Circle-"+i+"' has center = "+centers[i][1] );
        }
    }

    public void printClassicallness() {
        System.out.println( "This series is classical? (T/F) : " + data.isClassical() );
    }

    public void printNumElements() {
        System.out.println( "The number of elements in schottky group = " + schottky.getNumElements() );
    }

    public void printseriesEvaluableOrNot() {
        System.out.println( "Can we evaluate the series from the Schottky group? (T/F) : " + schottky.isSeriesEvaluable() );
    }

    public void printKPVectors() {
        System.out.println( "U = " + U );
        System.out.println( "V = " + V );
        System.out.println( "W = " + W );
        //System.out.println( "c = " + c );
    }

}
