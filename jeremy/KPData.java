package jeremy.tests;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.Schottky;
import de.jtem.blas.*;
import de.jtem.riemann.theta.*;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.File;
import java.io.IOException;

public class KPData {

    private SchottkyData data; // Now we can mess with this.
    private Schottky schottky;
    private ComplexVector U, V, W, Z, T;
    private Complex c;
    private ComplexMatrix periodMatrix;
    private Theta theta;

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
        Z = new ComplexVector( U.size() );
        T = new ComplexVector( U.size() );
        setTheta( schottky );
    }

    public Theta getTheta() {
        return theta;
    }

    public void setTheta( Schottky sk ) {
        theta = new Theta( sk.getPeriodMatrix() );
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

    public Complex KPSolutionAt( double x, double y, double t ) {
        Z.assignTimes( U, x );
        T.assignTimes( V, y ); Z.assignPlus( T );
        T.assignTimes( W, t ); Z.assignPlus( T );
        Complex result = theta.ddLogTheta( Z, U, U );
        // This doesn't work rn because c cannot be evaluated
        // result.assinPlus(c);
        result.assignTimes(2);
        return result;
    }

    public Complex[][][] KPSolutionOnGrid( int numxsteps, int numysteps, int numtsteps, double T) {
        // Make a grid on (2pi)x(2pi) with gridspacing (points)x(points) number of points
        double x = 0.0;
        double y = 0.0;
        double t = 0.0;

        double maxX = 2*Math.PI; double maxY = 2*Math.PI;

        double deltax = 2*Math.PI/numxsteps;
        double deltay = 2*Math.PI/numysteps;
        double deltat = T/numtsteps;

        Complex[][][] result = new Complex[numxsteps][numysteps][numtsteps];

        for( int i = 0; i < numxsteps; i += 1 ) {
            y = 0.0;
            for( int j = 0; j < numysteps; j += 1 ) {
                t = 0.0;
                for( int k = 0; k < numtsteps; k += 1 ) {
                    result[i][j][k] = KPSolutionAt( x, y, t );
                    t += deltat;
                }
                y += deltay;
            }
            x += deltax;
        }
        return result;
    }


    // Below we try to create a CSV file using PrintWriter, as seen in 
    // http://stackoverflow.com/questions/22824523/writing-to-csv-in-multiple-columns-inside-a-loop-java
    //
    public void writeKPSolutionOnGrid1( int numxsteps, int numysteps, int numtsteps, double T, PrintWriter coordinateWriter, PrintWriter solnWriter) throws IOException
    {

        // throws IOException means we simply won't catch it if there are any
        // problems with IO.

        // Make a grid on (2pi)x(2pi) with gridspacing (points)x(points) number of points
        double x = 0.0;
        double y = 0.0;
        double t = 0.0;

        double maxX = 2*Math.PI; double maxY = 2*Math.PI;

        double deltax = 2*Math.PI/numxsteps;
        double deltay = 2*Math.PI/numysteps;
        double deltat = T/numtsteps;

        Complex result;

        coordinateWriter.print("t,");
        coordinateWriter.print("x,");
        coordinateWriter.print("y");
        coordinateWriter.print("\n");

        solnWriter.print("Real part, Imag part");
        solnWriter.print("\n");

        // Want (i,j,k) to index (x,y,t). But, want the loop to go the other
        // way.  i.e. want to look at a certain time-step first.
        for( int k = 0; k < numtsteps; k += 1 ) // t loop
        {
            x = 0.0;
            for( int i = 0; i < numxsteps; i += 1 ) // x loop
            {
                y = 0.0;
                for( int j = 0; j < numysteps; j += 1 ) // y loop
                {
                    // First write to the coordinate file.
                    if ( t < 0.0001 ) {
                        // Only do this at the initial time. Otherwise we just
                        // keep writing the same thing. 
                        coordinateWriter.print( String.format("%f,",t) );
                        coordinateWriter.print( String.format("%f,",x) );
                        coordinateWriter.print( String.format("%f",y) );
                        coordinateWriter.print( "\n" );
                    }

                    // Then calculate solution and write to KP file
                    result = KPSolutionAt( x, y, t );
                    solnWriter.print(String.format("%f,", result.getRe()) );
                    solnWriter.print(String.format("%f", result.getIm()) );
                    solnWriter.print("\n");

                    // Then update y and move on.
                    y += deltay;
                }
                // Update x at the end of the y loop
                x += deltax;
            }
            // Update t at the end of the x loop
            t += deltat;
        }
        coordinateWriter.flush();
        solnWriter.flush();
    }
}
