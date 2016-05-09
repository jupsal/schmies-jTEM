// This file is to test what different groups do to genus one solutions. Can
// they have arbitrary angle with the axes and how do we arrange for that?


package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class Test4 {
    
    public static void main(String args[]) throws IOException{

        KPData data = example(0);
        data.printKPData();
        Complex[][][] gridSoln;

        // The options below are to be changed.
        int numxsteps = 100; int numysteps = 100;
        int numtsteps = 300; double T = 10;

        // Write to file
        File coordFile = new File("/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/jeremy/coords3.csv");
        File solnFile = new File("/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/jeremy/soln3.csv");
        PrintWriter coordWriter = new PrintWriter( new FileWriter(coordFile) );
        PrintWriter solnWriter = new PrintWriter( new FileWriter(solnFile) );
        data.writeKPSolutionOnGrid1( numxsteps, numysteps, numtsteps, T, coordWriter, solnWriter );
        coordWriter.close();
        solnWriter.close();
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
            // This example should test how we can rotate the soln. The only
            // thing that is different than example 1 is that we increase the
            // REAL part of A, B. The imaginary parts are still the same.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 5, 1 );
            B[0] = new Complex( 5, -1 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 3 ) {
            // This example should test how we can rotate the soln. The only
            // thing that is different than example 1 is that we increase the
            // IMAGINARY part of A, B. The real parts are still the same.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 1, 5 );
            B[0] = new Complex( 1, -5 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 4 ) {
            // This example should test how we can rotate the soln. The only
            // thing that is different than example 1 is that we increase the
            // radii of the circles, i.e. mu, by two orders of magnitude.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 1, 1 );
            B[0] = new Complex( 1, -1 );
            mu[0] = new Complex( 1./100, 0 );
            kdata = new KPData(A, B, mu);
        }
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

