// This file is to test what different groups do to genus one solutions. In
// particular, this test looks at schottky groups which give rise to solutions
// that are constant wrt y.


package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class Test4 {
    
    public static void main(String args[]) throws IOException{

        // Choose the appropriate local filestructure putting you at the
        // mainlevel of the git repo
        //
        // File structure for desktop in office. 
        String localFileStructure = "/home/jeremy/Documents/schmies-jTEM/"; 
        //String localFileStructure = "/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/"; 

        int numExamples = 6;
        for( int exampleNum = 0; exampleNum <= numExamples; exampleNum++)
        {
        KPData data = example(exampleNum);
        data.printKPData();
        Complex[][][] gridSoln;

        // The options below are to be changed.
        int numxsteps = 200; int numysteps = 200;
        int numtsteps = 1; double T = 1;

        // Define write filename
        String coordStringName = localFileStructure + "jeremy/plotting/data/Test4/coords" + Integer.toString(exampleNum) + ".csv";
        String solnStringName = localFileStructure + "jeremy/plotting/data/Test4/soln" + Integer.toString(exampleNum) + ".csv";
        String groupStringName = localFileStructure + "jeremy/plotting/data/Test4/group" + Integer.toString(exampleNum) + ".csv";

        // Create write file
        File coordFile = new File(coordStringName);
        File solnFile = new File(solnStringName);
        File groupFile = new File(groupStringName);

        // Define the writers
        PrintWriter coordWriter = new PrintWriter( new FileWriter(coordFile) );
        PrintWriter solnWriter = new PrintWriter( new FileWriter(solnFile) );
        PrintWriter groupWriter = new PrintWriter( new FileWriter(groupFile) );

        // Create KP data and write
        data.writeKPSolutionOnGrid1( numxsteps, numysteps, numtsteps, T, coordWriter, solnWriter );
        coordWriter.close();
        solnWriter.close();

        // Write group data
        data.writeGroupData( groupWriter );
        groupWriter.close();
        }
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
            // This one just changes the x-coordinate of A and B from example0
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 1, 0 );
            B[0] = new Complex( -1, 0 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 2 ) {
            // This one just changes mu from example0. By a factor of 1000.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 4, 0 );
            B[0] = new Complex( -4, 0 );
            mu[0] = new Complex( 1./10, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 3 ) {
            // This one is a whole new example. Now we have a group which has
            // centers on the imaginary axis and is symmetric wrt to the real
            // axis. You get the same thing! 1D solns constant wrt y
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 0, 4 ); // Interesting. These are also constant wrt y.
            B[0] = new Complex( 0, -4 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 4 ) {
            // This example should test how we can rotate the soln. The only
            // thing that is different than example 4 is that we decrease the
            // IMAGINARY part of A, B. The real parts are still the same.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 0, 1 );
            B[0] = new Complex( 0, -1 );
            mu[0] = new Complex( 1./10000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 5 ) {
            // This example should test how we can rotate the soln. The only
            // thing that is different than example 4 is that we increase the
            // radii of the circles, i.e. mu, by two orders of magnitude.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 0, 4 );
            B[0] = new Complex( 0, -4 );
            mu[0] = new Complex( 1./10, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 6 ) {
            // This example should help to determine, for stuff on the imaginary
            // axis, does period increase with larger radius circles or further
            // from imaginary axis?
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 0, 1 );
            B[0] = new Complex( 0, -1 );
            mu[0] = new Complex( 1./4, 0 );
            kdata = new KPData(A, B, mu);
        }
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

