// This file is to test what different groups do to genus one solutions. Can
// they have arbitrary angle with the axes and how do we arrange for that?


package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class Test5 {
    
    public static void main(String args[]) throws IOException{

        // Choose the appropriate local filestructure putting you at the
        // mainlevel of the git repo
        //
        // File structure for desktop in office. 
        //String localFileStructure = "/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/"; 
        // File structure for laptop. 
        String localFileStructure = "/home/jeremy/Documents/schmies-JTEM/"; 

        int numExamples = 6;
        for( int exampleNum = 0; exampleNum <= numExamples; exampleNum++)
        {
        //int exampleNum = 0;
        KPData data = example(exampleNum);
        data.printKPData();
        Complex[][][] gridSoln;

        // The options below are to be changed.
        int numxsteps = 200; int numysteps = 200;
        int numtsteps = 1; double T = 1;

        // Define write filename
        String coordStringName = localFileStructure + "jeremy/plotting/data/Test5/coords" + Integer.toString(exampleNum) + ".csv";
        String solnStringName = localFileStructure + "jeremy/plotting/data/Test5/soln" + Integer.toString(exampleNum) + ".csv";
        String groupStringName = localFileStructure + "jeremy/plotting/data/Test5/group" + Integer.toString(exampleNum) + ".csv";

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
            // Start with a group which is off both axes in an appropriate way.
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 1, 1 );
            B[0] = new Complex( 1, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 1 ) {
            // Increase the real part a little bit
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 2, 1 );
            B[0] = new Complex( 2, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 2 ) {
            // Increase the real part more
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( 4, 1 );
            B[0] = new Complex( 4, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 3 ) {
            // Flip example0 about the imaginary axis
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( -1, 1 );
            B[0] = new Complex( -1, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 4 ) {
            // Decrease the real part more
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( -2, 1 );
            B[0] = new Complex( -2, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 5 ) {
            // Decrease the real part even
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( -4, 1 );
            B[0] = new Complex( -4, -1 );
            mu[0] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else if( exampleNum == 6 ) {
            // Change the imaginary part
            genus = 1;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( -0.2, 1 );
            B[0] = new Complex( -0.2, -1 );
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

