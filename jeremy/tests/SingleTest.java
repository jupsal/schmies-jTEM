// This file is to test what different groups do to genus one solutions. In
// particular, this test looks at schottky groups with circles centered on the
// real axis.


package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class SingleTest {
    
    public static void main(String args[]) throws IOException{

        String localFileStructure = "/home/jeremy/Documents/schmies-jTEM/"; 

        for( int exampleNum = 0; exampleNum <= 0; exampleNum++)
        {
        //int exampleNum = 0;
        KPData data = example(exampleNum);
        data.printKPData();
        Complex[][][] gridSoln;

        // The options below are to be changed.
        int numxsteps = 200; int numysteps = 200;
        int numtsteps = 1; double T = 1;

        // Define write filename
        String coordStringName = localFileStructure+"jeremy/plotting/data/SingleTest/coords" + Integer.toString(exampleNum) + ".csv";
        String solnStringName = localFileStructure+"jeremy/plotting/data/SingleTest/soln" + Integer.toString(exampleNum) + ".csv";
        String groupStringName = localFileStructure+"jeremy/plotting/data/SingleTest/group" + Integer.toString(exampleNum) + ".csv";

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
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

