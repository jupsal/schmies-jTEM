package jeremy.tests;

import de.jtem.mfc.field.Complex;
import java.io.*;

public class Test3 {
    
    public static void main(String args[]) throws IOException{

        KPData data = example(3);
        data.printKPData();
        double x = 0.001;
        double y = 1.000;
        double t = 10;
        Complex soln;
        soln = data.KPSolutionAt( x, y, t );
        System.out.println(soln);
        Complex[][][] gridSoln;
        int numxsteps = 100; int numysteps = 100;
        int numtsteps = 300; double T = 10;
        //String filename = "/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/filename.txt";
        gridSoln = data.KPSolutionOnGrid(numxsteps,numysteps,numtsteps,T);

        // Try to do it with writing now.
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
        else if( exampleNum == 3 ) {
            // This is a genus 2 example of a KP2 solution as seen on p.163 of
            // Algebro-geometric approach. It is a UII uniformization.
            genus = 2;
            A = new Complex[genus];
            B = new Complex[genus];
            mu = new Complex[genus];
            A[0] = new Complex( -4, 0 );
            B[0] = new Complex( 4, 0 );
            mu[0] = new Complex( 1./10000, 0 );
            A[1] = new Complex( -1, 0 );
            B[1] = new Complex( 1, 0 );
            mu[1] = new Complex( 1./1000, 0 );
            kdata = new KPData(A, B, mu);
        }
        else {
            System.out.println( "Try a different example number, you used example "+exampleNum+" which does not exist." );
            throw new RuntimeException("Example does not exist");
        }
        return kdata;
    }
}

