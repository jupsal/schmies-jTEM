import mfc.field.Complex;
import riemann.schottky.*;
import riemann.theta.*;
import blas.*;

public class KP2 {
    Schottky schottky;
    Theta theta;
    ComplexVector U, V, W, Z, T;
    Complex c;
    public KP2( Complex [] A, double [] mu, double eps ) {
        SchottkyData data = new SchottkyData( A.length );
        for( int i=0; i<A.length; i++ ) {
            data.setA( i, A[i] );
            data.setB( i, A[i].conjugate() );
            data.setMu( i, mu[i], 0 );
        }
        if( !schottkyData.isClassical() )
            throw new IllegalArgumentException
                ( "schottky data is not classical");
        schottky = new Schottky( schottkyData, eps );
