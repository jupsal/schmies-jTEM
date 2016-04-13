import de.jtem.blas.*;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.*;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.theta.*;

public class GetSchottkyData {
    
    SchottkyData schottky;
    Complex[] A;
    Complex[] mu;

    public GetSchottkyData( Complex [] Ain, Complex [] muin ) {
        A = Ain;
        mu = muin;
        SchottkyData data = new SchottkyData( A.length );
    }
    public void printTest() {
        System.out.println("Length(A):" + A.length);
        for( int i=0; i<A.length; i++ ) {
            System.out.println("A[i]:"+A[i]);
            System.out.println("mu[i]:"+mu[i]);
        }
    }
}

