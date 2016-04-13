import de.jtem.blas.*;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.*;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.theta.*;

public class GetSchottkyData {
    
    SchottkyData schottky;

    public GetSchottkyData( SchottkyData schottky) {
        setSchottky(schottky);
    }

    public void setSchottky(SchottkyData v) {
        schottky = v;
        
    // Construct group?
    public GetSchottkyData( Complex [] A, Complex [] mu ) {
        SchottkyData data = new SchottkyData( A.length )
