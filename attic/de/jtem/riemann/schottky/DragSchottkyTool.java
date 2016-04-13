package de.jtem.riemann.schottky;


import java.awt.event.MouseEvent;
import java.util.Vector;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexValue;
import de.jtem.mfc.group.Moebius;
import de.jtem.moebiusViewer.AbstractPicker;
import de.jtem.moebiusViewer.tool.AbstractTool;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;

public class DragSchottkyTool extends AbstractTool
    implements  RealFunctionOfSeveralVariables {

    private static final long serialVersionUID = 1L;

    SchottkyShape selection;

    SchottkyData schottky;

    int index;

    char type;

    double alpha = 0;
    double beta  = 0;

    double radius;

    Complex A;
    Complex B;
    Complex C   = new Complex();
    Complex D   = new Complex();

    Complex tmp = new Complex();
    Complex mu  = new Complex();

    public DragSchottkyTool() {
	setLabel(  "drag fixpoints" );
    }

    public void mousePressed(MouseEvent e) {

	type = 'n';

	Vector dummy = new Vector();

	super.mousePressed( e );

	selection = (SchottkyShape)(getViewer().getSelection());

	schottky  = selection.getSchottky();

	AbstractPicker pickerContext = viewer.getPickerContext();

	pickerContext.setTransformStack( transformPath );

	pickerContext.initPicking ((int) e.getX (), (int) e.getY ());

	for( int i=0; i<schottky.numGenerators; i++ ) {

	    index = i;

	    A = schottky.fixpoint[i][0];
	    B = schottky.fixpoint[i][1];

	    pickerContext.point( A.re, A.im );

	    if( pickerContext.isPicked() ) {
		type = 'A';
		break;
	    }

	    pickerContext.point( B.re, B.im );

	    if( pickerContext.isPicked() ) {
		type = 'B';
		break;
	    }

	    selection.computeCAndD( i, C, D );

	    pickerContext.point( C.re, C.im );

	    if( pickerContext.isPicked() ) {
		tmp.assignMinus( D, B );
		beta = ComplexValue.argPositive(tmp);
		type = 'C';
		break;
	    }

	    pickerContext.point( D.re, D.im );

	    if( pickerContext.isPicked() ) {
		type = 'D';
		alpha = selection.alpha[i];
		break;
	    }
	}
    }

    Complex newPos = new Complex();

    void editSchottky() {

	if( type == 'n' )
	    return;

	newPos.assign( newPick.re, newPick.im );

	try {

	    switch( type ) {

	    case 'A':
		schottky.setA( index, newPos );
		break;

	    case 'B':
		schottky.setB( index, newPos );
		break;
	    case 'C':
		C.assign( newPos );
		tmp.assignMinus( C, A );
		alpha = selection.alpha[index] = ComplexValue.argPositive(tmp);

		radius = C.dist( schottky.center[index][0] );

		minimize();
		break;
	    case 'D':
		D.assign( newPos );
		tmp.assignMinus( D, B );
		beta = ComplexValue.argPositive(tmp);

		radius = D.dist( schottky.center[index][1] );

		minimize();
		break;
	    }

	    selection.fireMoebiusShapeChange();

	    getViewer().repaint();

	}
	catch( IllegalArgumentException exc ) {
	    java.awt.Toolkit.getDefaultToolkit().beep();
	}
    }

    public void mouseDragged(MouseEvent e) {

	super.mouseDragged( e );

	editSchottky();
    }



    public void mouseReleased(MouseEvent e) {

	super.mouseReleased( e );

	editSchottky();

	selection = null;
	schottky  = null;
	type = 'n';
    }


    public int getNumberOfVariables() {
	return 2;
    }

    static final double factor = 1e-3;

    public void setDoubleArrayParameter( double [] data, int offset ){

	schottky.getMu( index, mu );

	mu.re += data[ offset     ] * factor;
	mu.im += data[ offset + 1 ] * factor;
    }

    Moebius sigma = new Moebius();

    public double getDoubleValue() {
	double deltaArg;

	sigma.assign( A, B, mu );

	if( type == 'C' ) {

	    sigma.applyTo( C, D );
	    tmp.assignMinus( D, B );

	    deltaArg = Math.abs( beta - ComplexValue.argPositive(tmp) );


	} else if( type == 'D' ) {

	    sigma.applyInverseTo( D, C );
	    tmp.assignMinus( C, A );

	    deltaArg = Math.abs( alpha - ComplexValue.argPositive(tmp) );

	} else return 0;


	if( Math.abs( deltaArg - 2*Math.PI ) < deltaArg )
	    deltaArg -= 2*Math.PI;

	sigma.getC( tmp );

	double deltaAbs = radius - 1 / tmp.abs();

	return deltaArg * deltaArg + deltaAbs * deltaAbs;
    }

    public double eval( double [] p ) {
      setDoubleArrayParameter( p, 0 );
      return getDoubleValue();
    }

    void minimize() {

	double [] data = new double[2];

	de.jtem.numericalMethods.calculus.minimizing.Powell.search( data, 1e-9, this );

	setDoubleArrayParameter( data, 0 );

	sigma.getC( tmp );

//	SchottkyGroupElement tau = new SchottkyGroupElement();
//
//	tau.assign( sigma );
//
//	if( getDoubleValue() > 1e-9 || !schottky.isClassical() )
//	    throw new IllegalArgumentException();

	schottky.setMu( index, mu );

    }
}



























