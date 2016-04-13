/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2009, Technische Universität Berlin, jTEM
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

-	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

-	Redistributions in binary form must reproduce the above copyright notice, 
	this list of conditions and the following disclaimer in the documentation 
	and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.
**/

package de.jtem.riemann.theta;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

public class TestSiegelTransform {

    public static void main( String [] argb ) {

	int n = 2;

	ComplexMatrix B = new ComplexMatrix( 2 );

	B.set( 0, 0, -5.10972365633887, Math.PI );
	B.set( 0, 1, -4.24199777781055, 0       );

	B.set( 1, 0, -4.24199777781055, 0       );
	B.set( 1, 1, -5.78346380443502, Math.PI );

	B.print( "B" );

	/*
	IntegerMatrix A = TestModularGroup.getRandomInvertableIntegerMatrix(2);

	A.print( "A" );

	RealMatrix DA = new RealMatrix().set( A );
	ComplexMatrix CA = new ComplexMatrix(n).setRe( DA );

	Theta theta            = new Theta( B );
	Theta transformedTheta = new Theta( A.times( B ).assignTimes( A.transpose() ) );


	theta.setTol( 1e-16 );
	transformedTheta.setTol( 1e-16 );

	ComplexVector Z = new ComplexVector( n );


	for( int i=0; i<5; i++ ) {
	    Z.print(" ______ Z = " );
	    System.out.println( " trans theta: "
				+ transformedTheta.theta( Z ) );


	    System.out.println( "       theta: "
				+ theta.theta( CA.invert().times( Z )) );

	    Z.assignRandom().assignTimes( 3 );
	}
	*/

	/*
	IntegerMatrix b = TestModularGroup.getRandomSymmetricIntegerMatrix(2);

	b.print( "b" );

	ComplexMatrix cb = b.times( new Field.Complex( 0, 2*Math.PI ) );

	ComplexVector cv = cb.getDiagonal().assignDivide( 2 );

	Theta theta            = new Theta( B );
	Theta transformedTheta = new Theta( B.assignPlus( cb ) );


	theta.setTol( 1e-16 );
	transformedTheta.setTol( 1e-16 );

	ComplexVector Z = new ComplexVector( n );

	Field.Complex t = new Field.Complex();
	Field.Complex f = new Field.Complex();

	for( int i=0; i<5; i++ ) {
	    Z.print(" ______ Z = " );
	    transformedTheta.theta( Z, f, t );
	    System.out.println( " trans theta: " + t.assignTimes( f.exp() ) );
	    theta.theta( cv.plus( Z ), f, t );
	    System.out.println( "       theta: " + t.assignTimes( f.exp() ) );

	    Z.assignRandom().assignTimes( 3 );
	}
	*/

	ComplexMatrix invOfB = B.invert();


	Theta theta            = new Theta( B );
	Theta transformedTheta = new Theta( invOfB.times( 4*Math.PI*Math.PI  ) );

	theta.setTol( 1e-16 );
	transformedTheta.setTol( 1e-16 );

	ComplexVector Z = new ComplexVector( n );

	Complex t1 = new Complex();
	Complex f1 = new Complex();
	Complex t2 = new Complex();
	Complex f2 = new Complex();
	for( int i=0; i<5; i++ ) {

	    transformedTheta.thetaSum( invOfB.times( Z ).times( new Complex( 0, - 2*Math.PI ) ), t1 );
	    theta.           thetaSum( Z, t2 );

	    f2.assignPlus( Z.dotBilinear( invOfB.times( Z ) ).divide( 2 ),
			  B.divide( 2*Math.PI ).determinant().log().divide( 2 ) );

	    System.out.println( " trans theta: " + t1 );
	    System.out.println( "       theta: " + t2.times( f2.exp() ) );
	    System.out.println( "       theta: " + t1.divide( t2 ) );
	    /*
	    Z.print(" ______ Z = " );
	    System.out.println( " trans theta: " + t1.assignTimes( f1.exp() ) );
	    Z.print(" ______ Z = " );
	    transformedTheta._theta( invOfB.times( Z ).assignTimes( new Field.Complex( 0, -2*Math.PI ) ), f1, t1 );
	    System.out.println( " trans theta: " + t1.assignTimes( f1.exp() ) );
	    theta._theta( Z, f2, t2 );
	    f2.assignPlus( Z.times( invOfB.times( Z ) ).assignTimes( 0.5 ) );
	    t2.assignTimes( f2.exp() ).assignDivide( B.divide( 2*Math.PI ).determinant() );
	    System.out.println( "       theta: " + t2 );;
	    */
	    Z.assignRandom();
	    Z.assignTimes( 3 );
	}

    }
}
















