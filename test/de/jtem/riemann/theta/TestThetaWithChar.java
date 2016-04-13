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

class TestThetaWithChar {

    public static void main( String[] arg ) {

	ComplexMatrix B = new ComplexMatrix( 2 );

	B.set( 0, 0, new Complex( -3.996246, 0) );
	B.set( 1, 0, new Complex( -1.767355, 3.14159) );
	B.set( 0, 1, new Complex( -1.767355, 3.14159) );
	B.set( 1, 1, new Complex( -3.996246, 0) );

	ThetaWithChar th = new ThetaWithChar( B );

	
	System.out.println( "Radius   = " + th.theta.getRadius() );
	System.out.println( "numOfLatticePoints = " + th.theta.getNumOfLatticePoints() );

	System.out.println( "B =" + B );

	//th.getLatticePointsInEllipsoid().print();

	System.out.println( "==================================" );

	Complex result = new Complex();
	Complex th_z   = new Complex();
	Complex th_v   = new Complex();
	Complex factor = new Complex();

	ComplexVector z =  new ComplexVector( 2 );
	ComplexVector v =  new ComplexVector( 2 );

	z.set( 0, new Complex( -2.0104, 4.7123 ) );
	z.set( 1, new Complex( -2.0104, 1.5707 ) );

	v.set( 0, new Complex( 0,  0.3266 ) );
	v.set( 1, new Complex( 0, -0.3266 ) );

	th.dTheta( z, v, factor, th_z, th_v );

	System.out.println( "th_z    = " + th_z );
	System.out.println( "th_v    = " + th_v );
	System.out.println( "factor  = " + factor );
	System.out.println( "value z = " + Complex.exp(factor).times( th_z ) );
	System.out.println( "value v = " + Complex.exp(factor).times( th_v ) );

    }
}
















