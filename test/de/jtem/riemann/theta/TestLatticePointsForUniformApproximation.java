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

import de.jtem.blas.RealMatrix;

public class TestLatticePointsForUniformApproximation {


    public static void main( String[] arg ) {

	/*
	RealMatrix B = new RealMatrix( 1 );
	
	B.set( 0, 0, 1 );

	LatticePointsInEllipsoid s = new LatticePointsInEllipsoidForTheta( B, 1.1  );

	RealMatrix field = s.getField();

	RealMatrix croppedField = new RealMatrix( field );

	croppedField.resize( s.getFieldLength(), 1 );

	System.out.println( croppedField );
	*/

	RealMatrix B = new RealMatrix ( 2 );

	B.set( 0,0, 0.9510565162 ); 
	B.set( 0,1, 0.3632712640 );

	B.set( 1,0, 0.3632712640 ); 
	B.set( 1,1, 0.9510565162 );

	// use the different normalization

	B.assignTimes( 2 * Math.PI );
    /*
	B.set( 0, 0, 1 );
	B.set( 1, 1, 1 );
    */
	LatticePointsForUniformApproximation_old s = new LatticePointsForUniformApproximation_old( B, 3.43, true );

	RealMatrix field = s.getLatticePoints();

	RealMatrix croppedField = new RealMatrix ( field );

	croppedField.resize( s.getNumOfLatticePoints(), 2 );

	System.out.println( croppedField );
	/*

	RealMatrix B = new RealMatrix( 3 );

	B.set( 0, 0, 1 );
	B.set( 1, 1, 1 );
	B.set( 2, 2, 1 );

	LatticePointsInEllipsoid s = new LatticePointsInEllipsoidForTheta( B, 1.1 );

	RealVector center = new RealVector( 3 );

	center.set( 1, 10 );

	s.setCenter( center );
	
	RealMatrix field = s.getField();

	RealMatrix croppedField = new RealMatrix( field );

	croppedField.resize( s.getFieldLength(), 3 );

	System.out.println( croppedField );

	s.check();
	*/
    }
}
















