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

package de.jtem.numericalMethods.calculus.functionApproximation.bestFitting;

import java.util.List;

/**
 * Storage class for a multidimensional weighted sample.
 */
public class Sample implements Cloneable {

		final double [] x;
		final double [] y;
		
		final double weight;
		
		/**
		 * Creates sample with provided weight.
		 * @param x
		 * @param y
		 * @param weight
		 */
		public Sample( double [] x, double [] y, double weight ) {
			this.x = (double[])x.clone();
			this.y = (double[])y.clone();
			this.weight = weight;
		}
		
		/**
		 * Creates sample with weight equalt to 1.
		 * @param x
		 * @param y
		 */
		public Sample( double [] x, double [] y ) {
			this( x,y,1);
		}
		public boolean sameDimensions( Sample other ) {
			return other.x.length == x.length && other.y.length == y.length;
		}
			
		/**
		 * Converts a list of samples into an array.
		 * @param sampleList list of samples
		 * @return array of samples.
		 */
		public static Sample [] toArray( List sampleList ) {
			final int n =sampleList.size();
			Sample [] array = new Sample[n];
			for( int i=0; i<n; i++ ) {
				array[i] = (Sample)sampleList.get(i);
			}
			return array;
		}
		
		public static Sample [] toArray( double [][] x, double [][] y ) {
			if( x.length != y.length )
				throw new IllegalArgumentException( "list length do not match");
			final int n = x.length;
			Sample [] array = new Sample[n];
			for( int i=0; i<n; i++ ) {
				array[i] = new Sample( x[i], y[i] );
			}
			return array;
			
		}
}