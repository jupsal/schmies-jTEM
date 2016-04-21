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

import java.util.Arrays;

/**
 * Multi-dimension monomial: f(x <sub>0</sub> ,...,x <sub>n</sub>-> = x
 * <sub>0</sub><sup>k <sub>0</sub></sup> x <sub>n</sub><sup>k <sub>
 * n</sub></sup>.
 */
public class Monomial
	implements BasisFunction {

	final int index[];
	final int order;
	final int numberOfVariables;
	final double coeff;
	
	Monomial(int [] index, int numberOfVariables ) {
		this( index, numberOfVariables, 1 );
	}
	Monomial(int [] index, int numberOfVariables, double coeff ) {
		this.index = (int[])index.clone();
		Arrays.sort( this.index );
		order = index.length;
		this.numberOfVariables = numberOfVariables;
		this.coeff = coeff;
	}

	/**
	 * @param k
	 *                     porwers for the different variables.
	 */
	private static int[] index(int[] k) {
		int order = 0;
		for (int i = 0; i < k.length; i++) {
			order += k[i];
		}
		int[] index = new int[order];
		for (int i = 0, l = 0; i < k.length; i++) {
			for (int j = 0; j < k[i]; j++, l++) {
				index[l] = i;
			}
		}
		return index;
	}

	static private int order(int[] k) {
		int order = 0;
		for (int i = 0; i < k.length; i++) {
			order += k[i];
		}
		return order;
	}

	public int order( int k ) {
		int order=0;
		for( int i=0; i<index.length; i++ )
			if(index[i] == k)
				order++;
			return order;
	}
	
	public BasisFunction getPartialDerivative( int k ) {
		int parcialOrder = order(k);
		if( parcialOrder == 0 )
			return new Monomial( new int[0], numberOfVariables, 0 );
		int [] indexOfDerivative = new int[ index.length -1 ];
		int i=0;
		for( ; index[i] != k; i++ )
			indexOfDerivative[i] = index[i];
		for( i++; i<index.length; i++ )
			indexOfDerivative[i-1] = index[i];
		return new Monomial( indexOfDerivative, numberOfVariables, coeff * parcialOrder );
	}
	
	public String toString()  {
		if( order == 0 )
			return new Double(coeff).toString();
		
		StringBuffer sb=new StringBuffer(300);
		
		if( coeff != 1 )
			sb.append( " " + coeff + " *" );
		
		sb.append(" x[" + index[0] + "]" );
		for( int i=1;i<index.length; i++ ) {
			sb.append(" * x[" + index[i] + "]" );
		}
		return sb.toString();
	}
	
	public double eval(double[] x) {
		double result = coeff;
		for (int i = 0; i < order; i++) {
			result *= x[i];
		}
		return result;
	}

	public int getNumberOfVariables() {
		return numberOfVariables;
	}

	public boolean equals(Monomial other) {
		if (numberOfVariables != other.numberOfVariables
			|| order != other.order)
			return false;
		for (int i = 0; i < order; i++) {
			if (index[i] != other.index[i])
				return false;
		}
		return true;
	}

	public boolean equals(Object other) {
		try {
			return equals((Monomial) other);
		} catch (ClassCastException e) {
			return false;
		}
	}

	/**
	 * Creates a multi dimension monomial.
	 * 
	 * @param k powers for the different variables.
	 */
	public static Monomial create(int[] k) {
		return new Monomial(index(k), k.length);
	}

	public static class Constant extends Monomial {
		public Constant(int numberOfVariables) {
			super(new int[0], numberOfVariables);
		}
		public double eval(double[] x) {
			return 1;
		}
	}

	public static class Linear extends Monomial {
		final int i0;
		public Linear(int i0, int numberOfVariables) {
			super(new int[] { i0 }, numberOfVariables);
			this.i0 = i0;
		}
		public double eval(double[] x) {
			return x[i0];
		}
	}

	public static class Quadratic extends Monomial {
		final int i0, i1;
		public Quadratic(int i0, int i1, int numberOfVariables) {
			super(new int[] { i0, i1 }, numberOfVariables);
			this.i0 = i0;
			this.i1 = i1;
		}
		public double eval(double[] x) {
			return x[i0] * x[i1];
		}
	}
	
	public static class Cubic extends Monomial {
		final int i0, i1, i2;
		public Cubic(int i0, int i1, int i2, int numberOfVariables) {
			super(new int[] { i0, i1, i2 }, numberOfVariables);
			this.i0 = i0;
			this.i1 = i1;
			this.i2 = i2;
		}
		public double eval(double[] x) {
			return x[i0] * x[i1] * x[i2];
		}
	}
	
	public boolean isZero() {
		return coeff == 0;
	}
}
