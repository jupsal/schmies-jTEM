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

package de.jtem.numericalMethods.calculus.rootFinding;

import de.jtem.numericalMethods.algebra.linear.MatrixOperations;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;
import de.jtem.numericalMethods.algebra.linear.decompose.Householder;
import de.jtem.numericalMethods.algebra.linear.solve.RXB;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfSeveralVariablesWithJacobien;

/**
 * @author schmies
 */
public class Broyden {

	final static double DEFAULT_MAX_STEPSIZE = 1;
	final static double EPS = 1e-15;
	
	Broyden() {
	}

	public static class HitLocalMinimumException extends RuntimeException {
		public HitLocalMinimumException() {
			super();
		}
		public HitLocalMinimumException( String msg ) {
			super();
		}
	}
	
	public static void search( RealVectorValuedFunctionOfSeveralVariablesWithJacobien func,
			double x[] )
	{
		search( func, x, 200, 1e-8, 1e-10, 1e-15,  DEFAULT_MAX_STEPSIZE);
	}
	
	public static void search( RealVectorValuedFunctionOfSeveralVariablesWithJacobien func,
			double x[], int maxNumberOfIterations, 
			double acc, double precision  )
	{
		search( func, x, maxNumberOfIterations, acc, acc/100, precision, DEFAULT_MAX_STEPSIZE );
	}
	
	public static void search( RealVectorValuedFunctionOfSeveralVariablesWithJacobien func,
			double x[], int maxNumberOfIterations, 
			double acc, double accOfMin, double precision, double maxStepSize  )
	{
		
		final int n = func.getDimensionOfTargetSpace();
		

		double[] c=new double[n];
		double[] fvcold=new double[n];
		double[] g=new double[n];
		double[] p=new double[n];
		double[][] ts=new double[n][n];
		double[][] r=new double[n][n];
		double[][] r_=new double[n][n];
		double [][] q=new double[n][n];
		double[] s=new double[n];
		double[] t=new double[n];
		double[] w=new double[n];
		double[] xold=new double[n];
		double[] fvec=new double[n];
		double[] tmp1=new double[n];
		double[] tmp2=new double[n];
		final boolean [] check = new boolean[1];

		RealFunctionOfSeveralVariables fmin = Newton.normSqr( func, fvec );
		
		double f=fmin.eval(x);
		double test=0.0;
		
		for ( int i=0;i<n;i++)
			if ( Math.abs(fvec[i]) > test) test=Math.abs(fvec[i]);
			
		if (test<0.01*acc) {
			return;
		}
	
		double stpmax=maxStepSize*Math.max(Math.sqrt(VectorOperations.normSqr(x)),n);
	
		boolean restrt=true;
		
		for (int its=1;its<=maxNumberOfIterations;its++) {
			if (restrt) {

				func.eval(x,fvec,0,r);
			
				double det =Householder.decompose( r,q, tmp1, tmp2 );
				
				if (det==0)
					new RuntimeException("singular Jacobian in broydn");
	
			} else {
				
				VectorOperations.minus( x, xold, s);
				MatrixOperations.times( r, s, t );
				
				boolean skip=true;
				for (int i=0;i<n;i++) {
					double sum=0;
					for ( int j=0;j<n;j++) sum += q[i][j]*t[j];
					w[i]=fvec[i]-fvcold[i]-sum;
					if (Math.abs(w[i]) >= EPS*(Math.abs(fvec[i])+Math.abs(fvcold[i])))
						skip=false;
					else w[i]=0.0;
				}
				if (!skip) {
					MatrixOperations.times( w, q, t );
					
					VectorOperations.divide(s,VectorOperations.normSqr(s),s);
					
					// we have not qr update:qrupdt(r,qt,n,t,s);
					MatrixOperations.times(t,s, ts);
					MatrixOperations.plus( r, ts, r_);
					MatrixOperations.times( q, r_, r );
					double det = Householder.decompose( r, q, tmp1, tmp2 );
					
					if( det== 0 ) {
						 throw new RuntimeException("r singular");
					}
				}
			}
			
			MatrixOperations.times( fvec, q, tmp1 ); 
			MatrixOperations.times( tmp1, r, g );
			
			VectorOperations.assign(x,xold);
			VectorOperations.assign(fvec,fvcold);
			
			MatrixOperations.times(fvec, q ,p);
			VectorOperations.neg( p,p );
			
			RXB.solve(r,p);
			
			f = Newton.lnsrch(xold,f,g,p,x, precision, stpmax, check, fmin);
		
			test=0.0;
			for ( int i=0;i<n;i++)
				if (Math.abs(fvec[i]) > test) test=Math.abs(fvec[i]);
			if (test < acc) {
				return;
			}
			if (check[0]) {
				if (restrt) throw new HitLocalMinimumException();
				else {
					test=0.0;
					double den=Math.max(f,0.5*n);
					for (int i=0;i<n;i++) {
						double temp=Math.abs(g[i])*Math.max(Math.abs(x[i]),1.0)/den;
						if (temp > test) test=temp;
					}
					if (test < accOfMin ) {
						return;
					}
					else restrt=true;
				}
			} else {
				restrt=false;
				test=0.0;
				for (int i=0;i<n;i++) {
					double temp=(Math.abs(x[i]-xold[i]))/Math.max(Math.abs(x[i]),1.0);
					if (temp > test) test=temp;
				}
				if (test < precision) {
					return;
				}
			}
		}
		throw new RuntimeException("exeeded maximal number of iterations");
	}
	
}
