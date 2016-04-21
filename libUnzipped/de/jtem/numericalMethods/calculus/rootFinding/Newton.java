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
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.function.RealVectorValuedFunctionOfSeveralVariablesWithJacobien;

/**
 * @author schmies
 */
public class Newton {

	final static double DEFAULT_MAX_STEPSIZE = 1;
	Newton() {
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
			double acc, double accOfMIn, double precision, double maxStepSize  )
	{
			
		final int n = func.getDimensionOfTargetSpace();
		
		final double [][] fjac= new double[n][n];
		final double [][] q= new double[n][n];
		final double [] g=new double[n];
        final double [] p=new double[n];
		final double [] xold=new double[n];
		final double [] fvec=new double[n];
		double[] tmp1=new double[n];
		double[] tmp2=new double[n];
		final boolean [] check = new boolean[1];

		RealFunctionOfSeveralVariables fmin = normSqr( func, fvec );
		
		double f=fmin.eval(x);
		double test=0.0;
		
		for ( int i=0;i<n;i++)
			if ( Math.abs(fvec[i]) > test) test=Math.abs(fvec[i]);
			
		if (test<0.01*acc) {
			return;
		}
			
		double stpmax=maxStepSize*Math.max(Math.sqrt(VectorOperations.normSqr(x)),n);
		
		for (int its=1;its<=maxNumberOfIterations;its++) {
			func.eval(x,fvec,0,fjac);
	
			MatrixOperations.times( fvec, fjac, g );
			VectorOperations.assign(x,xold);

			double det =Householder.decompose( fjac,q, tmp1, tmp2 ); // jfac is upper right tri after
			
			if (det==0)
				new RuntimeException("singular Jacobian");
			
			MatrixOperations.times(fvec, q ,p);
			VectorOperations.neg( p,p );
			
			RXB.solve(fjac,p);
			
			f = lnsrch(xold,f,g,p,x, precision, stpmax, check, fmin);
		
			test=0.0;
			for (int i=0;i<n;i++)
				if (Math.abs(fvec[i]) > test) test=Math.abs(fvec[i]);
			if (test < acc) {
				return;
			}
			if (check[0]) {
				test=0.0;
				double den=Math.max(f,0.5*n);
				for (int i=0;i<n;i++) {
					double temp=Math.abs(g[i])*Math.max(Math.abs(x[i]),1.0)/den;
					if (temp > test) test=temp;
				}
				if( test < accOfMIn )
					throw new HitLocalMinimumException();
				return;
			}
			test=0.0;
			for (int i=0;i<n;i++) {
				double temp=(Math.abs(x[i]-xold[i]))/Math.max(Math.abs(x[i]),1.0);
				if (temp > test) test=temp;
			}
			if (test < precision) {
				return;
			}
		}
		throw new RuntimeException("exeed maximal number of iterations");
	}
	
	static RealFunctionOfSeveralVariables normSqr(
			final RealVectorValuedFunctionOfSeveralVariables F, final double [] v ) {
		
		return  new RealFunctionOfSeveralVariables () {
			
			final double [] values = v;
			
			public double eval(double[] x) {
				F.eval(x,values,0);
				return VectorOperations.normSqr(values) / 2;
			}

			public int getNumberOfVariables() {
				return F.getNumberOfVariables();
			}
			
		};
	}
	
	final static double ALF = 1.0e-4;

	public static double lnsrch(double xold[], double fold, double g[], double p[], double x[],
			double precision, double stpmax, boolean [] check, RealFunctionOfSeveralVariables func )
	{
		final int n = func.getNumberOfVariables();
		
		check[0]=false;
		double sum = 0;
		for ( int i=0;i<n;i++) 
			sum += p[i]*p[i];
		sum=Math.sqrt(sum);
		if (sum > stpmax)
			for (int i=0;i<n;i++) 
				p[i] *= stpmax/sum;
		double slope=0;
		for (int i=0;i<n;i++)
			slope += g[i]*p[i];
		double test=0.0;
		for ( int i=0;i<n;i++) {
			double temp=Math.abs(p[i])/Math.max(Math.abs(xold[i]),1.0);
			if (temp > test) test=temp;
		}
		double alamin=precision/test;
		double alam=1.0;
		double alam2=0, tmplam, f2=0, fold2=0;
		for ( int iter=0;; iter++) {
			for (int i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
			double f=func.eval(x);
			if (alam < alamin) {
				for (int i=0;i<n;i++) 
					x[i]=xold[i];
				check[0]=true;
				return f;
			} else if (f <= fold+ALF*alam*slope) return f;
			else {
				
				if (alam == 1.0) {
					 tmplam = -slope/(2.0*(f-fold-slope));
				} else {
					double rhs1 = f-fold-alam*slope;
					double rhs2=f2-fold2-alam2*slope;
					double a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
					double b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
					if (a == 0.0) tmplam = -slope/(2.0*b);
					else {
						double disc=b*b-3.0*a*slope;
						if (disc<0.0) 
							throw new RuntimeException("Roundoff problem in lnsrch.");
						else tmplam=(-b+Math.sqrt(disc))/(3.0*a);
					}
					if (tmplam>0.5*alam)
						tmplam=0.5*alam;
				}
			}
			alam2=alam;
			f2 = f;
			fold2=fold;
			alam=Math.max(tmplam,0.1*alam);
		}
	}
}
