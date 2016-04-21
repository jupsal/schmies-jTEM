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

/**
 * Newton-Raphson method's for finding roots of real and complex functions.
 *
 * This is defitnitly most famos root-finding methods for one dimensional problems.
 * It requires the function and its derivative. It is well
 * known that Newton?s method has only poor global convergence properties, but
 * <em>near</em> a root it converges quadratically. Therefor a common strategie
 * is to use a globally more stable method first and than, when you are <em>near</em>
 * to a root, to perform a few steps of Newton?s methods, which approxiamtely doubles 
 * the number of significant digits with each step.
 * <p>
 * Functions are defined using the inner interfaces 
 * {@link RealFunctionWithDerivative} and
 * {@link ComplexFunctionWithDerivative}.
 * <p>
 * Below the code which is needed to find the root of sin(x).
 * <pre>
 *     NewtonRaphson.RealFunctionWithDerivative rf = new NewtonRaphson.RealFunctionWithDerivative() {
 *	         public void eval(double x, double[] f, int offsetF, double[] df, int offsetDF) {
 *					f[offsetF] = Math.sin(x);
 *					df[offsetDF] = Math.cos(x);
 *			}
 *		};
 *
 *		double[] rootValues = new double[3];
 *
 *		NewtonRaphson.search(rf, 0.6, rootValues ); // start our serch in 0.6
 * </pre>
 *
 * @see de.jtem.numericalMethods.calculus.rootFinding.ExampleNewtonRaphson
 *
 * @author Schmies 
 */
public class NewtonRaphson {

	/** Default accuracy: <code>1e-12</code>. */
	public static double ACC = 1e-12;

	/** Default maximum number of steps performed in Newton?s method: <code>50</code>. */
	public static int MAX_NUM_OF_STEPS = 50;

	/** 
	 * Switches debug output on and off; the default is <code>false</code>.
	 * If this flag is <code>true</code>, the methods will list all
	 * the evaluated values.
	 */
	public static boolean DEBUG = false;

	/**
	 * Local interface defining a real function with derivative.
	 */
	public interface RealFunctionWithDerivative {

		/**
		 * Evaluates <code>f</code> and <code>f'</code> in <code>x</code>
		 * and returns the results in the provided arrays f and df, which may be null.
		 * @param x argument
		 * @param f value of f(z) on output, ( f[offsetF] = f(x)) )
		 * @param offsetF
		 * @param df vale of f'(z) on output, ( df[offsetDF] = re(df(x)) )
		 * @param offsetDF
		 */
		public void eval(double x, double[] f, int offsetF, double[] df, int offsetDF);
	}

	/**
	 * Local interface defining a complex function with derivative.
	 */
	public interface ComplexFunctionWithDerivative {

		/**
		 * Evaluates <code>f</code> and <code>f'</code> in <code>z=x+iy</code>
		 * and returns the results in the provided arrays f and df, which may be null.
		 * @param x real part of argument
		 * @param y imaginary part of argument
		 * @param f value of f(z) on output, ( f[offsetF] = re(f(z)), f[offsetF+1] = im(f(z)) )
		 * @param offsetF
		 * @param df vale of f'(z) on output, ( df[offsetDF] = re(df(z)), df[offsetDF+1] = im(df(z)) )
		 * @param offsetDF
		 */
		public void eval(double x, double y, double[] f, int offsetF, double[] df, int offsetDF);
	}

	private NewtonRaphson() {
	}

	/**
	 * Performs Newton?s method starting at <code>x</code>.
	 * The method returns the root <code>x<sub>r</sub></code> 
	 * of a real function <code>f</code>,
	 * in the first entry of <code>rootValues</code>. The second and third
	 * entry of <code>rootValues</code> are assigned with
	 * <code>f(x<sub>r</sub>)</code> and <code>g?(x<sub>r</sub>)</code>.
	 * The algorithm terminates if the approximation of the
	 * root value changes relatively less than {@link #ACC} or the
	 * function at this value drops below {@link #ACC}.
	 * @param f real function with derivative information
	 * @param x start position
	 * @param rootValues array of length <code>3</code> for the result
	 * @exception RuntimeException when {@link #MAX_NUM_OF_STEPS} was exeeded
	 */
	public static void search(final RealFunctionWithDerivative f, double x, final double[] rootValues) {

		search(f, x, rootValues, ACC, ACC, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, MAX_NUM_OF_STEPS, true);
	}

	/**
	 * Performs Newton?s method starting at <code>x</code>.
	 * The method returns the root <code>x<sub>r</sub></code> 
	 * of a real function <code>f</code>,
	 * in the first entry of <code>rootValues</code>. The second and third
	 * entry of <code>rootValues</code> are assigned with
	 * <code>f(x<sub>r</sub>)</code> and <code>g?(x<sub>r</sub>)</code>.
	 * The algorithm terminates if the approximation of the
	 * root value changes relatively less than <code>acc</code> or the
	 * function at this value drops below <code>acc</code>.
	 * @param f real function with derivative information
	 * @param x start position
	 * @param rootValues array of length <code>3</code> for the result
	 * @exception RuntimeException when {@link #MAX_NUM_OF_STEPS} was exeeded
	 */
	public static void search(final RealFunctionWithDerivative f, double x, final double[] rootValues, double acc) {

		search(f, x, rootValues, acc, acc, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, MAX_NUM_OF_STEPS, true);
	}

	/**
	 * Performs Newton?s method starting at <code>x</code> and
	 * restricting itself to the interval <code>[xMin,xMax]</code>.
	 * The method returns the root <code>x<sub>r</sub></code> 
	 * of a real function <code>f</code>,
	 * in the first entry of <code>rootValues</code>. The second and third
	 * entry of <code>rootValues</code> are assigned with
	 * <code>f(x<sub>r</sub>)</code> and <code>g?(x<sub>r</sub>)</code>.
	 * The algorithm terminates if the approximation of the
	 * root value changes relatively less than <code>accX</code> or the
	 * function at this value drops below <code>acc</code>.
	 * @param f real function with derivative information
	 * @param x start position in the interval <code>[xMin,xMax]</code>
	 * @param rootValues array of length <code>3</code> for the result
	 * @param acc accuracy of the zero (absolute)
	 * @param accX accuracy of root (relative)
	 * @param xMin lower bound of search interval
	 * @param xMax upper bound of search interval
	 * @param panicIfMaxNumOfStepsAreExeeded controls whether an exception is thrown if
	 * maximum number of steps are exeeded
	 * @exception RuntimeException when the algorithm leaves the search interval <code>[xMin,xMax]</code>
	 * @exception RuntimeException when <code>maxNumOfSteps</code> was exeeded and 
	 * <code>panicIfMaxNumOfStepsAreExeeded</code> is <code>true</code>
	 */
	public static void search(
		final RealFunctionWithDerivative f,
		double x,
		final double[] rootValues,
		final double acc,
		final double accX,
		final double xMin,
		final double xMax,
		final int maxNumOfSteps,
		final boolean panicIfMaxNumOfStepsAreExeeded) {

		if (DEBUG)
			System.out.println("NewtonRaphson DEBUG: " + x);

		for (int i = 0; i < maxNumOfSteps; i++) {

			f.eval(x, rootValues, 1, rootValues, 2);

			if (DEBUG) {
				System.out.println("NewtonRaphson :  x  = " + x);
				System.out.println("               f(x) = " + rootValues[1]);
				System.out.println("              f'(x) = " + rootValues[2]);
			}

			double dx = rootValues[1] / rootValues[2];

			rootValues[0] = x -= dx;

			if (x < xMin || x > xMax)
				throw new RuntimeException(" left initial brackets ");

			if (Math.abs(dx / x) < accX || Math.abs(rootValues[1]) < acc)
				return;
			if (Math.abs(dx) < acc)
				return;
		}

		if (panicIfMaxNumOfStepsAreExeeded)
			throw new RuntimeException(" exeed max number of steps !!! ");
	}

	/**
	 * Performs Newton?s method starting at <code>z=x+iy</code>.
	 * The method returns the root <code>z<sub>r</sub></code> 
	 * of a complex function <code>f</code>,
	 * in the first two entries of <code>rootValues</code>. The third to the sixth
	 * entry of <code>rootValues</code> are assigned with
	 * <code>re(g(z<sub>r</sub>), im(g(z<sub>r</sub>), re(g?(z<sub>r</sub>))</code> 
	 * and <code>g?(im(z<sub>r</sub>))</code>.
	 * The algorithm terminates if the approximation of the
	 * root value changes relatively less than {@link #ACC} or the
	 * function at this value drops below {@link #ACC}.
	 * @param f complex function with derivative information
	 * @param x real part of start position
	 * @param y imaginary part of start position
	 * @param rootValues array of length <code>6</code> for the result
	 * @exception RuntimeException when {@link #MAX_NUM_OF_STEPS} was exeeded
	 */
	public static void search(final ComplexFunctionWithDerivative f, double x, double y, final double[] rootValues) {

		search(
			f,
			x,
			y,
			rootValues,
			ACC,
			ACC,
			ACC,
			Double.NEGATIVE_INFINITY,
			Double.POSITIVE_INFINITY,
			Double.NEGATIVE_INFINITY,
			Double.POSITIVE_INFINITY,
			MAX_NUM_OF_STEPS,
			true);

	}

	/**
	 * Performs Newton?s method starting at <code>z=x+iy</code>.
	 * The method returns the root <code>z<sub>r</sub></code> 
	 * of a complex function <code>f</code>, 
	 * in the first two entries of <code>rootValues</code>. The third to the sixth
	 * entry of <code>rootValues</code> are assigned with
	 * <code>re(g(z<sub>r</sub>), im(g(z<sub>r</sub>), re(g?(z<sub>r</sub>))</code> 
	 * and <code>g?(im(z<sub>r</sub>))</code>.
	 * The algorithm terminates if the approximation of the
	 * root value changes relatively less than <code>acc</code> or the
	 * function at this value drops below <code>acc</code>.
	 * @param f complex function with derivative information
	 * @param x real part of start position
	 * @param y imaginary part of start position
	 * @param rootValues array of length <code>6</code> for the result
	 * @exception RuntimeException when {@link #MAX_NUM_OF_STEPS} was exeeded
	 */
	public static void search(final ComplexFunctionWithDerivative f, double x, double y, final double[] rootValues, double acc) {

		search(
			f,
			x,
			y,
			rootValues,
			acc,
			acc,
			acc,
			Double.NEGATIVE_INFINITY,
			Double.POSITIVE_INFINITY,
			Double.NEGATIVE_INFINITY,
			Double.POSITIVE_INFINITY,
			MAX_NUM_OF_STEPS,
			true);

	}

	/**
	  * Performs Newton?s method starting at <code>z=x+iy</code> and
	  * restricting itself to the rectangle <code>[xMin,xMax]x[yMin,yMax]</code>.
	  * The method returns the root <code>z<sub>r</sub></code> 
	  * of a complex function <code>f</code>,
	  * in the first two entries of <code>rootValues</code>. The thrid to the forth
	  * entry of <code>rootValues</code> are assigned with
	  * <code>re(g(z<sub>r</sub>), im(g(z<sub>r</sub>), re(g?(z<sub>r</sub>))</code> 
	  * and <code>g?(im(z<sub>r</sub>))</code>.
	  * The algorithm terminates if the approximation of the
	  * root value changes relatively less than <code>accX</code> in the 
	  * real component and <code>accY</code> in the imaginary or the
	  * function at this value drops below <code>acc</code>.
	  * @param f complex function with derivative information
	  * @param x start position in the interval <code>[xMin,xMax]</code>
	  * @param rootValues array of length <code>6</code> for the result
	  * @param acc accuracy of the zero (absolute)
	  * @param accX accuracy of real part of root (relative)
	  * @param accY accuracy of imaginary part of root (relative)
	  * @param xMin lower bound of real part of searach rectangle
	  * @param xMax upper bound of real part of searach rectangle
	  * @param yMin lower bound of imaginary part of searach rectangle
	  * @param yMax upper bound of imaginary part of searach rectangle
	  * @param panicIfMaxNumOfStepsAreExeeded controls whether an exception is thrown if
	  * maximum number of steps are exeeded
	  * @exception RuntimeException when the algorithm leaves the search rectangle
	  * <code>[xMin,xMax]x[yMin,yMax]</code>
	  * @exception RuntimeException when <code>maxNumOfSteps</code> was exeeded and 
	  * <code>panicIfMaxNumOfStepsAreExeeded</code> is <code>true</code>
	  */
	public static void search(
		final ComplexFunctionWithDerivative f,
		double x,
		double y,
		final double[] rootValues,
		final double acc,
		final double accX,
		final double accY,
		final double xMin,
		final double xMax,
		final double yMin,
		final double yMax,
		final int maxNumOfSteps,
		final boolean panicIfMaxNumOfStepsAreExeeded) {

		if (DEBUG)
			System.out.println("NewtonRaphson DEBUG: " + x + " + i " + y);

		for (int i = 0; i < maxNumOfSteps; i++) {

			rootValues[0] = x;
			rootValues[1] = y;

			f.eval(x, y, rootValues, 2, rootValues, 4);

			double denom = rootValues[4] * rootValues[4] + rootValues[5] * rootValues[5];

			double dx = (rootValues[2] * rootValues[4] + rootValues[3] * rootValues[5]) / denom;
			double dy = (rootValues[3] * rootValues[4] - rootValues[2] * rootValues[5]) / denom;

			if (DEBUG) {
				System.out.println("NewtonRaphson :  x  = " + x + " + i " + y);
				System.out.println("               f(x) = " + rootValues[2] + " + i " + rootValues[3]);
				System.out.println("              f'(x) = " + rootValues[4] + " + i " + rootValues[5]);

				System.out.println("              dx = " + dx + " + i " + dy);
			}

			x -= dx;
			y -= dy;

			if (x < xMin || x > xMax || y < yMin || y > yMax)
				throw new RuntimeException(" left initial brackets ");

			double absZ = Math.abs(x) + Math.abs(y);

			if (Math.abs(dx / absZ) < accX && Math.abs(dy / absZ) < accY || Math.abs(rootValues[2]) < acc && Math.abs(rootValues[3]) < acc)
				return;
		}

		if (panicIfMaxNumOfStepsAreExeeded)
			throw new RuntimeException(" exeed max number of steps !!! ");
	}
}
