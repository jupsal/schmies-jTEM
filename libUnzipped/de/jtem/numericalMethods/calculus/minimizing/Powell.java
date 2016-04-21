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

package de.jtem.numericalMethods.calculus.minimizing;

import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;


/**
 * This class represents the routine to find a minimum of multidimensional function f with powell methods in mutidimensional.
 * In short: we start at initial point p and then move to the minimum across lines biuld from point p and vector's of given basis.
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public final class Powell
    implements java.io.Serializable {

  private static final long serialVersionUID = 1L;

  /**
   * maximum number of function evaluation in search method.
   */
  static int ITMAX = 200;

  /**
   * Get the value of ITMAX.
   * @return Value of ITMAX.
   */
  /** @deprecated */
  public static int getITMAX() {
    return ITMAX;
  }

  /**
   * Set the value of ITMAX.
   * @param v  Value to assign to ITMAX.
   */
  /** @deprecated*/
  public static void setITMAX(int v) {
    ITMAX = v;
  }


  /**
   * Returns square of a.
   * @param a number to square.
   * @return a*a.
   */
  static final double sqr(double a) {
    return a * a;
  }

  /**
   * Return the standard basis of dimension dim.
   * @param dim dimension of basis.
   * @return standard basis.
   */
  public static double[][] getStandardBasis(int dim) {

    double[][] basis = new double[dim][dim];

    for (int i = 0; i < dim; i++)
      basis[i][i] = 1;

    return basis;
  }

  /**
   * Search the minimum of function f  with powell's method in standard basis
   *  and return the value of the minimum.
   * @param p starting piont
   * @param ftol precision tolerance.
   * @param f given function.
   * @return the minimum of function f.
   */
  public final static double search(double[] p, double ftol,
                             RealFunctionOfSeveralVariables f) {

    return search(p, getStandardBasis(p.length), ftol, f, ITMAX, null);
  }

  /**
   * Search the minimum of function f  with powell's method in standard basis
   * and return the value of the minimum.
   * @param p starting piont
   * @param ftol precision tolerance.
   * @param itmax maximal number of iterations
   * @param f given function.
   * @return the minimum of function f.
   */
  public final static double search(double[] p, double ftol, int maxIteration,
                             RealFunctionOfSeveralVariables f) {

    return search(p, getStandardBasis(p.length), ftol, f, maxIteration, null );
  }

  /**
   * Search the minimum of function f  with powell's method in standard basis
   * and return the value of the minimum.
   * @param p starting piont
   * @param ftol precision tolerance.
   * @param itmax maximal number of iterations
   * @param f given function.
   * @param info object holding process information.
   * @return the minimum of function f.
   * @see Info
   */
  public final static double search(double[] p, double ftol, int maxIteration,
                             RealFunctionOfSeveralVariables f, Info info) {

    return search(p, getStandardBasis(p.length), ftol, f, maxIteration, info );
  }

  /**
   * Search the minimum of function f  with powell's method in initial basis xi
   *  and return the value of the minimum.
   * @param p starting point.
   * @param x initial basis.
   * @param ftol precision tolerance.
   * @param f given function.
   * @param info some gebug information.
   * @return the minimum of function f.
   * @see Info
   */
  public static double search(double[] p,
                       double[][] xi,
                       double ftol,
                       RealFunctionOfSeveralVariables f,
                       int itMax, Info info) {

    final double[] aTuple = new double[2];

    final int n = p.length;

    final double[] pt = new double[n];
    final double[] ptt = new double[n];
    final double[] xit = new double[n];

    BrentOnLine brentOnLine = new BrentOnLine(p, xit, f);

    double fptt, fret = f.eval(p);

    if (info!= null) {
      String s = new String(" f(p) = " + fret + " , p = ");

      for (int i = 0; i < n; i++)
        s += p[i] + " ";

      info.setMessage(s);
      info.setMaxIter(itMax);
    }

    for (int j = 0; j < n; j++)
      pt[j] = p[j];

    //System.out.println( "value=" + fret );
    
    for (int iter = 1; ; iter++) {
      double fp = fret;
      int ibig = 0;
      double del = 0.0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
          xit[j] = xi[j][i];
        fptt = fret;

        fret = brentOnLine.search(2.0e-8); //(ftol);

        if (Math.abs(fptt - fret) > del) {
          del = Math.abs(fptt - fret);
          ibig = i;
        }
      }

      if (2.0 * Math.abs(fp - fret) <= ftol * (Math.abs(fp) + Math.abs(fret))) {
        if (info != null) {
          String s = new String("iter = " + iter + ", fret = " + fret +
                                ", fp = " +
                                fp + ", p = ");

          for (int i = 0; i < n; i++)
            s += p[i] + " ";

          info.addMessage(s);

          info.setCurrentIter(iter);
          info.printDebug();
        }
        return fret;
      }

      if (info!=null && iter == itMax) {
        info.setCurrentIter(iter);
        info.setMessage("Too many iterations in routine POWELL");
        info.printDebug();
        return fret;
      }

      for (int j = 0; j < n; j++) {
        ptt[j] = 2.0 * p[j] - pt[j];
        xit[j] = p[j] - pt[j];
        pt[j] = p[j];
      }

      fptt = f.eval(ptt);

      if (fptt < fp) {
        double t = 2.0 * (fp - 2.0 * fret + fptt) * sqr(fp - fret - del)
            - del * sqr(fp - fptt);
        if (t < 0.0) {

          fret = brentOnLine.search(2.0e-8); //(ftol);

          for (int j = 0; j < n; j++)
            xi[j][ibig] = xit[j];
        }
      }
      //System.out.println( "value=" + fret );
    }
  }

}
