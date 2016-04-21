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

package de.jtem.numericalMethods.calculus.interpolation;

/**
 * Interpolation with rational function using.
 */
public class RationalInterpolation {

  static double TINY = 1e-25;

  /**
   * Interpolates function given by sorted table (tx,ty) at x.
   * This function creates two arrays and envokes then
   * {@link #interpolate(doublex,double[],double[],double[],double[],double[])
   * interpolate(x,tx,ty,dy,tmp1,tmp2)}.
   * @param x position of interpolation
   * @param tx table for x-values, must be sorted
   * @param ty table for y-values
   * @param dy first entry of array carries the derivative on output, may be null
   * @return interpolated value at x
   */
  public static double interpolate( double x, double[] tx, double[] ty ) {
    double [] tmp1 = new double[ tx.length ];
    double [] tmp2 = new double[ tx.length ];
    return interpolate( x, tx, ty, null, tmp1, tmp2 );
  }

  /**
   * Interpolates function given by sorted table (tx,ty) at x.
   * This function creates two arrays and envokes then
   * {@link #interpolate(doublex,double[],double[],double[],double[],double[])
   * interpolate(x,tx,ty,dy,tmp1,tmp2)}.
   * @param x position of interpolation
   * @param tx table for x-values, must be sorted
   * @param ty table for y-values
   * @param dy first entry of array carries the derivative on output, may be null
   * @return interpolated value at x
   */
  public static double interpolate( double x, double[] tx, double[] ty,
                                    double[] dy ) {
    double [] tmp1 = new double[ tx.length ];
    double [] tmp2 = new double[ tx.length ];
    return interpolate( x, tx, ty, dy, tmp1, tmp2 );
  }

  /**
   * Interpolates function given by sorted table (tx,ty) at x.
   * This function do not create any objects and should therefore
   * be the choice if you want to evaluate a table many time.
   * @param x position of interpolation
   * @param tx table for x-values, must be sorted
   * @param ty table for y-values
   * @param dy first entry of array carries the derivative on output, may be null
   * @param tmp1 array for temp. data, length coincide with table
   * @param tmp2 array for temp. data, length coincide with table
   * @return interpolated value at x
   */
  public static double interpolate( double x, double[] xa, double[] ya,
                                    double[] dy,
                                    double[] tmp1, double[] tmp2 ) {
    final double[] c = tmp1;
    final double[] d = tmp2;

    final int n = xa.length;

    if (ya.length != n || c.length != n || d.length != n )
      throw new IllegalArgumentException("length of arrays do not coincide");

    int ns = 0;

    double hh = Math.abs(x - xa[0]);

    for ( int i = 0; i < n; i++) {
      final double h = Math.abs(x - xa[i]);
      if (h == 0.0) {
        if (dy != null)
          dy[0] = 0;
        return ya[i];
      }
      else if (h < hh) {
        ns = i;
        hh = h;
      }
      c[i] = ya[i];
      d[i] = ya[i] + TINY;
    }

    double dy_ = 0, y = ya[ns--];

    for ( int m = 0; m+1 < n; m++) {
      for ( int i = 0; i+1 <= n - (m+1); i++) {
        final double w = c[i + 1] - d[i];
        final double h = xa[i+1 + m ] - x;
        final double t = (xa[i] - x) * d[i] / h;
        double dd = t - c[i+1];
        if (dd == 0.0)
          throw new RuntimeException();
        dd = w / dd;
        d[i] = c[i+1] * dd;
        c[i] = t * dd;
      }

      dy_ = 2 * (ns+1) < (n - (m+1)) ? c[ns+1] : d[ns--];
      y += dy_;

    }

    if (dy != null)
      dy[0] = dy_;
    return y;
  }


}




