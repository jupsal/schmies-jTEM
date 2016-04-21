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

package de.jtem.numericalMethods.algebra.loopGroup;

import de.jtem.numericalMethods.algebra.linear.MatrixOperations;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;
import de.jtem.numericalMethods.algebra.linear.decompose.Householder;

/** Computes the Iwasawa factorization of a <code>GL_n(C)</code>-valued
    loop.
@author Markus Schmies
 */
public class IwasawaFactorization {

  /** Compute Iwasawa factorization
      Input: (xRe, xIm), R
      Output: (fRe, fIm)
  */
  public static void factor(
    final double [][][] xRe, final double [][][] xIm,
    double [][][] fRe, double [][][] fIm,
    final int blockSize ) {

    final int n = xRe.length;
    final int N = xRe[0][0].length / 2;
    final int R = blockSize;
    final int rows = 2*N+1+R;
    final int m = rows*n;
    final int r = R*n;

    final double [][] uRe = new double[ m ][ r ];
    final double [][] uIm = new double[ m ][ r ];
    final double []   dRe = new double[ r ];
    final double []   dIm = new double[ r ];
    final double []   pRe = new double[ m ];
    final double []   pIm = new double[ m ];
    final double [][] zRe = new double[ n ][ m ];
    final double [][] zIm = new double[ n ][ m ];
    final double [][] gRe = new double[ m ][ n ];
    final double [][] gIm = new double[ m ][ n ];

    factor( xRe, xIm, fRe, fIm,
             uRe, uIm, dRe, dIm, pRe, pIm, zRe, zIm, gRe, gIm, blockSize );
  }

  /** Compute Iwasawa factorization.
      Input: (xRe, xIm)
      Output: (fRe, fIm)
      The rest are temporaries.
  */
  public static void factor(
    final double [][][] xRe, final double [][][] xIm,
    final double [][][] fRe, final double [][][] fIm,
    final double [][]   uRe, final double [][]   uIm,
    final double []     dRe, final double []     dIm,
    final double []     pRe, final double []     pIm,
    final double [][]   zRe, final double [][]   zIm,
    final double [][]   gRe, final double [][]   gIm,
    int blockSize ) {

    /* Each block has dimensions (2*N+1+R, R), where
       R=blockSize, N=degree
       Note that he number of rows and the number of columns grows
       with R, so it makes sense to take R greater than N.
    */

    if (blockSize <= 0) {
      throw new IllegalArgumentException( "blockSize must be positive");
    }

    final int n = xRe.length;
    final int N = (xRe[0][0].length - 1)/ 2;
    final int R = blockSize;

    // Create the large matrix u consisting of n*n blocks.
    // x -> u
    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {
        setBlock( xRe[i][j], uRe, i, j, N, R );
        setBlock( xIm[i][j], uIm, i, j, N, R );
      }
    }

    // Perform QR-factorization of the large matrix u.
    // The result is encoded in u and d.
    // p is a temporary.
    double [] t = Householder.decompose( uRe, uIm, dRe, dIm, pRe, pIm );


    // Check whether the determinant is 0.
    if( t[0] == 0. && t[1] == 0. )
      throw new RuntimeException( "not solvable" );

    for( int i = 0; i < n; i++ ) {

      // set block x->z
      for( int j = 0; j < n; j++ ) {
        setBlock( xRe[j][i], zRe[i], j, N, R );
        setBlock( xIm[j][i], zIm[i], j, N, R );
      }

      // copy z to p
      VectorOperations.assign( zRe[i], pRe );
      VectorOperations.assign( zIm[i], pIm );

      // p = Ustar *  p
      Householder.qTimes( uRe, uIm, pRe, pIm, true );

      // project
      for( int k = R*n; k < pRe.length; k++ )
        pRe[k] = pIm[k] = 0;

      // p = U * p
      Householder.qTimes( uRe, uIm, pRe, pIm, false );

      // y = y - p
      VectorOperations.minus( zRe[i], pRe, zRe[i] );
      VectorOperations.minus( zIm[i], pIm, zIm[i] );

    }

    // Normalization

    MatrixOperations.transpose( zRe, gRe );
    MatrixOperations.transpose( zIm, gIm );

    double [] T = Householder.decompose( gRe, gIm, dRe, dIm, pRe, pIm );
    // use p temp.

    // Check whether the determinant is 0.
    if( T[0] == 0. && T[1] == 0. )
      throw new RuntimeException( "not solvable" );

    for( int i=0; i<n; i++ ) {

      assignZero( zRe[i] );
      assignZero( zIm[i] );

      zRe[i][i] = 1;

      Householder.qTimes( gRe, gIm, zRe[i], zIm[i], false );
    }

    // normalize QR decomposition such that diagonals of R are
    // positive reals.
    for( int i=0; i<n; i++ ) {

      final double absOfD = Math.sqrt( dRe[i]*dRe[i] + dIm[i]*dIm[i] );

      final double nRe = dRe[i] / absOfD;
      final double nIm = dIm[i] / absOfD;

      VectorOperations.times( zRe[i], zIm[i], nRe, nIm, zRe[i], zIm[i] );
    }

    // Extract the results from y into f
    for( int i=0; i<n; i++ ) {
      for( int j=0; j<n; j++ ) {
        getBlock( fRe[j][i], zRe[i], j, N, R );
        getBlock( fIm[j][i], zIm[i], j, N, R );
      }
    }

  }

  /** Matrix Laurent series multiplication b = star(f) * x
      Assumes degrees
      f: (-N, N),
      x: (-N, N),
      b: (0, N)
      Used to compute the positive part of the factorization.
   */
  public static void positivePart(
    double [][][] xRe, double [][][] xIm,
    double [][][] fRe, double [][][] fIm,
    double [][][] bRe, double [][][] bIm ) {

    int N = xRe[0][0].length;
    double [] tmpRe = new double [N+1];
    double [] tmpIm = new double [N+1];

    positivePart( xRe, xIm, fRe, fIm, bRe, bIm, tmpRe, tmpIm );
  }

  public static void positivePart(
    double [][][] xRe, double [][][] xIm, 
    double [][][] fRe, double [][][] fIm,
    double [][][] bRe, double [][][] bIm,
    double [] tmpRe, double [] tmpIm ) {

    times( fRe, fIm, xRe, xIm, bRe, bIm, tmpRe, tmpIm );
  }

  private static void times(
    double [][][] fRe, double [][][] fIm,
    double [][][] xRe, double [][][] xIm,
    double [][][] bRe, double [][][] bIm,
    double [] tmpRe, double [] tmpIm ) {

    int n = fRe.length;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        polynomialTimes( fRe[0][i], fIm[0][i],
                         xRe[0][j], xIm[0][j],
                         bRe[i][j], bIm[i][j] );

        for( int k=1; k<n; k++ ) {
          polynomialTimes( fRe[k][i], fIm[k][i],
                           xRe[k][j], xIm[k][j],
                           tmpRe,     tmpIm );
          
          VectorOperations.plus( bRe[i][j], tmpRe, bRe[i][j] );
          VectorOperations.plus( bIm[i][j], tmpIm, bIm[i][j] );
        }
      }
    }
  }

  /** Laurent series multiplication b = star(f) * x
      Assumes degrees
      f: (-N, N),
      x: (-N, N),
      b: (0, N)
      Used to compute the positive part of the factorization.
  */
  private static void polynomialTimes(
    double [] fRe, double [] fIm,
    double [] xRe, double [] xIm,
    double [] bRe, double [] bIm ) {

    if (fRe == bRe || xRe == bRe) {
      throw new IllegalArgumentException();
    }

    int N = (fRe.length - 1)/2;

    for ( int i = 0; i <= N; i++ ) {
      double br = 0.0;
      double bi = 0.0;
      for ( int j = 0; j <= 2*N-i; j++) {
        br += fRe[j] * xRe[i+j] + fIm[j] * xIm[i+j];
        bi += fRe[j] * xIm[i+j] - fIm[j] * xRe[i+j];
      }
      bRe[i] = br;
      bIm[i] = bi;
    }
  }


  /** Put coeff into m in block form.
    */

  private static void
  setBlock(
    final double [] coeff,
    final double [][] m,
    final int row,
    final int col,
    final int N,
    final int R ) {

    final int rows = 2*N+1+R;
    final int rowOffset = row * ( 2*N+1 + R );
    final int colOffset = col * R;

    for (int j = 0; j < R; j++) {

      int p0 = j+1;
      int p1 = Math.min(rows, 2*N+2+j);

      // initial zeros
      for (int i = 0; i < p0; i++) {
        m[i + rowOffset][j + colOffset] = 0.0;
      }

      // elements
      for (int i = p0; i < p1; i++) {
        m[i + rowOffset][j + colOffset] = coeff[i-j-1];
      }

      // final zeros
      for (int i = p1; i < rows; i++) {
        m[i + rowOffset][j + colOffset] = 0.0;
      }
    }
  }

  /** Put coeff into x in block form.
    */
  private static void
  setBlock(
    final double [] coeff,
    final double [] x,
    final int row,
    final int N,
    final int R ) {

    // elements
    System.arraycopy( coeff, 0, x, row*(2*N+1+R), 2*N+1 );

    // zeros
    for( int i = row*(2*N+1+R) + 2*N + 1, last = i + R; i<last; i++ )
      x[i] = 0;
  }

  /** Extract data from the block x into coeff.
    */
  private static void
  getBlock(
    final double [] coeff,
    final double [] x,
    final int row,
    final int N,
    final int R ) {
    System.arraycopy( x, row*(2*N+1+R), coeff, 0, 2*N+1 );
  }

  /** Set the entries of x to 0.
   */
  private static void
  assignZero( final double [] x ) {
    for( int i=0; i<x.length; i++ )
      x[i]=0;
  }
}
