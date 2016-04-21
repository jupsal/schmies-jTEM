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

import de.jtem.numericalMethods.algebra.linear.Inversion;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;
import de.jtem.numericalMethods.algebra.linear.decompose.Householder;
import de.jtem.numericalMethods.algebra.linear.solve.RXB;
import de.jtem.numericalMethods.algebra.polynomial.ComplexPolynomial;

/** Computes the Birkhoff factorization of a <code>GL_n(C)</code>-valued
    loop X = X1 X2.

    Depending on flags set, X1 may be a "negative" loop,
    extending to the exterior disk on CP^1 bounded by the loop containing
    infinity, and X2 a "positive" loop, extending to the interior disk
    containng 0, or vice-versa.

    Depending on flags set, either X1 or X2 may be chosen to be normalized
    to I at 0 or infinity.

    This algorithm does not apply to loops in the "big cell" (e.g. a
    non-constant monomial) and do not check for this condition.

    The class consists of two versions of the factorization
    method factor(), one to which
    temporary storage space is passed (for efficiency) and one which
    allocates its own storage (for convenience).

    @author Markus Schmies and Nick Schmitt
    @version 3/4/2003
*/

public class BirkhoffFactorization {

  public static final int NEG_POS = 0;
  public static final int POS_NEG = 1;
  public static final int FIRST = 0;
  public static final int SECOND = 1;

  /** Computes one of the factors of a Birkhoff factorization of a
      <code>GL_n(C)</code>-valued loop <code>X</code>.

      @param xRe The real part of <code>X</code> (input).
      @param xIm The imaginary part of <code>X</code> (input).
      @param yRe The real part of <code>XNeg</code> (output).
      @param yIm The imaginary part of <code>XNeg</code> (output).
      @param direction When set to NEG_POS, computes the Birkhof factorization
      X = Xneg Xpos, where Xneg is a negative loop and Xpos is a positive loop.
      When set to POS_NEG, computes the Birkhof factorization
      X = Xpos Xneg.
      @param factor When set to FIRST (resp. SECOND), computes in
      (yRe, yIm) the first (resp. second) factor of the Birkhoff factorization.
      @param normalize When set to FIRST, the first factor is normalized
      to be I at 0 (in the case direction = NEG_POS) or at infinity
      (in the case direction = POS_NEG).
   */
  public static void
  factor(
    final double [][][] xRe, final double [][][] xIm,
    final double [][][] yRe, final double [][][] yIm,
    int direction, int factor, int normalize ) {

    final int n = xRe.length;
    final int N = (xRe[0][0].length - 1) / 2;
    int blockSize = N+1;
    int m = blockSize*n;

    final double [][][]   tRe    = new double [ n ][ n ][ N+1 ];
    final double [][][]   tIm    = new double [ n ][ n ][ N+1 ];
    final double [][]   qRe    = new double [ m ][ m ];
    final double [][]   qIm    = new double [ m ][ m ];
    final double []     dRe    = new double [ m ];
    final double []     dIm    = new double [ m ];
    final double []     kRe    = new double [ m ];
    final double []     kIm    = new double [ m ];
    final double [][]   yyyRe    = new double [ n ][ m ];
    final double [][]   yyyIm    = new double [ n ][ m ];
    final double [][]   q0Re   = new double [ n ][ n ];
    final double [][]   q0Im   = new double [ n ][ n ];
    final double [][]   t0Re   = new double [ n ][ n ];
    final double [][]   t0Im   = new double [ n ][ n ];
    final double [][][] uRe    = new double [ n ][ n ][ N+1 ];
    final double [][][] uIm    = new double [ n ][ n ][ N+1 ];
    final double []     tmpmRe = new double [ 2*N+1 ];
    final double []     tmpmIm = new double [ 2*N+1 ];

    factor( xRe, xIm, yRe, yIm, tRe, tIm,
             qRe, qIm, dRe, dIm, kRe, kIm, yyyRe, yyyIm, q0Re, q0Im,
            t0Re, t0Im, uRe, uIm,
             tmpmRe, tmpmIm,
             direction, factor, normalize );
  }


/**
      @param qRe Temporary storage.
      @param qIm Temporary storage.
      @param dRe Temporary storage.
      @param dIm Temporary storage.
      @param kRe Temporary storage.
      @param kIm Temporary storage.
      @param yyyRe Temporary storage.
      @param yyyIm Temporary storage.
      @param q0Re Temporary storage.
      @param q0Im Temporary storage.
      @param tmpmRe Temporary storage.
      @param tmpmIm Temporary storage.
   */
  public static void
  factor( final double [][][] xRe,       final double [][][] xIm,
          final double [][][] yRe,    final double [][][] yIm,
          final double [][][] tRe,    final double [][][] tIm,
          final double [][]   qRe,       final double [][]   qIm,
          final double []     dRe,       final double []     dIm,
          final double []     kRe,       final double []     kIm,
          final double [][]   yyyRe,       final double [][]   yyyIm,
          final double [][]   q0Re,      final double [][]   q0Im,
          final double [][]   t0Re,      final double [][]   t0Im,
          final double [][][]   uRe,      final double [][][]    uIm,
          final double []     tmpmRe,    final double []     tmpmIm,
          int direction, int factor, int normalize ) {

    final int n = xRe.      length;
    final int N = (xRe[0][0].length - 1) / 2;
    int blockSize = N+1;

    if (blockSize <= N) {
      throw new IllegalArgumentException(
        "blockSize must be at least degree+1" );
    }

    final boolean negFlag = isNegative(direction, factor);
    final boolean transposeFlag = factor == SECOND;

    // put x into q in block form.
    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {
        if (!transposeFlag) {
          setBlock( xRe[i][j], qRe, i, j, N, blockSize );
          setBlock( xIm[i][j], qIm, i, j, N, blockSize );
        }
        else {
          setBlock( xRe[i][j], qRe, j, i, N, blockSize );
          setBlock( xIm[i][j], qIm, j, i, N, blockSize );
        }
      }
    }

    // Perform QR factorization of x.
    // The factorization is encoded in q and d.
    double [] t = Householder.decompose( qRe, qIm, dRe, dIm, kRe, kIm );
    // use x temp.

    // Check whether the determinant is 0.
    if( t[0] == 0. && t[1] == 0. )
      throw new RuntimeException( "not solvable" );

    // Solve the system A x = b via QR-factorization
    // Factor A = QR. Then Q R x = b, so R x = Q^* b
    // The following loop is computing Q^* b.
    for( int i = 0; i < n; i++ ) { // i is the column

      // the right-hand side of the equation
      assignZero( yyyRe[i] );
      assignZero( yyyIm[i] );

      if (negFlag) {
        yyyRe[i][ i*blockSize ] = 1.0;
      }
      else {
        yyyRe[i][ i*blockSize + blockSize-1 ] = 1.0;
      }

      Householder.qTimes( qRe, qIm, yyyRe[i], yyyIm[i], true);

    }

    // set diagonal for upper right matrix
    for( int i = 0; i<dRe.length; i++ ) { // i is a column
      qRe[i][i] = dRe[i];
      qIm[i][i] = dIm[i];
    }

    // Solve the linear system.
    // The result is in (tRe, tIm).
    for( int i = 0; i < n; i++ ) { // i is a column

      // Solve q k = y
      // Input: q, y
      // Output: k
      RXB.solve( qRe, qIm, kRe, kIm, yyyRe[i], yyyIm[i] );


      // put k into t
      int offset = negFlag ? 0 : blockSize-(N+1);
      if (!transposeFlag) {
        for( int j = 0; j < n; j++ ) {
          getBlock( tRe[j][i], kRe, j*blockSize+offset, N );
          getBlock( tIm[j][i], kIm, j*blockSize+offset, N );
        }
      }
      else {
        for( int j = 0; j < n; j++ ) {
          getBlock( tRe[i][j], kRe, j*blockSize+offset, N );
          getBlock( tIm[i][j], kIm, j*blockSize+offset, N );
        }
      }
    }

    double [][][] inRe = tRe;
    double [][][] inIm = tIm;
    double [][][] outRe = yRe;
    double [][][] outIm = yIm;

    if ( factor == FIRST && normalize == SECOND
         ||
         factor == SECOND && normalize == FIRST ) {

      // normalize t. Put result in u.
      normalize(
        tRe, tIm, uRe, uIm,
        t0Re, t0Im, q0Re, q0Im, dRe, dIm, yyyRe[0], yyyIm[0],
        factor==SECOND, !negFlag );

      inRe = uRe;
      inIm = uIm;
      outRe = yRe;
      outIm = yIm;
    }

    // Multiply
    // (qRe, qIm) are temporaries
    if (negFlag) {
      if (!transposeFlag) {
        blockMultiplyPosNeg(
          xRe, xIm, inRe, inIm, outRe, outIm, tmpmRe, tmpmIm );
      }
      else {
        blockMultiplyPosNegReverse(
          xRe, xIm, inRe, inIm, outRe, outIm, tmpmRe, tmpmIm );
      }
    }
    else {
      if (!transposeFlag) {
        blockMultiplyNegPos(
          xRe, xIm, inRe, inIm, outRe, outIm, tmpmRe, tmpmIm );
      }
      else {
        blockMultiplyNegPosReverse(
          xRe, xIm, inRe, inIm, outRe, outIm, tmpmRe, tmpmIm );
      }
    }
  }


  /**
 * @param direction
 * @param factor
 * @return
 */
public static boolean isNegative(int direction, int factor) {
	return (direction == NEG_POS && factor == FIRST
                             ||
                             direction == POS_NEG && factor == SECOND);
}


/** Special block Laurent series multiplication c = a * b
      Assumes lower and upper degrees
      a: (-n, n)
      b: (0, n)
      c: (-n, 0)
   */
  private static void
  blockMultiplyPosNeg(
    final double [][][] aRe,   final double [][][] aIm,
    final double [][][] bRe,   final double [][][] bIm,
    final double [][][] cRe,   final double [][][] cIm,
    final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe || bRe == cRe) {
      throw new IllegalArgumentException();
    }

    final int n = aRe.length;

    final int degOfA = aRe[0][0].length - 1;
    final int degOfB = bRe[0][0].length - 1;
    final int degOfC = cRe[0][0].length - 1;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        ComplexPolynomial.times( aRe[i][0], aIm[i][0], degOfA,
                                 bRe[0][j], bIm[0][j], degOfB,
                                 cRe[i][j], cIm[i][j], degOfC );

        for( int k=1; k<n; k++ ) {
          ComplexPolynomial.times( aRe[i][k], aIm[i][k], degOfA,
                                   bRe[k][j], bIm[k][j], degOfB,
                                   tmpRe,     tmpIm,     degOfC );
          
          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }

  private static void
  blockMultiplyPosNegReverse(
    final double [][][] aRe,   final double [][][] aIm,
    final double [][][] bRe,   final double [][][] bIm,
    final double [][][] cRe,   final double [][][] cIm,
    final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe || bRe == cRe) {
      throw new IllegalArgumentException();
    }

    final int n = aRe.length;

    final int degOfA = aRe[0][0].length - 1;
    final int degOfB = bRe[0][0].length - 1;
    final int degOfC = cRe[0][0].length - 1;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        ComplexPolynomial.times( aRe[0][j], aIm[0][j], degOfA,
                                 bRe[i][0], bIm[i][0], degOfB,
                                 cRe[i][j], cIm[i][j], degOfC );

        for( int k=1; k<n; k++ ) {
          ComplexPolynomial.times( aRe[k][j], aIm[k][j], degOfA,
                                   bRe[i][k], bIm[i][k], degOfB,
                                   tmpRe,     tmpIm,     degOfC );


          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }


  /** Special block Laurent series multiplication c = a * b
      Assumes lower and upper degrees
      a: (-n, n)
      b: (-n, 0)
      c: (0, n)
   */
  private static void
  blockMultiplyNegPos(
    final double [][][] aRe,   final double [][][] aIm,
    final double [][][] bRe,   final double [][][] bIm,
    final double [][][] cRe,   final double [][][] cIm,
    final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe || bRe == cRe) {
      throw new IllegalArgumentException();
    }

    final int n = aRe.length;
    final int N = bRe[0][0].length - 1;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        polynomialTimes( aRe[i][0], aIm[i][0],
                         bRe[0][j], bIm[0][j],
                         cRe[i][j], cIm[i][j], N );

        for( int k=1; k<n; k++ ) {

          polynomialTimes( aRe[i][k], aIm[i][k],
                           bRe[k][j], bIm[k][j],
                           tmpRe,     tmpIm,    N );

          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }

  private static void blockMultiplyNegPosReverse(
    final double [][][] aRe,   final double [][][] aIm,
    final double [][][] bRe,   final double [][][] bIm,
    final double [][][] cRe,   final double [][][] cIm,
    final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe || bRe == cRe) {
      throw new IllegalArgumentException();
    }

    final int n = aRe.length;
    final int N = bRe[0][0].length - 1;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        polynomialTimes( aRe[0][j], aIm[0][j],
                         bRe[i][0], bIm[i][0],
                         cRe[i][j], cIm[i][j], N );

        for( int k=1; k<n; k++ ) {

          polynomialTimes( aRe[k][j], aIm[k][j],
                           bRe[i][k], bIm[i][k],
                           tmpRe,     tmpIm,    N );


          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }


  /** Special Laurent series multipliction
      Assumes lower and upper degrees
      a: (-n, n)
      b: (-n, 0)
      c: (0, n)
      The negative part of a is not used.
      Assumes inputs are different.
   */
  private static void
  polynomialTimes(
    final double [] aRe, final double [] aIm,
    final double [] bRe, final double [] bIm,
    final double [] cRe, final double [] cIm, int N ) {

    if (aRe == cRe || bRe == cRe) {
      throw new IllegalArgumentException();
    }

    double cr, ci;

    for (int i = 0; i <= N; i++) {
      cr = 0.0;
      ci = 0.0;
      for (int j = i; j <= N; j++) {
        cr += aRe[N+j] * bRe[N+i-j] - aIm[N+j] * bIm[N+i-j];
        ci += aRe[N+j] * bIm[N+i-j] + aIm[N+j] * bRe[N+i-j];
      }
      cRe[i] = cr;
      cIm[i] = ci;
    }
  }


  /** Copies coeff into m
      (coeff[0] ... coeff[2N+1]) represent (a[-N], ... a[N])
      That is, a[p] = coeff[N+p]
   */
  private static void
  setBlock(
    final double [] coeff,
    final double [][] m,
    final int row, final int col,
    final int N, final int blockSize ) {

    final int rowOffset = row * blockSize;
    final int colOffset = col * blockSize;
    int R = blockSize;

    for( int j = 0; j < blockSize; j++ ) {
      int p0 = Math.max( 0, j-N );
      int p1 = Math.min( R, N+1+j );

      // initial zeros
      for( int i = 0; i < p0; i++ ) {
        m[i+rowOffset][j+colOffset] = 0.0;
      }
      // coefficients
      for( int i = p0; i < p1; i++ ) {
        m[i+rowOffset][j+colOffset] = coeff[N+i-j];
      }
      // final zeros
      for ( int i = p1; i < R; i++ ) {
        m[i+rowOffset][j+colOffset] = 0.0;
      }
    }
  }

  /** Extract x into coeff.
   */
  private static void
  getBlock(
    final double [] coeff, final double [] x,
    final int offset,
    final int N ) {
    System.arraycopy( x, offset, coeff, 0, N+1 );
  }

  /** Set all the entries of x to 0
   */
  private static void assignZero(
    final double [] x ) {
    for( int i = 0; i < x.length; i++ ) {
      x[i] = 0.0;
    }
  }

  public static void
  normalize(
    final double [][][] tRe,    final double [][][] tIm,
    final double [][][] uRe,    final double [][][] uIm,
    final double [][]   t0Re,      final double [][]   t0Im,
    final double [][]   q0Re,      final double [][]   q0Im,
    final double []     dRe,       final double []     dIm,
    final double []     pRe,       final double []     pIm,
    boolean reverse, boolean neg ) {

    if (tRe == uRe) {
      throw new IllegalArgumentException();
    }

    final int n = tRe.length;
    final int N = tRe[0][0].length - 1;

    // Put the constant term of t into t0
    int c = ( ! neg ? 0 : N );

    extractCoefficient( tRe, tIm, c, t0Re, t0Im );

    // compute the inverse of t0 in place
    Inversion.compute( t0Re, t0Im, q0Re, q0Im, dRe, dIm, pRe, pIm );
    // q0, d, p used temp.

    if (! reverse) {
      // u = t * q0
      times1( tRe, tIm, q0Re, q0Im, uRe, uIm, pRe, pIm );
    }
    else {
      // u = q0 * t
      times1( q0Re, q0Im, tRe, tIm, uRe, uIm, pRe, pIm );
    }
  }

  /** Copy the kth coefficient of b into bk */
  private static void
  extractCoefficient(
    final double [][][] bRe,
    final double [][][] bIm,
    final int k,
    final double [][] bkRe,
    final double [][] bkIm ) {

    final int n = bRe.length;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {
        bkRe[i][j] = bRe[i][j][k];
        bkIm[i][j] = bIm[i][j][k];
      }
    }
  }

  /** c = a * b
  */
  private static void
  times1( final double [][][] aRe,   final double [][][] aIm,
         final double [][]   bRe,   final double [][]   bIm,
         final double [][][] cRe,   final double [][][] cIm,
         final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe) {
      throw new IllegalArgumentException("a==c");
    }

    final int n = aRe.length;

    final int deg = aRe[0][0].length;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        VectorOperations.times( aRe[i][0], aIm[i][0],
                                bRe[0][j], bIm[0][j], cRe[i][j], cIm[i][j] );

        for( int k=1; k<n; k++ ) {

          VectorOperations.times( aRe[i][k], aIm[i][k],
                                  bRe[k][j], bIm[k][j], tmpRe, tmpIm );

          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }

  /** c = b * a
  */
  private static void
  times1(
    final double [][]   bRe,   final double [][]   bIm,
    final double [][][] aRe,   final double [][][] aIm,
    final double [][][] cRe,   final double [][][] cIm,
    final double []     tmpRe, final double []     tmpIm ) {

    if (aRe == cRe) {
      throw new IllegalArgumentException("a==c");
    }

    final int n = aRe.length;

    final int deg = aRe[0][0].length;

    for( int i = 0; i < n; i++ ) {
      for( int j = 0; j < n; j++ ) {

        VectorOperations.times( aRe[0][j], aIm[0][j],
                                bRe[i][0], bIm[i][0], cRe[i][j], cIm[i][j] );

        for( int k=1; k<n; k++ ) {

          VectorOperations.times( aRe[k][j], aIm[k][j],
                                  bRe[i][k], bIm[i][k], tmpRe, tmpIm );

          VectorOperations.plus( cRe[i][j], tmpRe, cRe[i][j] );
          VectorOperations.plus( cIm[i][j], tmpIm, cIm[i][j] );
        }
      }
    }
  }

}
