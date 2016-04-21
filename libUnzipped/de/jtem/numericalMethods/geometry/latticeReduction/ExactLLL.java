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

package de.jtem.numericalMethods.geometry.latticeReduction;

public final class ExactLLL {

  /* bei umstellung ist zu beachten:

     1. indexierung geht von 1..n und musy zu 0..n-1           geaendert werden
     2. pivot-element ist markiert durch 0 und musz zu -1      geaendert werden
     3. D[x] durch D[1+x] ersaetzen
  */

    private static long exactDiv( long a, long b ) {
    
	if( a % b != 0 )
	    throw new RuntimeException( "nonbzero remainder");
	
	return a / b;
    }
    
    //*  rounds a/d to nearest integer, breaking ties
    //    by rounding towards zero.  Assumes d > 0.
    private static long balDiv( long a, long b ) {
	final long r = a % b;
	final long q = a / b;
	
	if( r+r > b || r+r == b && q < 0 )
	    return q+1;
	
	return q;
    }
    
    //* c = (x*c1 + y*c2)/z
    private static long mulAddDiv( long c1, long c2, long x, long y, long z ) {
	return (x*c1 + y*c2)/z;
    }

    //* c = (x*c1 - y*c2)/z
    private static long mulSubDiv( long c1, long c2, long x, long y, long z ) {
	return (x*c1 - y*c2)/z;
    }

    //* c = (x*c1 - y*c2)/z
    private static void mulSubDiv( long [] c, long [] c1, long [] c2, long x, long y, long z ) {
	int n = c.length;
	
	if( c1.length != n || c2.length != n )
	    throw new IllegalArgumentException( "length mismatch" );
	
	for( int i=0; i<n; i++ )
	    c[i] = ( x*c1[i] - y*c2[i] ) / z;
    }
    
    //* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2)
    private static void rowTransform( long [] c1, long [] c2, long x, long y, long u, long v ) {
	int n = c1.length;
	
	if( c2.length != n )
	    throw new IllegalArgumentException( "length mismatch" );
	
	for( int i=0; i<n; i++ ) {
	    
	    long t1 = x * c1[i] + y * c2[i];
	    long t2 = u * c1[i] + v * c2[i];
	    
	    c1[i] = t1;
	    c2[i] = t2;
	}
    }
    
    //* (M[i][k], M[j][k]) = ( x*M[i][k] + y*M[j][k], u*M[i][k] + v*M[j][k] )
    private static void rowTransform( long [][] M, int i, int j, int k, 
				      long x, long y, long u, long v ) {
	
    	    long t1 = x * M[i][k] + y * M[j][k];
	    long t2 = u * M[i][k] + v * M[j][k];

	    M[i][k] = t1;
	    M[j][k] = t2;
    }
  
    //* c = c1 + y*c2
    private static long mulAdd( long c1, long c2, long y ) {
	return c1 + y*c2;
    }

    // c = c1 - x*c2
    private static long mulSub( long c1, long c2, long y ) {
	return c1 - y*c2;
    }


    //* c = c1 - x*c2
    private static  void mulSub( long [] c, long [] c1, long [] c2, long x ) {
	int n = c.length;

	if( c1.length != n || c2.length != n )
	    throw new IllegalArgumentException( "length mismatch" );
	
	for( int i=0; i<n; i++ )
	    c[i] = c1[i] - x*c2[i];
    }
      
    //* test if a*d1^2 > b*(d0*d2 + lam^2)    
    private static boolean swapTest( long d0, long d1, long d2, long lam, long a, long b ) {
	return a * d1 * d1 > b * ( d0 * d2 + lam * lam );
    }
   
    private static void reduce( int k, int l, 
				long [][] B, int [] P, long [] D, long [][] lam, long [][] U ) {

	if( P[l] == -1 ) return;

	if( 2 * Math.abs( lam[k][P[l]] ) < D[1+P[l]] ) return;

	long r = balDiv( lam[k][P[l]], D[1+P[l]] );

	mulSub( B[k], B[k], B[l], r );

	if( U != null ) 
	    mulSub( U[k], U[k], U[l], r );

	for( int j = 0; j <= l-1; j++)
	    if( P[j] != -1 )
		lam[k][P[j]] = mulSub( lam[k][P[j]], lam[l][P[j]], r );

	lam[k][P[l]] = mulSub( lam[k][P[l]], D[1+P[l]], r);
    }

    private static void swap( long [] a, long [] b ) {
	long [] tmp = a;

	b = a;
	b = tmp;
    }

  private static final long gcd(long dvr,long dvd) {
    long rem;
    while (dvr!=0) {
      rem=dvd%dvr;
      dvd=dvr;
      dvr=rem;
    }
    return dvd;
  }

  private static final long xgcd( long a, long b, long [] s ) {

    long r0, r1 = a, r2 =   b%a;
    long x0, x1 = 1, x2 = - b/a;
    long y0, y1 = 0, y2 =   1;

    do {
      r0 = r1; r1 = r2;
      x0 = x1; x1 = x2;
      y0 = y1; y1 = y2;

      long q = r0 / r1;

      r2 = r0 % r1;

      x2 = x0 - q * x1;
      y2 = y0 - q * y1; 
      
    } while ( r2 != 0 );
    
    s[0] = x1;
    s[1] = y1; 

    return r1;
  }

    /* swaps vectors k-1 and k;  assumes P(k-1) != -1
       returns true if vector k-1 need to be reduced after the swap...
       this only occurs in 'case 2' when there are linear dependencies */

    private static boolean swap( int k, 
				 long [][] B, int [] P, long [] D, long [][] lam, long [][] U, 
				 int m, boolean verbose ) {


	if( P[k] != -1 ) {
	    if( verbose ) System.out.println( "swap case 1: " + k );
	    	    
	    swap(B[k-1], B[k]);

	    if( U != null ) swap( U[k-1], U[k] );
   
	    for( int j = 0; j <= k-2; j++)
		if( P[j] != -1 ) {
		    long tmp = lam[k-1][P[j]];
		    lam[k-1][P[j]] = lam[k][P[j]];
		    lam[k]  [P[j]] = tmp;
		} 

	    for( int i = k+1; i < m; i++) {
		long t1 = mulAddDiv(lam[i][P[k]-1], lam[i][P[k]], lam[k][P[k]-1], D[1+P[k]-2], D[1+P[k]-1] ); 
		long t2 = mulSubDiv(lam[i][P[k]-1], lam[i][P[k]], D[1+P[k]],   lam[k][P[k]-1], D[1+P[k]-1] );
		
		lam[i][P[k]-1] = t1;
		lam[i][P[k]  ] = t2;
	    }

	    D[1+P[k]-1] = mulAddDiv( D[1+P[k]], lam[k][P[k]-1], D[1+P[k]-2], lam[k][P[k]-1], D[1+P[k]-1]);
	    
	    return false;
	
	} else if( lam[k][P[k-1]] != 0 ) {

	    if( verbose ) System.out.println( "swap case 2: " + k );

	    long [] xy = new long[2];

	    long e = xgcd( lam[k][P[k-1]], D[1+P[k-1]], xy );

	    long x = xy[0];
	    long y = xy[1];

	    long t1 = exactDiv( lam[k][P[k-1]], e );
	    long t3 = exactDiv( D[1+P[k-1]],    e );
	    long t2 = -t3;

	    rowTransform(B[k-1], B[k], t1, t2, y, x);

	    if (U != null) 
	      rowTransform( U[k-1], U[k], t1, t2, y, x );

	    for( int j = 0; j <= k-2; j++)
	      if( P[j] != -1 )
		  rowTransform( lam, k-1, k, P[j], t1, t2, y, x );
  
	    t2 *= t2;

	    D[1+P[k-1]] = exactDiv( D[1+P[k-1]], t2);

	    for( int i = k+1; i < m; i++)
		if( P[i] != -1 ) {
		    D[1+P[i]] = exactDiv( D[1+P[i]], t2);
		    for( int j = i+1; j < m; j++) {
			lam[j][P[i]] = exactDiv( lam[j][P[i]], t2);
		    }
		}

	    for( int i = k+1; i < m; i++) {
		lam[i][P[k-1]] = exactDiv( lam[i][P[k-1]], t3);
	    }

	    int tmp = P[k-1];
	    P[k-1] = P[k];
	    P[k  ] = tmp;

	    return true;
	
	} else {

	    if( verbose ) System.out.println( "swap case 3: " + k );

	    swap(B[k-1], B[k]);

	    if( U != null ) swap( U[k-1], U[k] );
   
	    for ( int j = 0; j <= k-2; j++)
		if( P[j] != -1 ) {
		    long tmp = lam[k-1][P[j]];
		    lam[k-1][P[j]] = lam[k][P[j]];
		    lam[k  ][P[j]] = tmp;
		}  

	    int tmp = P[k-1];
	    P[k-1] = P[k];
	    P[k  ] = tmp;

	    return false;
	}
    }

    private static long innerProduct( long [] v, long [] w ) {
	long p = 0;

	if( v.length != w.length )
	    throw new IllegalArgumentException( "length mismatch" );

	for( int i=0; i<v.length; i++ )
	    p+= v[i]*w[i];
	
	return p;
    }

    //* returns the altered s
    private static int incrementalGS( long [][] B, int [] P, long [] D, long [][] lam, int s, int k ) {

	long n = B[0].length;
	long m = B.   length;

	for( int j = 0; j <= k-1; j++) {
	    int posj = P[j];
	    
	    if( posj == -1 ) continue;
	    
	    long u = innerProduct( B[k], B[j] );
	    
	    for( int i = 0; i <= posj-1; i++) {
		long t1 = D[1+i]*u;
		long t2 = lam[k][i] * lam[j][i];
		
		u = (t1 - t2) / D[1+i-1];
	    }
	    
	    lam[k][posj] = u;
	}

	long u = innerProduct( B[k], B[k] );
	
	for( int i = 0; i < s; i++) {
	    long t1 = D[1+i] * u;
	    long t2 = lam[k][i] * lam[k][i];
	    
	    u = (t1 - t2) / D[1+i-1];
	}
	
	if (u == 0) {
	    P[k] = -1;
	}
	else {
	    P[k] = s;
	    D[1+s] = u;

	    s++;
	}

	return s;
    }

    private static int  []   P   = new int [10];
    private static long []   D   = new long[11];
    private static long [][] lam = new long[10][10];
    
    private static long [] reduce( long [][] B, long [][] U, long a, long b, boolean verbose ) {

	int n = B[0].length;
	int m = B.   length;

	boolean force_reduce = true;

	if( P.length < m ) {
	    P = new int[m];
	    D = new long[m+1];

	} else for( int i=0; i<m; i++ ) {
	    P[i] = 0;
	    D[i] = 0;
	    for( int j=0; j<m; j++ )
		lam[i][j] = 0;
	}

	D[m] = 0;
	D[0] = 1;

	if( U.length != m || U[0].length != m )
	    throw new IllegalArgumentException( "U has wrong length" );

	if( U != null )
	    for( int i=0; i<m; i++ )
		for( int j=0;j<m; j++ )
		    U[i][j] = i==j ? 1 : 0;


	int s     = 0;
	int k     = 0;
	int max_k =-1;

	while( k < m ) {
	    if (k > max_k) {
		s = incrementalGS(B, P, D, lam, s, k);
		max_k = k;
	    }

	    if( k == 0 ) {
		force_reduce = true; 
		k++;
		continue;
	    }

	    if( force_reduce )
		for( int j = k-1; j >= 0; j-- )
		    reduce(k, j, B, P, D, lam, U);

	    if( P[k-1] != -1 && 
		( P[k] == -1 || swapTest(D[1+P[k]], D[1+P[k]-1], D[1+P[k]-2], lam[k][P[k]-1], a, b))) {
		force_reduce = swap(k, B, P, D, lam, U, max_k, verbose);
		k--;
		
	    } else {
		force_reduce = true;
		k++;
	    }
	}
	
	long [] finalD = new long[s+1];
	
	System.arraycopy( D, 0, finalD, 0, s+1 );

	return finalD;
    }

    public static int reduce( long [][] B, long[][] U ) {
	return reduce( B, U, 3, 4, false ).length - 1;
    }

    public static int reduce( long [][] B, long[][] U, long [] det ) {

	long [] D = reduce( B, U, 3, 4, false );
	
	int s = D.length - 1;

	det[0] = D[s];

	return s;
    }

    public static int reduce( long [][] B, long[][] U, long a, long b, long [] det ) {

	if (a <= 0 || b <= 0 || a > b || b/4 >= a) 
	    throw new IllegalArgumentException( "wrong constants");

	long [] D = reduce( B, U, a, b, false );
	
	int s = D.length - 1;

	det[0] = D[s];

	return s;
    }

    public static int reduce( long [][] B, long a, long b, long [] det ) {
	return reduce( B, null, a, b, det );
    }

    /* dimension of image of the lattice
       returns the det of the lattice
     */
    private static long image( long [][] B, long [][] U, long [] det, boolean verbose ) {
  
	int n = B[0].length;
	int m = B.   length;

	boolean force_reduce = true;

	if( P.length < m ) {
	    P = new int[m];
	    D = new long[m+1];

	} else for( int i=0; i<m; i++ ) {
	    P[i] = 0;
	    D[i] = 0;
	    for( int j=0; j<m; j++ )
		lam[i][j] = 0;
	}

	D[m] = 0;
	D[0] = 1;

	if( U.length != m || U[0].length != m )
	    throw new IllegalArgumentException( "U has wrong length" );

	if( U != null )
	    for( int i=0; i<m; i++ )
		for( int j=0;j<m; j++ )
		    U[i][j] = i==j ? 1 : 0;

	int s     = 0;
	int k     = 0;
	int max_k =-1;

	while (k < m) {
	    if (k > max_k) {
		incrementalGS(B, P, D, lam, s, k);
		max_k = k;
	    }

	    if (k == 0) {
		force_reduce = true; 
		k++;
		continue;
	    }

	    if (force_reduce)
		for( int j = k-1; j >= 0; j-- ) 
		    reduce(k, j, B, P, D, lam, U);

	    if (P[k-1] != -1 && P[k] == -1 ) {
		force_reduce = swap(k, B, P, D, lam, U, max_k, verbose);
		k--;
	    }
	    else {
		force_reduce = true;
		k++;
	    }
	}

	det[0] = D[s];
	
	return s;
    }

    /* dimension of image of the lattice
       returns the det of the lattice
    */
    public static long image( long [][] B, long [][] U, long [] det ) {
  
	return image( B, U, det, false );
    }

    static long [] det = new long[1];

    public static long image( long [][] B, long [][] U ) {

	return image( B, U, det, false );
    }

    public static long image( long [][] B ) {

	return image( B, null, det, false );
    }






    /*

long LatticeSolve(vec_long& x, const mat_long& A, const vec_long& y, long reduce)
{
   long n = A.NumRows();
   long m = A.NumCols();

   if (y.length() != m)
      Error("LatticeSolve: dimension mismatch");

   if (reduce < 0 || reduce > 2)
      Error("LatticeSolve: bad reduce parameter");

   if (IsZero(y)) {
      x.SetLength(n);
      clear(x);
      return 1;
   }

   mat_long A1, U1;
   long det2;
   long im_rank, ker_rank;

   A1 = A;

   im_rank = image(det2, A1, U1);
   ker_rank = n - im_rank;

   mat_long A2, U2;
   long new_rank;
   long i;

   A2.SetDims(im_rank + 1, m);
   for (i = 1; i <= im_rank; i++)
      A2(i) = A1(ker_rank + i);

   A2(im_rank + 1) = y;

   new_rank = image(det2, A2, U2);

   if (new_rank != im_rank || 
      (U2(1)(im_rank+1) != 1  && U2(1)(im_rank+1) != -1))
      return 0;

   vec_long x1;
   x1.SetLength(im_rank);

   for (i = 1; i <= im_rank; i++)
      x1(i) = U2(1)(i);

   if (U2(1)(im_rank+1) == 1)
      negate(x1, x1);

   vec_long x2, tmp;
   x2.SetLength(n);
   clear(x2);
   tmp.SetLength(n);

   for (i = 1; i <= im_rank; i++) {
      mul(tmp, U1(ker_rank+i), x1(i));
      add(x2, x2, tmp);
   }

   if (reduce == 0) {
      x = x2;
      return 1;
   }
   else if (reduce == 1) {
      U1.SetDims(ker_rank+1, n);
      U1(ker_rank+1) = x2;
      image(det2, U1);

      x = U1(ker_rank + 1);
      return 1;
   }
   else if (reduce == 2) {
      U1.SetDims(ker_rank, n);
      LLL(det2, U1);
      U1.SetDims(ker_rank+1, n);
      U1(ker_rank+1) = x2;
      image(det2, U1);

      x = U1(ker_rank + 1);
      return 1;
   }

   return 0;
} 

*/




}
