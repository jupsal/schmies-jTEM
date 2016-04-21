/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
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

package de.jtem.blas;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.Field;

abstract class AbstractMatrix implements Serializable,Cloneable
{
    private static final long serialVersionUID = 1L;

    // are entries equal ?    a "=" b  <=>  |a-b|^2 < EPSILON
    static double EPSILON = 1.e-28;

    public int numRows, numCols;

    /** return total number of entries  */
    public final int size()
    {
	return numRows * numCols;
    }

    /** return total number of entries  */
    public final int getNumEntries() {
	return numRows * numCols;
    }

    /** return number of rows  */
    public final int getNumRows() {
	return numRows;
    }
    
    /** set number of rows,  set possible new Entries to zero */
    public final void setNumRows( int aNumRows) {
	resize( aNumRows, numCols );
    }

    /** return number of columns  */
    public final int getNumCols() {
	return numCols;
    }
 
    /**  set number of columns,  set possible new Entries to zero */
    public final void setNumCols( int aNumCols) {
	resize( numRows, aNumCols );
    }

    /** throw <b>IllegalArgumentException</b> if <i>anIndex</i> exceeds number of rows */
    public final void checkRowIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= numRows )
	    throw new IllegalArgumentException("row index "+anIndex+" is not in [0,"+numRows+"[");
    }

    /** throw <b>IllegalArgumentException</b> if <i>anIndex</i> exceeds number of columns */
    public final void checkColIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= numCols )
	    throw new IllegalArgumentException("col index "+anIndex+" is not in [0,"+numCols+"[");
    }

    /** throw <b>IllegalArgumentException</b> if <i>numOfRows</i> exceeds number of rows */
    public final void checkNumRows( int numOfRows ) {
	if( numOfRows != numRows ) 
	    throw new IllegalArgumentException("wrong shape: num of rows is: "+numRows+", but should: "+numOfRows );
    }

    /** throw <b>IllegalArgumentException</b> if <i>numOfCols</i> exceeds number of columns */
    public final void checkNumCols( int numOfCols ) {
	if( numOfCols != numCols ) 
	    throw new IllegalArgumentException("wrong shape: num of cols is: "+numCols+", but should: "+numOfCols );
    }

    /** throw <b>IllegalArgumentException</b> if matrix' shape is different from <i>numOfRows*numOfCols</i> */
    public final void checkShape( int numOfRows, int numOfCols ) {
	if( numOfRows != numRows ) 
	    throw new IllegalArgumentException("wrong shape: num of rows is: "+numRows+", but should: "+numOfRows );
	if( numOfCols != numCols ) 
	    throw new IllegalArgumentException("wrong shape: num of cols is: "+numCols+", but should: "+numOfCols );
    }

    /** throw <b>IllegalArgumentException</b> if <i>M</i>'s shape is different from <i>this</i> */
    public final void checkShape( AbstractMatrix M ) {
	checkShape( M.numRows, M.numCols );
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( AbstractMatrix M ) {
	return M.numRows == numRows && M.numCols == numCols;
    }


    /** throw <b>IllegalArgumentException</b> if <i>M</i>'s shape is different from <i>this</i> */
    public final void checkShape( double[][] M ) {
	checkArray( M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( double[][] M ) {
	return hasMatrixShape(M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> matches a newNumRows*newNumCols shaped rectangle, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( double[][] M, int newNumRows, int newNumCols ) 
    {
      boolean result = (M.length == newNumRows);
      if (result) 
	for (int i=0;i<newNumRows;i++) result &= (M[i].length == newNumCols);
      return result;
    }

    /** return <tt>true</tt> if <i>M</i> is rectangular, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( double[][] M ) 
    {
      if (M.length > 0) return hasMatrixShape(M,M.length,M[0].length);
      return true;
    }

    /**  check whether shape of <i>M</i> matches a newNumRows*newNumCols shaped rectangle  */
    public final static void checkArray( double[][] M, int newNumRows, int newNumCols ) 
    {
      if (!hasMatrixShape(M,newNumRows,newNumCols))
	throw new IllegalArgumentException("Array has wrong size !");
    }

    /**  check whether shape of <i>M</i> mathces <i>this</i> matrix */
    public final void checkArray( double[][] M ) 
    {
      if (M.length > 0) checkArray(M, M.length, M[0].length);
    }

    /** throw <b>IllegalArgumentException</b> if <i>M</i>'s shape is different from <i>this</i> */
    public final void checkShape( int[][] M ) {
	checkArray( M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( int[][] M ) {
	return hasMatrixShape(M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> matches a newNumRows*newNumCols shaped rectangle, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( int[][] M, int newNumRows, int newNumCols ) 
    {
      boolean result = (M.length == newNumRows);
      if (result) 
	for (int i=0;i<newNumRows;i++) result &= (M[i].length == newNumCols);
      return result;
    }

    /** return <tt>true</tt> if <i>M</i> is rectangular, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( int[][] M ) 
    {
      if (M.length > 0) return hasMatrixShape(M,M.length,M[0].length);
      return true;
    }

    /**  check whether shape of <i>M</i> matches a newNumRows*newNumCols shaped rectangle  */
    public final static void checkArray( int[][] M, int newNumRows, int newNumCols ) 
    {
      if (!hasMatrixShape(M,newNumRows,newNumCols))
	throw new IllegalArgumentException("Array has wrong size !");
    }

    /**  check whether shape of <i>M</i> mathces <i>this</i> matrix */
    public final void checkArray( int[][] M ) 
    {
      if (M.length > 0) checkArray(M, M.length, M[0].length);
    }

    /** return <tt>true</tt> if <i>this</i> is a square matrix, <tt>false</tt> o/w */
    public final boolean isSquared() {
	return numRows == numCols;
    }

    /**  check whether <i>this</i> is a square matrix */
    public final void checkSquared() {
	if (!isSquared()) throw new IllegalArgumentException("matrix is NOT squared");
    }

    /** @deprecated use hasMatrixShape( double[][] M )*/
    public final static boolean isRectangular( double[][] M ) 
    {
      return hasMatrixShape(M);
    }
    /** @deprecated please tell me if this occurs */
    public final static void makeRectangular( double [][] M ) 
    {
      int aNumRows = M.length;
      
      if (aNumRows>0) {
	int aNumCols = M[0].length;
	int maxLength = aNumCols;
	
	for( int i = 1; i<aNumRows; i++ )
	  if( M[i].length > maxLength )
	    maxLength = M[i].length; 
	
	if( maxLength != aNumCols )
	  for( int i = 0; i<aNumRows; i++ )
	    if( M[i].length < maxLength ) {
	      double [] newM = new double[maxLength];
	      
	      System.arraycopy( M[i], 0, newM, 0, M[i].length );
	      
	      M[i] = newM;
	    }
	
      }
    }

    abstract void newSize( int newNumRows, int newNumCols, boolean copyValues );


    /** resize matrix to <i>newNumRows*newNumCols</i>, but ignore the values.*/
    public void newSize( int newNumRows, int newNumCols ) {
	newSize( newNumRows, newNumCols, false );	
    }
    /** resize matrix to <i>newNumRowsAndCols*newNumRowsAndCols</i>, but ignore the values.*/
    public final void newSize( int newNumRowsAndCols ) {
	newSize( newNumRowsAndCols, newNumRowsAndCols, false );
    }

    /** resize matrix to <i>aMatrix</i>'s shape, but ignore the values. */
    public final void newSize( AbstractMatrix aMatrix ) {
	newSize( aMatrix.numRows, aMatrix.numCols );
    }
    
    /** resize matrix to <i>newNumRows*newNumCols</i>
         w/o changing its (surviving) entries, possibly adding zero entries */
    public void resize( int newNumRows, int newNumCols ) {
	newSize( newNumRows, newNumCols, true );
    }

    /** resize matrix to a square of size <i>newNumRowsAndCols*newNumRowsAndCols</i>
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public final void resize( int newNumRowsAndCols ) {
	newSize( newNumRowsAndCols, newNumRowsAndCols, true );
    }

    /** resize matrix to <i>aMatrix</i>'s shape
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public final void resize( AbstractMatrix aMatrix ) {
	newSize( aMatrix.numRows, aMatrix.numCols, true );
    }
    
    /** resize matrix to <i>newNumRows*newNumCols</i>
        w/o changing its (surviving) entries, possibly adding zero entries          */
    public void setSize( int newNumRows, int newNumCols ) {
        newSize(newNumRows,newNumCols, true);
    }

    /** resize matrix to a square of size <i>newNumRowsAndCols*newNumRowsAndCols</i>
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public void setSize( int newNumRowsAndCols ) {
	newSize( newNumRowsAndCols, newNumRowsAndCols, true );
    }

    /** resize matrix to <i>aMatrix</i>'s shape
        w/o changing its (surviving) entries, possibly adding zero entries               */
    public final void setSize( AbstractMatrix aMatrix ) {
	newSize( aMatrix.numRows, aMatrix.numCols, true );
    }


    /** return the total sum of squares of entries' length */
    public abstract double normSqr();

    /** return <tt>true</tt> if <i>this</i> has only zero entries */
    public boolean isZero() {
      return (normSqr() <= numRows*numCols*EPSILON);
    }
    
    /** check whether <i>M</i>'s shape is different from <i>this</i> */
    public final void checkShape( final Complex[][] M ) {
	checkArray( M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( final Complex[][] M ) {
	return hasMatrixShape(M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> matches a newNumRows*newNumCols shaped rectangle, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( final Complex[][] M, final int newNumRows, final int newNumCols ) 
    {
      boolean result = (M.length == newNumRows);
      if (result) 
	for (int i=0;i<newNumRows;i++) result &= (M[i].length == newNumCols);
      return result;
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( final ComplexConstant[][] M ) {
	return hasMatrixShape(M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> matches a newNumRows*newNumCols shaped rectangle, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( final ComplexConstant[][] M, final int newNumRows, final int newNumCols ) 
    {
      boolean result = (M.length == newNumRows);
      if (result) 
	for (int i=0;i<newNumRows;i++) result &= (M[i].length == newNumCols);
      return result;
    }

    /** return <tt>true</tt> if <i>M</i> and <i>this</i> have same shape, <tt>false</tt> o/w */
    public final boolean hasSameShape( final Field.Complex[][] M ) {
	return hasMatrixShape(M,numRows,numCols);
    }

    /** return <tt>true</tt> if <i>M</i> matches a newNumRows*newNumCols shaped rectangle, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( final Field.Complex[][] M, final int newNumRows, final int newNumCols ) 
    {
      boolean result = (M.length == newNumRows);
      if (result) 
	for (int i=0;i<newNumRows;i++) result &= (M[i].length == newNumCols);
      return result;
    }

    /** return <tt>true</tt> if <i>M</i> is rectangular, <tt>false</tt> o/w */
    public final static boolean hasMatrixShape( final Complex[][] M ) 
    {
      if (M.length > 0) return hasMatrixShape(M,M.length,M[0].length);
      return true;
    }

    /**  check whether shape of <i>M</i> matches a newNumRows*newNumCols shaped rectangle  */
    public final static void checkArray( final Complex[][] M, final int newNumRows, final  int newNumCols ) 
    {
      if (!hasMatrixShape(M,newNumRows,newNumCols))
	throw new IllegalArgumentException("Array has wrong size !");
    }

   /**  check whether shape of <i>M</i> matches a newNumRows*newNumCols shaped rectangle  */
    public final static void checkArray( final Field.Complex[][] M, final int newNumRows, final  int newNumCols ) 
    {
      if (!hasMatrixShape(M,newNumRows,newNumCols))
	throw new IllegalArgumentException("Array has wrong size !");
    }

    /**  check whether shape of <i>M</i> matches a newNumRows*newNumCols shaped rectangle  */
    public final static void checkArray( final ComplexConstant[][] M, final int newNumRows, final  int newNumCols ) 
    {
      if (!hasMatrixShape(M,newNumRows,newNumCols))
	throw new IllegalArgumentException("Array has wrong size !");
    }

   /**  check whether shape of <i>M</i> is rectangular  */
    public final void checkArray( final Complex[][] M ) 
    {
      if (M.length > 0) checkArray(M, M.length, M[0].length);
    }

    /**  check whether shape of <i>M</i> is rectangular  */
    public final void checkArray( final Field.Complex[][] M ) 
    {
      if (M.length > 0) checkArray(M, M.length, M[0].length);
    }

    /**  check whether shape of <i>M</i> is rectangular  */
    public final void checkArray( final ComplexConstant[][] M ) 
    {
      if (M.length > 0) checkArray(M, M.length, M[0].length);
    }

    /** return a clone */
    public abstract Object clone();

    /** print this matrix */
    public void print( String title ) {
	
	System.out.println( title );

	System.out.println( toString() );
    }
}
