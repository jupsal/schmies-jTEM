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

import de.jtem.mfc.field.Field;
import de.jtem.numericalMethods.algebra.linear.MatrixOperations;

/**
 * This class represents general integer matrices.
 * Following the number package's philosophy, operations on int matrices 
 * generally occur three-fold: For instances a,b,c of this class, 
 * <table> 
 * <tr><td valign="up"><tt>a.plus(b)</tt></td>
 *     <td>creates and returns a new instance containing the sum
 *         <i>a+b</i></td></tr> 
 * <tr><td valign="up"><tt>a.assignPlus(b)</tt></td><td>stores the sum 
 *         <i>a+b</i> in <i>a</i></td></tr>
 * <tr><td valign="up"><tt>c.assignPlus(a,b)</tt></td><td>stores the sum
 *         <i>a+b</i> in c</td></tr> 
 * </table> <p>
 * Non-commutative operations are <em>right operations</em>: 
 * <table>
 * <tr><td valign="up"></td><td><tt>a.assignTimes()</tt></td>
 *     <td> assigns <i>a</i> to <i>a*b</i>,</td></tr>
 * <tr><td valign="up">use</td><td><tt>a.assignTimes(b,a)</tt></td>
 *     <td> to assign <i>a</i> to <i>b*a</i></td></tr> 
 * </table> <p>
 * The array containing the matrixes content can be accessed directly using
 * <tt>re()</tt> to get and <tt>re(double[][])</tt>, to set <em>reference (!)</em>.
 *
 * @author	marina,vitali
 */
public final class IntegerMatrix extends AbstractMatrix {
    private static final long serialVersionUID = 1L;
    public int[] [] re;

    /** create a zero matrix */
    public IntegerMatrix() {
    	this( 0,0);
    }

    /** create a matrix with all entries equal to zero */
    public IntegerMatrix( int newNumRows, int newNumCols ) {
        numRows = newNumRows;
        numCols = newNumCols;
	re = new int[numRows] [numCols];
    }

    /** create a matrix and initialize all entries with <i>initValue</i> */
    public IntegerMatrix(int numrows, int numcols, int initValue) {
        this(numrows, numcols);
        assign(initValue);
    }

    /** create copy of matrix A */
    public IntegerMatrix( IntegerMatrix A ) {
        this(A.numRows, A.numCols);
	MatrixOperations.assign( A.re, re );
    }

    /** create a quadratic matrix with all entries equal to zero */
    public IntegerMatrix(int numrows) {
        this(numrows, numrows);
    }

    /**
     * create a Matrix consisting only of one row or column: If <i>covariant</i> is <tt>true</tt>, the matrix has one row
     * equal to <i>V</i>, otherwise it has one column equal to <i>V</i>
     */
    public IntegerMatrix( final IntegerVector V, boolean covariant) {
        this(covariant ? 1 : V.size(), covariant ? V.size() : 1);
        if (covariant) {
            System.arraycopy(V.re, 0, re[0], 0, numCols);
        } else {
            final int[] VRe = V.re;
            for (int i = 0; i < numRows; i++)
                re[i] [0] = VRe[i];
        }
    }

    /** create a quadratic matrix with its diagonal initialized to the vector <i>diag</i> */
    public IntegerMatrix( final IntegerVector diag ) {
        this( diag.size() );
	setDiagonal ( diag );
    }

    /** create a matrix with initialized with <i>m</i> if possible */
    IntegerMatrix(int[] [] m, boolean performArrayCheck) {
        assign(m, performArrayCheck);
    }

    /** create a matrix with initialized with <i>m</i> if possible */
    public IntegerMatrix( final int[] [] m )  {
        this(m, true);
    }

    void newSize(int newNumRows, int newNumCols, boolean copyValues) {
        int[] [] newRe = null;
        if (newNumCols != numCols) {
            newRe = new int[newNumRows] [newNumCols];
            if (copyValues) {
                int minNumCols = newNumCols < numCols ? newNumCols : numCols;
                for (int i = 0; i < newNumRows && i < numRows; i++) {
                    System.arraycopy(re[i], 0, newRe[i], 0, minNumCols);
                }
            }
        } else if (newNumRows != numRows) {
            newRe = new int[newNumRows] [];
            int i;
            for (i = 0; i < newNumRows && i < numRows; i++) {
                newRe[i] = re[i];
            }
            for ( ; i < newNumRows; i++) {
                newRe[i] = new int[newNumCols];
            }
        } else
            return;
        re = newRe;
        numCols = newNumCols;
        numRows = newNumRows;
    }

    /** return entry at <i>rowIndex,colIndex</i> */
    public int get(int rowIndex, int colIndex) {
        return re[rowIndex] [colIndex];
    }

    /** store <i>c</i> in entry <i>rowIndex,colIndex</i> */
    public void set(int rowIndex, int colIndex, int c) {
        re[rowIndex] [colIndex] = c;
        
    }

    /** assign all entries to <i>c</i> */
    public void assign(int c) {
        MatrixOperations.assign(re, c);
    }

    /** assign matrix to <i>M</i> */
    public void assign( final IntegerMatrix M ) {
	if (this == M) return ;
	newSize(M.numRows, M.numCols);
	MatrixOperations.assign(M.re,re);
    }
    
    void assign( int[] [] m, boolean performArrayCheck ) {
        if (m.length == 0) {
            newSize(0, 0);
        } else {
            if (performArrayCheck)
                checkArray(m);
            newSize(m.length, m[0].length);
	    MatrixOperations.assign(m,re);
	}
    }

    /** assign matrix to <i>M</i> */
    public void assign( final int[][] m ) {
        assign(m, true);
    }

    /** set <i>M</i> to <i>this</i> matrix */
    public void get( RealMatrix M ) {
        checkShape(M);
        M.newSize(M.numRows, M.numCols); // set all entries to zero
	MatrixOperations.assign(re,M.re);
    }

    /** The following two methods return/set REFERENCES only ! */

    /** return the REFERENCE to <i>this</i>' array */
    public int[] [] re() {
        return re;
    }

    /** set <i>this</i>' array by REFERENCING <i>M</i> */
    public void re(int[] [] M) {
        checkShape(M);
        re = M;
    }

    /** create and return new <b>de.jtem.blas.IntegerVector</b> containing row <i>aRowNum</i> */
    public IntegerVector getRow(int aRowNum) {
        return getRow(aRowNum, null);
    }

    /**
     * store row <i>aRowNum</i> in <i>V</i> and return <i>V</i>. If <i>V</i> equals null an IntegerVector
     * containing the row is created.
     */
    public IntegerVector getRow( int aRowNum, IntegerVector V ) {
        checkRowIndex(aRowNum);
        if (V == null)
            V = new IntegerVector (numCols);
        else
            V.newSize(numCols);
        System.arraycopy(re[aRowNum], 0, V.re, 0, numCols );
        return V;
    }

   /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( int aRowNum, final IntegerVector V ) {
        checkRowIndex(aRowNum);
        if (V.size() != numCols)
            throw new IllegalArgumentException("size of vector does not match with number of colums");
        System.arraycopy(V.re, 0, re[aRowNum], 0, numCols);
    }

     /** store <i>V</i> in row <i>aRowNum</i> */
    public void setRow( int aRowNum, final int [] V ) {
        checkRowIndex(aRowNum);
        if (V.length != numCols)
            throw new IllegalArgumentException("size of vector does not match with number of colums");
        System.arraycopy(V, 0, re[aRowNum], 0, numCols);
    }

   
    /** create and return new <b>de.jtem.blas.IntegerVector</b> containing column <i>aColNum</i> */
    public IntegerVector getCol( int aColNum ) {
        return getCol(aColNum, null);
    }

    /**
     * store column <i>aColNum</i> in <i>V</i> and return <i>V</i>. If <i>V</i> equals null an IntegerVector
     * containg the column is created.
     */
    public IntegerVector getCol( int aColNum, IntegerVector V ) {
        checkColIndex(aColNum);
        if (V == null)
            V = new IntegerVector (numRows);
        else
            V.newSize(numRows);
        final int[] VRe = V.re;
        for (int i = 0; i < numRows; i++) {
            VRe[i] = re[i] [aColNum];
        }
        return V;
    }

    /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( int aColNum, final IntegerVector V ) {
        checkColIndex(aColNum);
        if (V.size() != numRows)
            throw new IllegalArgumentException("size of vector does not match with number of rows");
        final int[] VRe = V.re;
	MatrixOperations.assignCol( re, VRe, aColNum );        
    }

    /** store <i>V</i> in column <i>aColNum</i> */
    public void setCol( int aColNum, final int[] V ) {
        checkColIndex(aColNum);
        if (V.length != numRows)
            throw new IllegalArgumentException("size of vector does not match with number of rows");
	MatrixOperations.assignCol( re, V, aColNum );        
    }

    /** create and return a block of size to <i>(numOfRows,numOfCol)</i> from <i>(rowIndex,colIndex)</i> */
    public IntegerMatrix getBlock(int rowIndex, int colIndex, int numOfRows, int numOfCols) {
        IntegerMatrix M = new IntegerMatrix (numOfRows, numOfCols);
        getBlock(rowIndex, colIndex, numOfRows, numOfCols, M);
        return M;
    }

    /** store a block of size to <i>(numOfRows,numOfCol)</i> starting at <i>(rowIndex,colIndex)</i> in <i>M</i> */
    public IntegerMatrix getBlock(int rowIndex, int colIndex, int numOfRows, int numOfCols, IntegerMatrix M) {
        checkRowIndex(rowIndex);
        checkRowIndex(rowIndex + numOfRows - 1);
        checkColIndex(colIndex);
        checkColIndex(colIndex + numOfCols - 1);
        M.newSize(numOfRows, numOfCols);

        final int[] [] MRe = M.re;
        for (int i = rowIndex, I = 0; I < numOfRows; i++, I++)
            System.arraycopy(re[i], colIndex, MRe[I], 0, numOfCols);
        return M;
    }

    /** store <i>M</i> in this matrix beginning at <i>(rowIndex,colIndex)</i> */
    public void setBlock(int rowIndex, int colIndex, IntegerMatrix M) {
        setBlock(rowIndex, colIndex, M, M.numRows, M.numCols);
    }

    /** store a block from M of size <i>(numOfRows,numOfCols)</i> in this matrix beginning at (0,0) */
    public void setBlock(IntegerMatrix M, int numOfRows, int numOfCols) {
        setBlock(0, 0, M, numOfRows, numOfCols);
    }

    /** store a block from M of size <i>(numOfRows,numOfCols)</i> in this matrix beginning at <i>(rowIndex,colIndex)</i> */
    public void setBlock(int rowIndex, int colIndex, IntegerMatrix M, int numOfRows, int numOfCols) {
        checkRowIndex(rowIndex);
        checkRowIndex(rowIndex + numOfRows - 1);
        checkColIndex(colIndex);
        checkColIndex(colIndex + numOfCols - 1);
        M.checkRowIndex(numOfRows - 1);
        M.checkColIndex(numOfCols - 1);
        int[] [] MRe = M.re;
        for (int i = 0, j = rowIndex; i < numOfRows; i++, j++)
            System.arraycopy(MRe[i], 0, re[j], colIndex, numOfCols);
    }

    /** assign this matrix to the adjoint of M */
    public void assignAdjoint(IntegerMatrix M, int ithRow, int jthCol) {
        M.checkRowIndex(ithRow);
        M.checkColIndex(jthCol);
        newSize(M.numRows - 1, M.numCols - 1);
        for (int i = 0; i < ithRow; i++) {
            System.arraycopy(M.re[i], 0, re[i], 0, jthCol);
            System.arraycopy(M.re[i], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1);
        }
        for (int i = ithRow; i < numRows; i++) {
            System.arraycopy(M.re[i + 1], 0, re[i], 0, jthCol);
            System.arraycopy(M.re[i + 1], jthCol + 1, re[i], jthCol, M.numCols - jthCol - 1);
        }
        if ((ithRow + jthCol) % 2 == 1)
            assignTimes(-1);
    }

    public IntegerMatrix getAdjoint(int ithRow, int jthCol) {
        IntegerMatrix adjoint = new IntegerMatrix (numRows - 1, numCols - 1);
        adjoint.assignAdjoint(this, ithRow, jthCol);
        return adjoint;
    }

    /** stores <i>this</i>' diagonal in <i>diag</i> */
    public void getDiagonal( IntegerVector diag ) {
        int numDiag = Math.min(numRows, numCols);
        if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
	MatrixOperations.getDiagonal( re , diag.re ) ;
    }

    /** create and return a <b>de.jtem.blas.IntegerVector</b> containing <i>this</i>' diagonal */
    public IntegerVector getDiagonal() {
        IntegerVector diag = new IntegerVector (Math.min(numRows, numCols));
        getDiagonal(diag);
        return diag;
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( IntegerVector diag ) {
        int numDiag = Math.min(numRows, numCols);
        if (diag.size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
	MatrixOperations.assignDiagonal( re, diag.re );        
    }

    /** stores <i>diag</i> in <i>this</i>' diagonal while other entries remain unchanged */
    public void setDiagonal( int [] diag ) {
        int numDiag = Math.min(numRows, numCols);
        if (diag.length != numDiag) throw new IllegalArgumentException("vector has wrong length");
	MatrixOperations.assignDiagonal( re, diag );        
    }

    /** assign diagonal to @param diag while other entries remain unchanged */
    public void setDiagonal( final int diag ) {
        int numDiag = Math.min(numRows, numCols);
	MatrixOperations.assignDiagonal( re, diag );        
    }

    /** create and return an array <b>int[][]</b> containing <i>this</i>' entries */
    public int[][] toArray() {
        final int[][] res = new int[numRows][numCols];
	MatrixOperations.assign(re,res);
        return res;
    }

    /** store <i>this</i>' entries in <i>M</i> */
    public void toArray(int[] [] M) {
        checkShape(M);
        for (int i = 0; i < numRows; i++)
            System.arraycopy(re[i], 0, M[i], 0, numCols);
    }

    /** create and return a new matrix containing <i>this</i>' entries */
    public IntegerMatrix copy() {
        return new IntegerMatrix (this);
    }

    /** return a clone */
    public Object clone() { return copy(); }

    /** create and return an identity matrix of size <i>numrows*numrows</i> */
    public static IntegerMatrix id(int numrows) {
        final IntegerMatrix ID = new IntegerMatrix (numrows);
	ID.setDiagonal(1);
        return ID;
    }

    /** assign <i>this</i> to an identity matrix of size <i>numrows*numrows</i> */
    public void assignId(int numrows) {
        newSize(numrows, numrows);
        assignZero();
        for (int i = 0; i < numrows; i++) re[i] [i] = 1;
     }

    /** assign <i>this</i> to an identity matrix of size <i>numRows*numRows</i> */
    public void assignId() {
        assignId(numRows);
    }

    /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public IntegerMatrix plus( final IntegerMatrix b ) {
	checkShape(b);
	IntegerMatrix result = new IntegerMatrix ( numRows, numCols );
	MatrixOperations.plus( re, b.re, result.re);
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public IntegerMatrix minus( final IntegerMatrix b ) {
	checkShape(b);
	IntegerMatrix result = new IntegerMatrix ( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

     /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public RealMatrix plus( final RealMatrix b ) {
	checkShape(b);
	RealMatrix result = new RealMatrix ( numRows, numCols );
	MatrixOperations.plus( re, b.re, result.re);
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public RealMatrix minus( final RealMatrix b ) {
	checkShape(b);
	RealMatrix result = new RealMatrix ( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

     /** add <i>b</i> and <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix plus( final ComplexMatrix b ) {
	checkShape(b);
	ComplexMatrix result = new ComplexMatrix ( numRows, numCols );
	MatrixOperations.plus( re, b.re, result.re);
	return result;
    }

    /** subtract <i>b</i> from <i>this</i> and create and return a matrix containing the result */
    public ComplexMatrix minus( final ComplexMatrix b ) {
	checkShape(b);
	ComplexMatrix result = new ComplexMatrix ( numRows, numCols );
	MatrixOperations.minus( re, b.re, result.re );
	return result;
    }

   /** add <i>a</i> and <i>b</i> and store result in <i>this</i> matrix */
    public void assignPlus( final IntegerMatrix a, final IntegerMatrix b ) {
	b.checkShape(a);
	newSize( a );
	MatrixOperations.plus(a.re,b.re,re);
    }

    /** subtract <i>b</i> from <i>a</i> and store result in <i>this</i> matrix */
    public void assignMinus( final IntegerMatrix a, final IntegerMatrix b ) {
        a.checkShape(b);
	newSize( a );
	MatrixOperations.minus(a.re,b.re,re);
    }

    /** add <i>b</i> and <i>this</i> and store the result in <i>this</i> matrix */
    public void assignPlus( final IntegerMatrix b) {
        checkShape(b);
	MatrixOperations.plus(re,b.re,re);
    }

    /** subtract <i>b</i> from <i>this</i> and store the result in <i>this</i> matrix */
    public void assignMinus( final IntegerMatrix b) {
        checkShape(b);
	MatrixOperations.minus(re,b.re,re);
    }

    /* **************************************************************************************** */

    /* **************       matrix-scalar multiplication      ********************************* */

    /* **************************************************************************************** */

    /** multiply <i>this</i> and <i>a</i> and create and return a matrix containing the result */
    public IntegerMatrix times( final int a ) {
	IntegerMatrix result = new IntegerMatrix (numRows, numCols);
	MatrixOperations.times(re,a,result.re);
	return result;
    }

    /** multiply <i>this</i> and <i>a</i> and create and return a matrix containing the result */
    public RealMatrix times( final double a ) {
	RealMatrix result = new RealMatrix (numRows, numCols);
	MatrixOperations.times(re,a,result.re);
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times( final Field.Complex a ) {
	ComplexMatrix result = new ComplexMatrix(numRows, numCols);
	result.assignTimes( this, a );
	return result;
    }

    /** multiply <i>this</i> and <i>a</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final int a ) {
	MatrixOperations.times(re,a,re);
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final IntegerMatrix A, final int b ) {
        newSize(A.numRows, A.numCols);
        MatrixOperations.times(A.re, b, re);        
    }
    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignTimes( final int b,final IntegerMatrix A ) {
	assignTimes( A, b);
    }

     /** divides <i>this</i> by <i>a</i> and create and return a matrix containing the result */
    public IntegerMatrix divide( final int a ) {
	IntegerMatrix result = new IntegerMatrix (numRows, numCols);
	MatrixOperations.divide(re,a,result.re);
	return result;
    }

    /** divide <i>this</i> and <i>a</i> and create and return a matrix containing the result */
    public RealMatrix divide( final double a ) {
	RealMatrix result = new RealMatrix (numRows, numCols);
	MatrixOperations.times(re,(double)(1.0/a),result.re);
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix divide( final Field.Complex a ) {
	ComplexMatrix result = new ComplexMatrix(numRows, numCols);
	MatrixOperations.divide( re, a.getRe(), result.re );
	MatrixOperations.divide( re, a.getIm(), result.im );
	return result;
    }

    /** multiply <i>this</i> and <i>a</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final int a ) {
	MatrixOperations.divide(re,a,re);
    }

    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final IntegerMatrix A, final int b ) {
        newSize(A.numRows, A.numCols);
        MatrixOperations.divide(A.re, b, re);        
    }
    /** multiply <i>A</i> and <i>b</i> and store the result in <i>this</i> matrix */
    public void assignDivide( final int b,final IntegerMatrix A ) {
        newSize(A.numRows, A.numCols);
        MatrixOperations.divide( b, A.re, re);        
    }


    /* **************************************************************************************** */

    /* **************       matrix-matrix multiplication      ********************************* */

    /* **************************************************************************************** */

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public IntegerMatrix times(final IntegerMatrix A) {
	IntegerMatrix result = new IntegerMatrix (numRows, A.numCols);
	result.assignTimes(this, A);
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public RealMatrix times(final RealMatrix A) {
	RealMatrix result = new RealMatrix (numRows, A.numCols);
	result.assignTimes(this, A);
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and create and return a matrix containing the result */
    public ComplexMatrix times(final ComplexMatrix A) {
	ComplexMatrix result = new ComplexMatrix(numRows, A.numCols);
	result.assignTimes(this, A);
	return result;
    }

    /** multiply <i>this</i> and <i>A</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A) {
	assignTimes(this,A);
    }

    /** multiply <i>A</i> and <i>B</i> and store the result in <i>this</i> matrix */
    public void assignTimes(final IntegerMatrix A, final IntegerMatrix B) {
        if (A.numCols != B.numRows)
            throw new IllegalArgumentException("matrices have wrong sizes");
	int [][] tmpA = A.re;
	int [][] tmpB = B.re;
	if(A==this) tmpA = MatrixOperations.copy(re);
	if(B==this) tmpB = MatrixOperations.copy(re);	
        newSize(A.numRows, B.numCols);
	MatrixOperations.times(tmpA,tmpB, re);
    }

    /* *** end of matix-matrix multiplication *** */

    /** multiply <i>this</i> and <i>V</i> and create and return a vector containing the result */
    public IntegerVector times(final IntegerVector V) {
        IntegerVector res = new IntegerVector (numRows);
        res.assignTimes(this, V);
        return res;
    }

    /** mod all entries by <i>b</i> and create and return a matrix containing the result */
    public IntegerMatrix mod( int b ) {
	IntegerMatrix result = new IntegerMatrix (this);
	MatrixOperations.mod(result.re,b,result.re);
	return result;
    }

    /** mod all entries by <i>b</i> */
    public void assignMod(int b) {
	MatrixOperations.mod(re,b,re);
    }

    /** mod all entries of <i>M</i> by <i>b</i> and store the result in <i>this</i> matrix */
    public void assignMod(final IntegerMatrix M, int b) {
        assign(M);
	assignMod(b);
    }

    /** return <tt>true</tt> if <i>this</i> is symmetric */
    public boolean isSymmetric() {
        if (!isSquared())
            return false;
        for (int i = 0; i < numRows; i++) {
            int[] Rre = re[i];
            for (int j = 0; j < i; j++) {
                int dummy = Rre[j] - re[j] [i];
                if (dummy * dummy > EPSILON)
                    return false;
            }
        }
        return true;
    }

    /** return the sum of the squares of all entries */
    public double normSqr() {
	return MatrixOperations.normSqr(re);
    }

    /** return <tt>true</tt> if all entries are zero */
    public static double normSqr(final IntegerMatrix M) {
        return MatrixOperations.normSqr(M.re);
    }

    /** return the sum of the squares of all entries */
    public double norm() {
	return Math.sqrt(MatrixOperations.normSqr(re));
    }

    /** return <tt>true</tt> if all entries are zero */
    public static double norm(final IntegerMatrix M) {
        return Math.sqrt(MatrixOperations.normSqr(M.re));
    }

    /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero(final RealMatrix M) {
        return M.isZero();
    }

    /** return <tt>true</tt> if @param M equals the identity */
    public static boolean isId(final IntegerMatrix M) {
        return M.isId();
    }

    /** return <tt>true</tt> if this equals the identity */
    public boolean isId() {
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++) {
                final int[] ithRow = re[i];
                if (ithRow[j] != (i == j ? 1 : 0))
                    return false;
            }
        return true;
    }

    /** return <tt>true</tt> if <i>this</i> equals <i>m</i> */
    public boolean equals( final IntegerMatrix m ) {

	if (m.numRows != numRows || m.numCols != numCols) return false;
            final int[] [] mRe = m.re;
            for (int i = 0; i < numRows; i++) {
                final int[] rowRe = re[i];
                final int[] rowMRe = mRe[i];
                for (int j = 0; j < numCols; j++)
                    if (rowRe[j] != rowMRe[j]) return false;
            }
	    return true;
    }

    /** return <tt>true</tt> if <i>this</i> equals <i>o</i> */
    public boolean equals(Object o) {
        try {
	    return equals( (IntegerMatrix)o );
        }
        catch (ClassCastException ex) {
            return false;
        }
    }

   /** creates a <b>String</b> containing this' matrix' entries */
    public String toString() {
        StringBuffer sb = new StringBuffer(300);
        sb.append("(");
        for (int i = 0; i < numRows; i++) {
            sb.append('(');
            if (numCols > 0) {
                sb.append(re[i] [0]);
                for (int j = 1; j < numCols; j++) {
                    sb.append(", ");
                    sb.append(re[i] [j]);
                }
            }
            sb.append(')');
            if (i < numRows - 1) sb.append(",\n");
        }
        sb.append(')');
        return sb.toString();
    }

    /** transpose <i>M</i> and create and return a matrix containing the result */
    public static IntegerMatrix transpose( final IntegerMatrix M ) {
        IntegerMatrix res = new IntegerMatrix (M.numCols, M.numRows);
        MatrixOperations.transpose( M.re, res.re );
        return res;
    }

    /** transpose <i>this</i> matrix create and return a matrix containing the result */
    public IntegerMatrix transpose() {
        IntegerMatrix res = new IntegerMatrix (numCols, numRows);
        MatrixOperations.transpose( re, res.re );
        return res;
    }

    /** transpose <i>M</i> and store the result in <i>this</i> matrix */
    public void assignTranspose( final IntegerMatrix M ) {
      	
	IntegerMatrix tmpM = (  M.isSquared() || M != this ) ? M : new IntegerMatrix ( this );
	newSize( tmpM.numCols, tmpM.numRows );

	MatrixOperations.transpose( tmpM.re, re );
    }

    /** transpose <i>this</i> matrix */
    public void assignTranspose() {
        assignTranspose(this);
    }

    public int determinant() {
        checkSquared();
        switch (numRows) {
            case 0:
                return 0;
            case 1:
                return re[0] [0];
            case 2:
                return re[0] [0] * re[1] [1] - re[0] [1] * re[1] [0];
            case 3:
                return re[0] [0] * re[1] [1] * re[2] [2] + re[0] [1] * re[1] [2] * re[2] [0] + re[0] [2] * re[1] [0] *
                    re[2] [1] - re[2] [0] * re[1] [1] * re[0] [2] - re[2] [1] * re[1] [2] * re[0] [0] - re[2] [2] *
                    re[1] [0] * re[0] [1];
            default: {
                    RealMatrix dummy = new RealMatrix (numRows);
		    dummy.assign(this);
                    return (int)dummy.determinant();
                }
        }
    }

    /**
     * Computes and returns the inverse of <em>this</em> matrix. <p />
     * @return	A <code>IntegerMatrix</code> containing the result.
     */
    public IntegerMatrix invert() {
        IntegerMatrix inv = new IntegerMatrix (this);
        inv.assignInvert();
        return inv;
    }

    /**
     * Computes and returns the inverse of matrix <em>M</em>. <p />
     * @param	M   A <code>IntegerMatrix</code> to invert.
     * @return	A <code>IntegerMatrix</code> containing the result.
     */
    public static IntegerMatrix invert(final IntegerMatrix M) {
        IntegerMatrix inv = new IntegerMatrix (M);
        inv.assignInvert();
        return inv;
    }

    /**
     * Computes the inverse of matrix <em>M</em> and stores the result in <em>this</em> matrix. <p />
     * @param	M   A <code>IntegerMatrix</code> to invert.
     */
    public void assignInvert(final IntegerMatrix M) {
        assign(M);
        assignInvert();
    }

    /**
     * Inverts <em>this</em> matrix.
     */
    public void assignInvert() {
        checkSquared();
        RealMatrix dummy = new RealMatrix (numRows);
        dummy.assign(this);
        dummy.assignInvert();
        IntegerMatrix inv = new IntegerMatrix (numRows);
	MatrixOperations.round( dummy.re, inv.re );
        for (int i = 0; i < numRows; i++) {
            final int[] invRow = inv.re[i];
            final double[] dummyRow = dummy.re[i];           
            // we check the result
            for (int j = 0; j < numRows; j++) {
                int sum = 0;
                for (int k = 0; k < numRows; k++)
                    sum += invRow[k] * re[k] [j];
                if (sum != (j == i ? 1 : 0))
                    throw new RuntimeException("matrix not invertable");
            }
        }
        assign(inv);        
    }

    /** assign all entries to zero */
    public void assignZero() {
    	MatrixOperations.assign(re, 0);        
    }

    public static java.util.Random iRandom = new java.util.Random();

    /** assign all entries randomly */
    public void assignRandom() {
        assignRandom(1024);        
    }

    public void assignRandom(int n) {
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++)
                re[i] [j] = iRandom.nextInt(n);        
    }

    /** returns round of @param A */
    public static IntegerMatrix round(final RealMatrix A) {
        IntegerMatrix result = new IntegerMatrix (A.numRows, A.numCols);
        MatrixOperations.round(A.re, result.re);
        return result;
    }

    /** assign <em>this</em> to round of @param A */
    public void assignRound(final RealMatrix A) {
        newSize(A.numRows, A.numCols);
        MatrixOperations.round(A.re, re);
    }

    /** returns floor of @param A */
    public static IntegerMatrix floor( final RealMatrix A ) {
        IntegerMatrix result = new IntegerMatrix (A.numRows, A.numCols);
        MatrixOperations.floor(A.re, result.re);
        return result;
    }

    /** assign <em>this</em> to floor of @param A */
    public void assignFloor( final RealMatrix A ) {
        newSize(A.numRows, A.numCols);
        MatrixOperations.floor(A.re, re);
    }

    /** creates int matrix and assigns it with neg of <code>v</code> */
    public final static IntegerMatrix neg( final IntegerMatrix im ) {
	IntegerMatrix result = new IntegerMatrix ( im.numRows ,im.numCols );
	MatrixOperations.neg( im.re, result.re );
	return result;
    }

    /** assign neg of <code>v</code> */
    public final void assignNeg( final IntegerMatrix im ) {
	newSize( im.numRows,im.numCols,false );
	MatrixOperations.neg( im.re, re );
    }

    /** assign this to neg of <code>this</code> */
    public final void assignNeg() {
	MatrixOperations.neg( re, re );
    }

}
