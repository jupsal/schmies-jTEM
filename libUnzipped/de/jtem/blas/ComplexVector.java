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

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.Field;
import de.jtem.numericalMethods.algebra.linear.MatrixOperations;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;

/**
 * This class represents general complex vectors.
 *
 * Following the number package's philosophy, operations on complex vectors
 * generally occur three-fold: For instances a,b,c of this class,
 * <table>
 * <tr><td valign="up"><tt>a.plus(b)</tt></td><td>creates and returns
 * a new instance containing the sum <i>a+b</i></td></tr>
 * <tr><td valign="up"><tt>a.assignPlus(b)</tt></td><td>stores the
 * sum <i>a+b</i> in <i>a</i></td></tr>
 * <tr><td valign="up"><tt>c.assignPlus(a,b)</tt></td><td>stores the
 * sum <i>a+b</i> in c</td></tr>
 * </table>
 * <p>
 *
 * @author	marina
 */
public final class ComplexVector extends AbstractVector {

    private static final long serialVersionUID = 1L;

    public double[] re;
    public double[] im;

    /** create a zero vector */
    public ComplexVector() {
	re=new double[0];
	im=new double[0];
	length=0;
    }

    /** create a vector with all entries equal to zero */
    public ComplexVector(int size){
	re=new double[size];
	im=new double[size];
	length=size;
    }

    /** create a vector and initialize all entries with <i>initValue</i> */
    public ComplexVector( int size, final Field.Complex initValue ){
	this( size );
	assign( initValue );
    }

    /** create a vector and initialize all entries with <i>initValue</i> */
    public ComplexVector( final int size, final double initValue){
	this( size );
	assign( initValue );
    }

    /** create a vector and initialize all entries with the <i>initRe,initIm</i> */
    public ComplexVector( final int size, final double initRe, final double initIm )    {
	this( size );
	assign( initRe, initIm );
    }

    /** create copy of vector V */
    public ComplexVector( final ComplexVector V )    {
	re = (double[])V.re.clone();
	im = (double[])V.im.clone();
	length = V.length;
    }

    /** create copy of vector V */
    public ComplexVector( final RealVector V )    {
	re = (double[])V.re.clone();
	length = V.length;
	im = new double[length];
    }

    /** create copy of vector V */
    public ComplexVector( final IntegerVector V )    {
	this( V.length );
	VectorOperations.assign( V.re(), re );
    }

    /** create a vector with initialized with <i>b</i> */
    public ComplexVector( final Complex[] V)    {
	this( V.length );
	for(int i=0; i<length; i++)	    {
	    re[i]=V[i].re;
	    im[i]=V[i].im;
	}
    }

     /** create a vector with initialized with <i>b</i> */
    public ComplexVector( final ComplexConstant[] V)    {
	this( V.length );
	for(int i=0; i<length; i++)	    {
	    re[i]=V[i].getRe();
	    im[i]=V[i].getIm();
	}
    }

    /** create a vector with initialized with <i>b</i> */
    public ComplexVector( final Field.Complex[] V)    {
	this( V.length );
	for(int i=0; i<length; i++)	    {
	    re[i]=V[i].getRe();
	    im[i]=V[i].getIm();
	}
    }

   /** create a vector with initialized with <i>Re</i> and <i>Im</i> if possible */
    public ComplexVector( final double[] Re, final double[] Im){
        if (Re.length==0 && Im.length==0) newSize(0);
	else {
	    re = (double[])Re.clone();
	    im = (double[])Im.clone();
	    length = Re.length;
	}
    }

    /** assigns <i>this</i> with <i>reV</i> and <i>imV</i> if possible */
    public void assign( final RealVector reV, final RealVector imV )    {
	assign( reV.re, imV.re);
    }

    /** assigns <i>this</i> with <i>re</i> and <i>im</i> if possible */
    public void assign( final double [] re, final double [] im ){
	VectorOperations.checkShape( re, im);
	newSize( re.length );
	VectorOperations.assign( re, this.re );
	VectorOperations.assign( im, this.im );

    }

    void newSize( final int newSize, final boolean copyValues ) {

	if( newSize != length ) {
	    double [] newRe = new double[newSize];
	    double [] newIm = new double[newSize];

	    int minSize = newSize < length ? newSize : length;

	    if( copyValues ) {
		System.arraycopy( re, 0, newRe, 0, minSize );
		System.arraycopy( im, 0, newIm, 0, minSize );
	    }

	    length = newSize;

	    re = newRe;
	    im = newIm;
	}
    }

    /** Returns real part of entry at <i>index</i> */
    public double getRe( final int index ) {
	return re[index];
    }

    /** Sets real part of entry at <i>index</i> */
    public void setRe( final int index, final double aValue ) {
	re[index] = aValue;
    }

    /** Return imaginary part of entry at <i>index</i> */
    public double getIm( final int index ) {
	return im[index];
    }


    /** Sets imaginary part of entry at <i>index</i> */
    public void setIm( final int index, final double aValue ) {
	im[index] = aValue;
    }


    /** create and return new <b>Field.Complex</b> containing entry at <i>index</i> */
    public Complex get( final int index )     {
	return new Complex( re[index], im[index] );
    }

    /** stores <b>Field.Complex</b> containing entry at <i>index</i> */
    public void get( final int index , Complex c)     {
      c.re = re[index];
      c.im = im[index];
    }

    /** store <i>c</i> in entry <i>index</i> */
    public void set( final int index, final ComplexConstant c  ) {
	re[index]=c.getRe();
	im[index]=c.getIm();
    }

    /** store <i>c</i> in entry <i>index</i> */
  public void set( final int index, final Field.Complex c  ) {
      re[index]=c.getRe();
      im[index]=c.getIm();
  }

  /** store <i>c</i> in entry <i>index</i> */
    public void set( final int index, final Complex c  ) {
        re[index]=c.re;
        im[index]=c.im;
    }


    /** store <i>(aRe,aIm)</i> in entry <i>index</i> */
    public void set( final int index, final double aRe, final double aIm  ) {
	re[index]=aRe; im[index]=aIm;
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final ComplexVector V ) {
	newSize(V.length);

	System.arraycopy( V.re, 0, re, 0, length );
	System.arraycopy( V.im, 0, im, 0, length );
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final Complex [] V ) {
	newSize(V.length);
	for(int i=0; i<length; i++){
		re[i]=V[i].re;
		im[i]=V[i].im;
	    }
    }

     /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final ComplexConstant [] V ) {
	newSize(V.length);
	for(int i=0; i<length; i++){
		re[i]=V[i].getRe();
		im[i]=V[i].getIm();
	    }
    }

      /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final Field.Complex [] V ) {
	newSize(V.length);
	for(int i=0; i<length; i++){
		re[i]=V[i].getRe();
		im[i]=V[i].getIm();
	    }
    }

  /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final int [] V ) {
	newSize(V.length);
	VectorOperations.assign( V, re );
	VectorOperations.assign( im, 0 );
    }

    /** assign <i>this<i/> vector to <i>V</i> */
    public void assign( final double [] V ) {
	newSize(V.length);
	VectorOperations.assign( V, re );
	VectorOperations.assign( im, 0 );
    }

    /** assign all entries to <i>(Re,Im)</i> */
     public void assign( final double Re, final double Im )  {
	VectorOperations.assign(re,Re);
	VectorOperations.assign(im,Im);
    }

    /** assign all entries to <i>b</i> */
    public void assign( final Field.Complex b )     {
        assign( b.getRe(), b.getIm() );
    }

    /** assign all entries to (<i>r</i>,0) */
    public void assign( final double r ) {
        VectorOperations.assign(re,r);
        VectorOperations.assign(im,0);
    }

    /** assign vector to <i>V</i> */
    public void assign( final RealVector V ) {
	newSize( V.length );
	System.arraycopy( V.re, 0, re, 0, length );
	VectorOperations.assign(im,0.0);
    }

      /** assign vector to <i>V</i> */
    public void assign( final IntegerVector V ) {
	newSize( V.length );
	VectorOperations.assign(V.re,re);
	VectorOperations.assign(im,0.0);
    }

    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal( final ComplexMatrix diag ){
	diag.setDiagonal( this );
    }

    /** stores <i>this</i> in <i>diag</i>' diagonal while other entries remain unchanged */
    public void assignDiagonal( Complex[][] diag ){
        int numDiag = Math.min(diag.length, diag[0].length);
        if (size() != numDiag) throw new IllegalArgumentException("vector has wrong length");
	for (int i =0; i<numDiag; i++ ){
	    diag[i][i].re = re[i];
	    diag[i][i].im = im[i];
	}
    }

    /** stores <i>this</i> in i. row of <i>row</i> while other entries remain unchanged.
     * i goes from 0 to length-1
     */
    public void assignRow( ComplexMatrix row, int i ){
	row.setRow( i, this );
    }

    /** stores <i>this</i> in i. row of <i>row</i>  while other entries remain unchanged
     * i goes from 0 to length-1
     */
    public void assignRow( Complex[][] row, int i ){
	if( row.length != length ) throw new IllegalArgumentException("sizes does not match");
       	if( row[0].length <= i ) throw new IllegalArgumentException("wrong second parameter");
	for(int j=0; j<row.length; j++ ){
	    row[i][j].re = re[j];
	    row[i][j].im = im[j];
	}
    }

    /** stores <i>this</i> in i. column of <i>col</i>  while other entries remain unchanged.
     * i goes from 0 to length-1
     */
    public void assignCol( ComplexMatrix col, int i ){
	col.setCol( i, this );
    }

    /** stores <i>this</i> in i. column of <i>col</i>  while other entries remain unchanged.
     * i goes from 0 to length-1
     */
    public void assignCol( Complex[][] col, int i ){
       	if( col.length <= i ) throw new IllegalArgumentException("wrong second parameter");
       	if( col[0].length != length ) throw new IllegalArgumentException("sizes does not match");
	for(int j=0; j<col[0].length; j++ ){
	    col[j][i].re = re[j];
	    col[j][i].im = im[j];
	}
    }

    /** The following four methods return/set REFERENCES only ! */

    /** return a REFERENCE to <i>this'</i> real part */
    public double[] re()     {
	return re;
    }

    /** return a REFERENCE to <i>this'</i> imaginary part */
    public double[] im() {
      return im;
    }

    /** set <i>this'</i> real part by REFERENCING <i>V</i>  */
    public void re( final double[] V ) {
      if (V.length != length)
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      re = V;
    }

    /** set <i>this'</i> imaginary part by REFERENCING <i>V</i> */
    public void im( final double[] V ) {
      if (V.length != length)
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      im = V;
    }

    /** set <i>this</i>' real part to <i>V</i>, leave imaginary part unchanged */
    public void setRe( final double[] V ) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      VectorOperations.assign(V,re);
    }

    /** set <i>this</i>' imaginary part to <i>V</i>, leave real part unchanged */
    public void setIm( final double[] V ) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      System.arraycopy( V, 0, im, 0, length );

    }

     /** set <i>V</i> to <i>this</i>' real part  */
    public void getRe( final double[] V ) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      VectorOperations.assign(re,V);
    }

    /** set <i>V</i> to <i>this</i>' imaginary part  */
    public void getIm( final double[] V ) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      System.arraycopy( im, 0, V, 0, length );
    }

    /** create and return new <b>de.jtem.blas.RealVector</b> containing <i>this</i>' real part */
    public RealVector getRe() {
      return new RealVector ( re );
    }

    /** create and return new <b>de.jtem.blas.RealVector</b> containing <i>this</i>' imaginary part */
    public RealVector getIm() {
      return new RealVector ( im );
    }

    /** store <i>this</i>' real part in <i>V</i> if possible */
    public void getRe( final RealVector V ) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      V.assign(re);
    }

    /** store <i>this</i>' real part in <i>V</i> if possible */
    public void getIm(RealVector V) {
      if( V.length != length )
	throw new IllegalArgumentException("array does not have the same length as 'this'");
      V.assign(im);
    }

    /** set <i>this</i>' real part to <i>V</i>, leave imaginary part unchanged */
    public void setRe( final RealVector V ) {
	setRe(V.re);
    }

    /** set <i>this</i>' imaginary part to <i>V</i>, leave real part unchanged */
    public void setIm( final RealVector V ) {
	setIm(V.re);
    }


    /** create and return a block from <i>anIndex</i> to the bottom */
    public ComplexVector getBlock( final int anIndex ) {
	return getBlock( anIndex, length-anIndex );
    }

    /** store a block from an <i>anIndex</i> in <i>V</i> */
    public void getBlock( final int anIndex, final ComplexVector V ) {
	getBlock( anIndex, length-anIndex, V);
    }

    /** create and return a block of size <i>aSize</i> beginning at <i>anIndex</i> */
    public ComplexVector getBlock( final int anIndex, final int aSize ) {
	ComplexVector V = new ComplexVector(aSize);
	getBlock(anIndex,aSize,V);
	return V;
    }

    /** store a block of size <i>aSize</i> beginning at <i>anIndex</i> in <i>V</i> */
    public void getBlock( final int anIndex, final int aSize, final ComplexVector V ) {

	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );

	V.newSize( aSize );

	System.arraycopy( re, anIndex, V.re, 0, aSize );
	System.arraycopy( im, anIndex, V.im, 0, aSize );
    }

    /** store <i>V</i> as a block in <i>this</i> vector, beginning at <i>anIndex</i> */
    public void setBlock( final int anIndex, final ComplexVector V ) {
	setBlock( anIndex, V, V.length );
    }

    /** store a block of size <i>aSize</i> beginning at <i>V</i>'s top
        in <i>this</i> vector, start at <i>anIndex</i> */
    public void setBlock( final int anIndex, final ComplexVector V, final int aSize ) {

	checkIndex( anIndex );
	checkIndex( anIndex + aSize - 1 );

	V.checkIndex( aSize - 1 );

	System.arraycopy( V.re, 0, re, anIndex, aSize );
	System.arraycopy( V.im, 0, im, anIndex, aSize );
    }


    /** create and return an array <b>Field.Complex[]</b> containing <i>this</i>' entries */
    public Complex[] toArray()
    {
	Complex[] res=new Complex[length];
	for(int i=0; i<length; i++){
	    res[i]=new Complex(re[i],im[i]);
	}
	return res;
    }

    /** store <i>this</i>' entries in <i>res</i>*/
    public void toArray( final Complex[] res )    {
      checkLength(res.length);
      for(int i=0; i<length; i++)
	  if( res[i] != null )
	      res[i].assign(re[i],im[i]);
	  else
	      res[i] = new Complex(re[i],im[i]);
    }

   /** create and return a new vector containing <i>this</i>' entries */
    public ComplexVector copy(){
	return new ComplexVector(this);
    }

    /** return a clone */
    public Object clone() { return copy(); }

    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public ComplexVector plus( final ComplexVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignPlus(this, V);
	return result;
    }

    /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public ComplexVector plus( final RealVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignPlus(this, V);
	return result;
    }

     /** add <i>V</i> and <i>this</i> and create and return a vector containing the result */
    public ComplexVector plus( final IntegerVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignPlus(this, V);
	return result;
    }

   /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public ComplexVector minus( final ComplexVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignMinus(this, V);
	return result;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public ComplexVector minus( final RealVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignMinus(this, V);
	return result;
    }

    /** subtract <i>V</i> from <i>this</i> and create and return a vector containing the result */
    public ComplexVector minus( final IntegerVector V ){
	ComplexVector result = new ComplexVector( length );
	result.assignMinus(this, V);
	return result;
    }

    /** add <i>V</i> and <i>this</i> and store result in <i>this</i> vector */
    public void assignPlus( final ComplexVector V ){
	checkShape(V);
	VectorOperations.plus(re,V.re,re);
	VectorOperations.plus(im,V.im,im);
    }

     /** add <i>V</i> and <i>this</i> and store result in <i>this</i> vector.
      * imag. entry of <i>this</i> will be left untouched.
     */
    public void assignPlus( final RealVector V ){
	checkShape(V);
	VectorOperations.plus(re,V.re,re);
    }

      /** add <i>V</i> and <i>this</i> and store result in <i>this</i> vector.
      * imag. entry of <i>this</i> will be left untouched.
     */
    public void assignPlus( final IntegerVector V ){
	checkShape(V);
	VectorOperations.plus(re,V.re,re);
    }

  /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final ComplexVector V ){
	checkShape(V);
	VectorOperations.minus(re,V.re,re);
	VectorOperations.minus(im,V.im,im);
    }

  /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V ){
	checkShape(V);
	VectorOperations.minus(re,V.re,re);
    }

  /** subtract <i>V</i> from <i>this</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V ){
	checkShape(V);
	VectorOperations.minus(re,V.re,re);
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final ComplexVector V, final ComplexVector W )   {
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.plus( W.im, V.im, im );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final RealVector V, final RealVector W )   {
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.assign( im,0 );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final IntegerVector V, final IntegerVector W )   {
	V.checkShape( W );
	newSize( V.length );
	int[] res = new int[length];
	VectorOperations.plus( W.re, V.re, res );
	VectorOperations.assign( res, re );
	VectorOperations.assign( im,0 );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final IntegerVector V, final ComplexVector W )   {
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.assign( W.im, im );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final IntegerVector V, final RealVector W )   {
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.assign( im,0 );
    }

     /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final RealVector V, final IntegerVector W )   {
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.assign( im,0 );
    }

   /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final ComplexVector V, final IntegerVector W )    {
	assignPlus( W, V );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public void assignPlus( final RealVector V, final ComplexVector W ){
	assignPlus( W, V );
    }

    /** add <i>W</i> and <i>V</i> and store result in <i>this</i> vector */
    public final void assignPlus( final ComplexVector V, final RealVector W ){
	V.checkShape(W);
	newSize( V.length );

	VectorOperations.plus( W.re, V.re, re );
	VectorOperations.assign(V.im, im );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final ComplexVector V, final ComplexVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.minus( V.im, W.im, im );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V, final RealVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.assign( im,0);
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final IntegerVector W ){
	V.checkShape( W );
	newSize( V.length );
	int[] res = new int[length];
	VectorOperations.minus( V.re, W.re, res );
	VectorOperations.assign(res,re);
	VectorOperations.assign( im,0);
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V, final ComplexVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.neg( W.im, im );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final ComplexVector V, final RealVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.assign( V.im, im );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final ComplexVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.neg( W.im, im );
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final ComplexVector V, final IntegerVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.assign( V.im, im );
    }


    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final IntegerVector V, final RealVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.assign( im,0);
    }

    /** subtract <i>W</i> from <i>V</i> and store result in <i>this</i> vector */
    public void assignMinus( final RealVector V, final IntegerVector W ){
	V.checkShape( W );
	newSize( V.length );
	VectorOperations.minus( V.re, W.re, re );
	VectorOperations.assign( im,0);
    }

    /** multiplies real entries of <i>V</i> with <i>W</i> and returns it back .
     */
    public static double dotRe( final ComplexVector V, final ComplexVector W ){
	V.checkShape(W);
	return VectorOperations.times(V.re,W.re);
    }

     /** multiplies real entries of <i>V</i> with <i>this</i> and returns it back .
     */
    public double dotRe( final ComplexVector V ){
	checkShape(V);
	return VectorOperations.times(V.re,re);
    }

    /** multiplies imag. entries of <i>V</i> with <i>W</i> and returns it.
      */
    public static double dotIm( final ComplexVector V, final ComplexVector W ){
	V.checkShape(W);
	return VectorOperations.times(V.im,W.im);
    }


     /** multiplies imag. entries of <i>V</i> with <i>this</i> and returns it.
      */
    public double dotIm( final ComplexVector V){
	checkShape(V);
	return VectorOperations.times(V.im,im);
    }

    /** multiplies entries of <i>V</i> with <i>W</i> pairwise and stores
     *	the sum of these products.
     */
    public static void dotBilinear(  final ComplexVector V, final RealVector W,  final Complex VW ) {
	V.checkShape(W);

	VW.assign( VectorOperations.times(V.re,W.re),
		   VectorOperations.times(V.im,W.re));
    }

    /** multiplies entries of <i>V</i> with <i>W</i> pairwise and stores
     *	the sum of these products.
     */
    public static void dotBilinear( final ComplexVector V, final IntegerVector W, final Complex VW ) {
	V.checkShape(W);
	VW.assign( VectorOperations.times( V.re, W.re ),
		   VectorOperations.times( V.im, W.re ) );
    }

     /** multiplies entries of <i>V</i> with <i>W</i> pairwise and stores
     *	the sum of these products.
     */
    public static void dotBilinear(  final RealVector W, final ComplexVector V,  final Complex VW ) {
	V.checkShape(W);

	VW.assign( VectorOperations.times(W.re,V.re),
		   VectorOperations.times(W.re,V.im));
    }

    /** multiplies entries of <i>V</i> with <i>W</i> pairwise and stores
     *	the sum of these products.
     */
    public static void dotBilinear(  final IntegerVector W, final ComplexVector V, final Complex VW ) {
	V.checkShape(W);
	VW.assign( VectorOperations.times( W.re, V.re ),
		   VectorOperations.times( W.re, V.im ) );
    }

  /**@deprecated */
   public static void dot( final ComplexVector V, final RealVector W, final Complex VW ) {
      dotBilinear( V, W, VW);
    }

  /**@deprecated */
    public static void dot( final ComplexVector V, final IntegerVector W, final Complex VW ) {
      dotBilinear( V, W, VW);
    }

  /**@deprecated */
   public static void dot(  final RealVector W, final ComplexVector V, final Complex VW ) {
      dotBilinear( V, W, VW);
    }

  /**@deprecated */
    public static void dot(  final IntegerVector W, final ComplexVector V, final Complex VW ) {
      dotBilinear( V, W, VW);
    }


   /** multiplies entries of <i>V</i> with <i>W</i> pairwise and returns
     *	the sum of these products.
     */
    public static Complex dotBilinear( final ComplexVector V,  final ComplexVector W ) {

	V.checkShape( W );
	final double [] VRe = V.re;
	final double [] VIm = V.im;
	final double [] WRe = W.re;
	final double [] WIm = W.im;
	final int dim = V.length;
	double rr=0, ii=0;
	for(int i=0; i<dim; i++) {
	    final double r1=VRe[i];
	    final double r2=WRe[i];
	    final double i1=VIm[i];
	    final double i2=WIm[i];
	    rr += r1*r2 - i1*i2;
	    ii += r1*i2 + r2*i1;
	}
	return new Complex( rr, ii );
    }



    /** multiplies entries of <i>V</i> with <i>W</i> pairwise and stores
	the sum of these products.
     */
    public static void dotBilinear( final ComplexVector V,
			    final ComplexVector W, final Complex VW ) {

	V.checkShape( W );

	final double [] VRe = V.re;
	final double [] VIm = V.im;
	final double [] WRe = W.re;
	final double [] WIm = W.im;
	final int dim = V.length;

	double rr=0, ii=0;

	for(int i=0; i<dim; i++) {
	    final double r1=VRe[i];
	    final double r2=WRe[i];
	    final double i1=VIm[i];
	    final double i2=WIm[i];

	    rr += r1*r2 - i1*i2;
	    ii += r1*i2 + r2*i1;
	}

	VW.assign(rr,ii);
    }

  /**@deprecated */
    public static void dot( final ComplexVector V, final ComplexVector W,
			                              final Complex VW ) {
      dotBilinear( V, W, VW );
    }


    /** multiplies entries of <i>V</i> with <i>this</i> pairwise and returns
     *	the sum of these products.
     */
    public  Complex dotBilinear( final ComplexVector V ) {
	checkShape( V );
	final double [] VRe = V.re;
	final double [] VIm = V.im;
	final int dim = V.length;
	double rr=0, ii=0;
	for(int i=0; i<dim; i++) {
	    final double r1=VRe[i];
	    final double r2=re[i];
	    final double i1=VIm[i];
	    final double i2=im[i];
	    rr += r1*r2 - i1*i2;
	    ii += r1*i2 + r2*i1;
	}
	return new Complex( rr, ii );
    }

  /**@deprecated */
   public  Complex dot( final ComplexVector V ) {
     return dotBilinear( V );
   }



    /** multiply all entries with b and create and return a vector containing the result */
     public ComplexVector times( final Complex b )    {
	ComplexVector result = new ComplexVector(this);
	result.assignTimes(b);
	return result;
     }

    /** multiply all entries with r and create and return a vector containing the result */
    public ComplexVector times( final double r ){
	ComplexVector result = new ComplexVector(this);
	result.assignTimes(r);
	return result;
    }

    /** multiply all entries with b */
    public void assignTimes( final Field.Complex b )    {
	VectorOperations.times( re, im, b.getRe(), b.getIm(), re, im);
    }

    /** multiply all entries with r */
    public void assignTimes( final double r ) {
	VectorOperations.times(re,r,re);
	VectorOperations.times(im,r,im);
    }

     /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final double c ) {
	newSize( V.length );

	VectorOperations.times( V.re, c, re );
	VectorOperations.assign( im, 0 );
    }

     /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final double c,final IntegerVector V ) {
	assignTimes( V,c );
    }

     /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealVector V, final double c ) {
	newSize( V.length );

	VectorOperations.times( V.re, c, re );
	VectorOperations.assign( im, 0 );
    }

     /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final double c , final RealVector V ) {
	assignTimes( V,c );
    }

   /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealVector V, final Field.Complex c ) {
	newSize( V.length );
	VectorOperations.times( V.re, c.getRe(), re );
	VectorOperations.times( V.re, c.getIm(), im );
    }

    /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final Field.Complex c ,final RealVector V ) {
	assignTimes( V,c );
    }

  /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes(  final Field.Complex c, final IntegerVector V) {
	assignTimes( V,c );
    }

  /** multiply <i>V</i> with <i>c</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final Field.Complex c ) {
	newSize( V.length );

	VectorOperations.times( V.re, c.getRe(), re );
	VectorOperations.times( V.re, c.getIm(), im );
    }

  /** multiply all entries of <i>V</i> with <i>b</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexVector V, final Field.Complex b ){
	newSize( V.length );
	VectorOperations.times(V.re,V.im, b.getRe(), b.getIm(), re, im );
    }

    /** multiply all entries of <i>V</i> with <i>b</i> and store the result in <i>this</i> vector */
    public void assignTimes( final Field.Complex b, final ComplexVector V )    {
	assignTimes(V, b);
    }

   /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexVector V, final double r ) {
	newSize( V.length );
	VectorOperations.times(V.re, r, re );
	VectorOperations.times(V.im, r, im );
    }

    /** multiply all entries of <i>V</i> with <i>r</i> and store the result in <i>this</i> vector */
    public void assignTimes( final double r, final ComplexVector V ) {
	assignTimes(V, r);
    }

    /** divide all entries by b and create and return a vector containing the result */
    public ComplexVector divide( final Field.Complex b )    {
	ComplexVector result = new ComplexVector(this);
	result.assignDivide(b);
	return result;
    }

    /** divide all entries by b and create and return a vector containing the result */
    public ComplexVector divide( final double r ){
	ComplexVector result = new ComplexVector(this);
	result.assignDivide(r);
	return result;
    }

    /** divide all entries by b */
    public void assignDivide( final Field.Complex b ){
	assignDivide( this, b);
    }

    /** divide all entries by r */
    public void assignDivide( final double r ){
	assignDivide( this, r );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final ComplexVector V, final Field.Complex b ){
	newSize( V.length );
	VectorOperations.divide( V.re, V.im, b.getRe(), b.getIm(), re, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final ComplexVector V, final double r ){
	newSize( V );
	VectorOperations.divide( V.re , r, re );
	VectorOperations.divide( V.im , r, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final RealVector V, final Field.Complex b ){
	newSize( V.length );
	VectorOperations.divide( V.re, V.re, b.getRe(), b.getIm(), re, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final RealVector V, final double r ){
	newSize( V );
	VectorOperations.divide( V.re , r, re );
	VectorOperations.assign( im, 0. );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final IntegerVector V, final Field.Complex b ){
	newSize( V.length );
	VectorOperations.divide( V.re, b.getRe(), re );
	VectorOperations.divide( V.re, b.getIm(), im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final IntegerVector V, final double r ){
	newSize( V );
	VectorOperations.divide( V.re , r, re );
	VectorOperations.assign( im, 0. );
    }
    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide(  final Field.Complex b, final ComplexVector V ){
	newSize( V.length );
	VectorOperations.divide(b.getRe(), V.re, re );
	VectorOperations.divide(b.getIm(), V.im, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final double r,final ComplexVector V ){
	newSize( V );
	VectorOperations.divide( r, V.re , V.im, re, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide(  final Field.Complex b,final RealVector V ){
	newSize( V.length );
	VectorOperations.divide(b.getRe(), V.re, re );
	VectorOperations.divide(b.getIm(), V.re, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide(  final double r ,final RealVector V){
	newSize( V );
	VectorOperations.divide(  r, V.re , re );
	VectorOperations.assign( im, 0. );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide( final Field.Complex b ,final IntegerVector V ){
	newSize( V.length );
	VectorOperations.divide( b.getRe(), V.re,  re );
	VectorOperations.divide( b.getIm(), V.re, im );
    }

    /** divide all entries of <i>V</i> by <i>b</i> and store the result in <i>this</i> vector */
    public void assignDivide(  final double r,final IntegerVector V ){
	newSize( V );
	VectorOperations.divide(  r, V.re , re );
	VectorOperations.assign( im, 0. );
    }

    /** multiply  <i>this</i> and <i>M</i> and create and return a vector containing the result */
    public ComplexVector times( final ComplexMatrix M ) {
	ComplexVector result = new ComplexVector( M.numCols );
	result.assignTimes( this, M );
	return result;
    }

    /** multiply  <i>this</i> and <i>M</i> and create and return a vector containing the result */
    public ComplexVector times( final RealMatrix M ) {
	ComplexVector result = new ComplexVector( M.numCols );
	result.assignTimes( this, M );
	return result;
    }

    /** multiply  <i>this</i> and <i>M</i> and create and return a vector containing the result */
    public ComplexVector times( final IntegerMatrix M ) {
	ComplexVector result = new ComplexVector( M.numCols );
	result.assignTimes( this, M );
	return result;
    }

    /** multiply  <i>M</i> and <i>V</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexMatrix M, final IntegerVector V ) {
	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );

	MatrixOperations.times( M.re, V.re, re );
	MatrixOperations.times( M.im, V.re, im );
    }

    /** multiply  <i>M</i> and <i>V</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexMatrix M, final RealVector V ) {
	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );

	MatrixOperations.times( M.re, V.re, re );
	MatrixOperations.times( M.im, V.re, im );
    }

    /** multiply  <i>M</i> and <i>V</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexMatrix M, final ComplexVector V ) {

	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	final double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	final double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;
	newSize( M.numRows );
	MatrixOperations.times(M.re, M.im, VRe, VIm, re, im);
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexVector V, final ComplexMatrix M ) {
	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;

	newSize( M.numCols );
	MatrixOperations.times( VRe, VIm, M.re, M.im, re, im );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexVector V, final RealMatrix M ) {
	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;

	newSize( M.numCols );
	MatrixOperations.times( VRe, M.re, re );
	MatrixOperations.times( VIm, M.re, im );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final ComplexVector V, final IntegerMatrix M ) {
	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;

	newSize( M.numCols );
	MatrixOperations.times( VRe, M.re, re );
	MatrixOperations.times( VIm, M.re, im );
    }

     /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerMatrix M, final ComplexVector V ) {
	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;

	newSize( M.numRows );
	MatrixOperations.times( M.re,VRe, re );
	MatrixOperations.times(  M.re,VIm, im );
    }

     /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealMatrix M, final ComplexVector V ) {
	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	double[] VRe = ( this != V ) ? V.re : (double[])V.re.clone() ;
  	double[] VIm = ( this != V ) ? V.im : (double[])V.im.clone() ;

	newSize( M.numRows );
	MatrixOperations.times( M.re,VRe, re );
	MatrixOperations.times(  M.re,VIm, im );
    }

   /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealVector V, final ComplexMatrix M ) {
	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	MatrixOperations.times( V.re, M.im, im );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final ComplexMatrix M ) {

	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");
	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	MatrixOperations.times( V.re, M.im, im );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealVector V, final RealMatrix M ) {

	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealVector V, final IntegerMatrix M ) {

	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes(  final IntegerMatrix M,final RealVector V ) {

	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final RealMatrix M ) {

	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final IntegerVector V, final IntegerMatrix M ) {

	if( V.length != M.numRows )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numCols );
	MatrixOperations.times( V.re, M.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes(final IntegerMatrix M,final IntegerVector V ) {

	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealMatrix M,  final RealVector V ) {

	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** multiply  <i>V</i> and <i>M</i> and store the result in <i>this</i> vector */
    public void assignTimes( final RealMatrix M,  final IntegerVector V ) {

	if( V.length != M.numCols )
	    throw new IllegalArgumentException("sizes do not Match");

	newSize( M.numRows );
	MatrixOperations.times( M.re, V.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** conjugate <i>V</i> and create and return a vector containing the result */
    public static ComplexVector conjugate( final ComplexVector V ) {
	ComplexVector result = new ComplexVector( V.length );
	result.assignConjugate( V );
	return result;
    }

    /** conjugate <i>this</i> vector and create and return a vector containing the result */
    public ComplexVector conjugate() {
	return conjugate( this );
    }

    /** conjugate <i>V</i> and store the result in <i>this</i> vector */
    public void assignConjugate( final ComplexVector V ) {
	newSize( V.length );

	VectorOperations.assign( V.re, re );
	VectorOperations.neg(    V.im, im );
    }

    /** conjugate <i>this</i> vector */
    public void assignConjugate() {

	VectorOperations.neg( im, im );
    }

    /* return norm of <i>V</i> */
    public static double norm( final ComplexVector V )    {
	return Math.sqrt(V.normSqr());
    }

    /* return norm of <i>this</i> vector */
    public double norm(){
	return Math.sqrt(normSqr());
    }

    /* return square of norm of <i>V</i> */
    public static double normSqr( final ComplexVector V ){
	return VectorOperations.normSqr( V.im ) + VectorOperations.normSqr( V.re );
    }

    /* return square of norm of <i>this</i> vector */
    public double normSqr(){
	return normSqr(this);
    }

    /** return <tt>true</tt> if all entries are zero */
    public static boolean isZero(ComplexVector V) {
      return V.isZero();
    }

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>v</i> */
    public boolean equals ( final ComplexVector v ){
	return equals( v, EPSILON );
    }

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>v</i> */
    public boolean equals ( final ComplexVector v , final double eps ){
      if( v == null || v.length!=length) return false;
      for(int i=0; i<length; i++) {
        if( Math.abs( re[i]-v.re[i] ) > eps ||
            Math.abs( im[i]-v.im[i] ) > eps   )
          return false;
      }
      return true;
    }

    /** return <tt>true</tt> if <i>this</i> vector is equal to <i>o</i> */
    public boolean equals( final Object o){
      if( o==null )
        return false;
	try
	    {
		return equals((ComplexVector)o);
	    }
	catch(ClassCastException ex)
	    {
		return false;
	    }
    }

     /** creates a <b>String</b> containing this' matrix' entries */
    public String toString(){
	StringBuffer sb=new StringBuffer(300);
	sb.append("(");
	if(length>0)
	    {
		sb.append(re[0]);
		if(im[0]>=0.) sb.append('+');
		sb.append(im[0]);
		sb.append('i');
		for(int j=1; j<length; j++)
		    {
			sb.append(", ");
			sb.append(re[j]);
			if(im[j]>=0.) sb.append('+');
			sb.append(im[j]);
			sb.append('i');
		    }
	    }
	sb.append(')');
	return sb.toString();
    }

    /** normalize <i>this</i> vector */
    public void assignNormalize(){
	double n = norm();
	if (n<Double.MIN_VALUE) throw new IllegalArgumentException("Vector has length zero");
	assignDivide(n);

    }

    /** normalize <i>V</i> and store the result in <i>this</i> vector */
    public void assignNormalize( ComplexVector V ){
       	newSize( V );
	double n = V.norm();
	VectorOperations.divide( V.re , n, re );
	VectorOperations.divide( V.im , n, im );
    }

    /** normalize <i>this</i> vector and create and return a vector containing the result */
    public ComplexVector normalize(){
	ComplexVector res = new ComplexVector(this);
	res.assignNormalize();
	return res;
    }

    /** rounds <i>this</i> */
    public void assignRound(){
	assignRound( this );
    }

    /** round <i>V</i> and store the result in <i>this</i> vector */
    public void assignRound( ComplexVector V ){
	newSize( V.size() );
	VectorOperations.round( V.re, re );
	VectorOperations.round( V.im, im );
    }

    /** round <i>V</i> and store the result in <i>this</i> vector */
    public void assignRound( RealVector V ){
	newSize( V.size() );
	VectorOperations.round( V.re, re );
	VectorOperations.assign( im, 0 );
    }

      /** round <i>V</i> and store the result in <i>this</i> vector */
    public void assignRound( IntegerVector V ){
	newSize( V.size() );
	VectorOperations.assign( V.re, re );
	VectorOperations.assign( im, 0 );
    }

    /** round <i>this</i> vector and create and return a vector containing the result */
    public ComplexVector round(){
	ComplexVector res = new ComplexVector(this);
	res.assignRound();
	return res;
    }

    /** round <i>V</i> and create and return a vector containing the result */
    public static ComplexVector round( ComplexVector V ){
	return V.round();
    }

    /** floors <i>this</i> */
    public void assignFloor(){
	assignFloor( this );
    }

    /** floor <i>V</i> and store the result in <i>this</i> vector */
    public void assignFloor( ComplexVector V ){
	newSize( V.size() );
	VectorOperations.floor( V.re, re );
	VectorOperations.floor( V.im, im );
    }

    /** floor <i>V</i> and store the result in <i>this</i> vector */
    public void assignFloor( RealVector V ){
	newSize( V.size() );
	VectorOperations.floor( V.re, re );
	VectorOperations.assign( im, 0 );
    }

      /** floor <i>V</i> and store the result in <i>this</i> vector */
    public void assignFloor( IntegerVector V ){
	newSize( V.size() );
	VectorOperations.assign( V.re, re );
	VectorOperations.assign( im, 0 );
    }

  /** floor <i>this</i> vector and create and return a vector containing the result */
    public ComplexVector floor(){
	ComplexVector res = new ComplexVector(this);
	res.assignFloor();
	return res;
    }

    /** floor <i>V</i> and create and return a vector containing the result */
    public static ComplexVector floor( ComplexVector V ){
	return V.floor();
    }

    /** assign all entries to zero */
    public void assignZero(){
	VectorOperations.assign(re,0.0);
	VectorOperations.assign(im,0.0);
    }

    /** assign all entries randomly */
    public void assignRandom(){
       for( int i = 0; i<length; i++ ){
         re[i] = Math.random();
         im[i] = Math.random();
       }
    }


    /** assign all entries randomly */
    public void assignRandom( final double range){
       for( int i = 0; i<length; i++ ){
         re[i] = Math.random()*range;
         im[i] = Math.random()*range;
       }
    }

   /** creates int vector and assigns it with neg of <code>v</code> */
    public final static ComplexVector neg( final ComplexVector v ) {
	ComplexVector result = new ComplexVector ( v.length );
	VectorOperations.neg( v.re, result.re );
	VectorOperations.neg( v.im, result.im );
	return result;
    }

    /** assign neg of <code>v</code> */
    public final void assignNeg( final ComplexVector v ) {
	newSize( v.length );
	VectorOperations.neg( v.re, re );
	VectorOperations.neg( v.im, im );
    }

     /** assign neg of <code>v</code> */
    public final void assignNeg( final RealVector v ) {
	newSize( v.length );
	VectorOperations.neg( v.re, re );
	VectorOperations.assign( im, 0.0 );
    }

     /** assign neg of <code>v</code> */
    public final void assignNeg( final IntegerVector v ) {
	newSize( v.length );
	VectorOperations.neg( v.re, re );
	VectorOperations.assign( im , 0.0 );
    }

    public final void assignNeg() {
	VectorOperations.neg( re, re );
	VectorOperations.neg( im, im );
    }

}
