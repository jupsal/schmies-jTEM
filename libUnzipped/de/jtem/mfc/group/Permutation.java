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

package de.jtem.mfc.group;




/**
 * Permutation.java 
 *
 * <p>This class deals with the group of permutations of finite sets as
 * objects. It uses the static methods on <code>int[]</code> contained in
 * {@link de.jtem.numericalMethods.algebra.group.Permutation}. The permutation that
 * sends <code>{A,B,C,D}</code> to <code>{D,B,A,C}</code> is encoded as
 * <code>{3,1,0,2}</code> and its <b>inverse</b> <code>{2,1,3,0}</code> is
 * the function from index to index</p> <p>Functions which are sensible to
 * that notation have the <code>Fun</code> form which adopt the functional
 * notation, for example {@link #times} and {@link #timesFun}.</p> 
 * 
 * Created: Thu Nov 21 14:42:03 2002
 *
 * @author <a href="mailto:mercat@sfb288.math.tu-berlin.de">Christian Mercat</a>
 * @version 1.0
 */

public class Permutation {


    private static final long serialVersionUID = 1L;

    /**
     * The data, a permutation of the indices <code>0..n-1</code> where
     * <code>n=s.length</code>. We don't use the functional notation,
     * <code>s[0]</code> is <b>not</b> the index of the image of the first
     * element, <code>{@link #inverse}(s)[0]</code> is.
     *
     */
    private int[] s;

    /** For intermediate computations. Has to be the same size as s.
     */
    private int[] t;

    private void exchangeST(){
	final int[] u = s;
	s=t;
	t=u;
    }
    /**
     * Constructs a new empty <code>Permutation</code> instance.
     *
     */
    public Permutation() {
	assign(new int[0]);	
    }
    /**
     * Constructs a copy of the <code>Permutation</code> instance.
     *
     * @param p a <code>Permutation</code> value
     */
    public Permutation(Permutation p) {
	assign(p);
    }
    /**
     * Constructs a <code>Permutation</code> instance from an
     * <code>int[]</code>.
     *
     * @param s an <code>int[]</code> value
     */
    public Permutation(int[] s) {
	assign(s);
    }
    /**
     * Constructs a <code>Permutation</code> instance from an array of cycles.
     *
     * @param cycles an <code>int[][]</code> value
     */
    public Permutation(int[][] cycles) {
	assign(cycles);
    }


    /**
     * Constructs a new identity <code>Permutation</code> instance of
     * length <code>len</code>..
     *
     * @param len The length of the identity permutation.
     */
    public Permutation(final int len) {
	assignIdentity(len);
    }

    /**
     * Performs {@link #getPermutation permutation}[i]=i on all the
     * entries. Doesn't change its length.
     *
     */
    public void assignIdentity(){
	de.jtem.numericalMethods.algebra.group.Permutation.identity(s);
    }

    /**
     * <code>{@link #assign assign}({@link #identity identity(len)})</code>
     *
     * @param len an <code>int</code> value
     */
    public void assignIdentity(final int len){
	assign(de.jtem.numericalMethods.algebra.group.Permutation.identity(len));
    }

     /**
     * Replaces {@link #getPermutation permutation} by a random permutation
     * of the same length.
     *
     */
    public void assignRandom(){
	de.jtem.numericalMethods.algebra.group.Permutation.random(s);
    }

    /**
     * Replaces {@link #getPermutation permutation} by a random permutation
     * of the given length.  <code>{@link #assign}({@link
     * de.jtem.numericalMethods.algebra.group.Permutation#random
     * random(len)})</code>.
     *
     * @param len an <code>int</code> value
     */
    public void assignRandom(final int len){
	assign(de.jtem.numericalMethods.algebra.group.Permutation.random(len));
    }
  
 
    /**
     * Returns a copy of {@link #getPermutation permutation}.
     * @return the Permutation value.
     */
    public int[] getPermutation() {
	int[] result = new int[s.length];
	System.arraycopy(s,0,result,0,s.length);
	return result;
    }

    /**
     * Sets the permutation.
     * @param newPermutation The new Permutation value.
     */
    public void setPermutation(int[] newPermutation) {
	s=newPermutation;
 	if ((t==null)||(t.length!=s.length)) {
	    t = new int[s.length];
	}
   }

    /**
     * Hands out a copy.
     *
     * @return a <code>Permutation</code> value
     */
    public Permutation copy(){
	return new Permutation(this);
    }
    
    /**
     * Sets the permutation.
     *
     * @param s an <code>int[]</code> value
     */
    public void assign(int[] s){
	setPermutation(s);
	
    }
    /**
     * Sets the permutation.
     *
     * @param s a <code>Permutation</code> value
     */
    public void assign(Permutation s){
	setPermutation(s.getPermutation());
    }
    /**
     * Assigns this to the permutation <b>of the same length</b>, having
     * these cycles. The fixed point don't have to be included. If you don't
     * know the length of the permutation, the fixed point have to be
     * included and you should call <code>{@link #assign}({@link
     * de.jtem.numericalMethods.algebra.group.Permutation#fromCycles}(cycles))</code>.
     *
     * @param cycles an <code>int[][]</code> value
     */
    public void assign(final int[][] cycles){
	de.jtem.numericalMethods.algebra.group.Permutation.fromCycles(cycles,s);
    }

    
    /**
     * Returns the <code>cycles</code> of this permutation. Not optimised.
     *
     * @return The <code>int[][]</code> giving the array of cycles.
     */
    public  int[][] getCycles() {
	return de.jtem.numericalMethods.algebra.group.Permutation.cycles(s);
    }

    /**
     * Assign this to the permutation <b>of the same length</b> having these
     * cycles. Cycles of length <code>1</code> don't have to be included.
     *
     * @param cycles an <code>int[][]</code> value
     */
    public void setCycles(int[][] cycles){
	assign(cycles);
    }

    public void previous() {
	t=s;
	de.jtem.numericalMethods.algebra.group.Permutation.previous(t);
	assign(t);
    }
    public void next() {
	t=s;
	de.jtem.numericalMethods.algebra.group.Permutation.next(t);
	assign(t);
    }
    /**
     * Returns the inverse permutation of this. 
     *
     * @return a <code>Permutation</code> value
     */
    public Permutation inverse(){
	de.jtem.numericalMethods.algebra.group.Permutation.inverse(s,t);
	return new Permutation(t);
    }
  
 
    /**
     * Each permutation can be written as a product of transpositions. Each
     * cycle in the {@link #cycle} decomposition gives rise to an independant
     * product of transposition. The cycle {2,3,...,n,1} can be written as
     * the composition of n-1 transpositions: (1,n)(1,n-1)...(1,2). This
     * method returns the minimal total number of required transpositions.
     *
     * @return an <code>int</code> value
     */
    public int numTranspos(){
	return de.jtem.numericalMethods.algebra.group.Permutation.numTranspos(s,t); }
	
   /**
     * Each permutation can be written as a product of transpositions of
     * <b>consecutive</b> elements. This method returns the minimal total
     * number of required transpositions, sum of the elements in the {@link
     * #inversions} vector.
     *
     * @return an <code>int</code> value
     */
    public int numInversions(){	
	return de.jtem.numericalMethods.algebra.group.Permutation.numInversions(s); }
	
     /**
     * Describe <code>parity</code> method here.
     *
     * @return an <code>int</code> value
     */
    public int parity(){ 
	return de.jtem.numericalMethods.algebra.group.Permutation.parity(s);}

    /**
     * Gives the order of this permutation. It is the least number
     * <code>d</code> such that <code>p^d=Id</code>.
     *
     * @return an <code>int</code> value
     */
    public int order(){	
	return de.jtem.numericalMethods.algebra.group.Permutation.order(s,t); }



    /**
     * Assigns this to its inverse.
     *
     */
    public void assignInvert(){
	de.jtem.numericalMethods.algebra.group.Permutation.inverse(s,t);
	exchangeST();
    }


    /**
     * Assigns this to this composed with the inverse of p.
     *
     * @param p an <code>int[]</code> value
     */
    public void assignDivide(int[] p){
	de.jtem.numericalMethods.algebra.group.Permutation.divide(s,p,t);
	exchangeST();
    }
    /**
     *  Assigns this to this composed with the inverse of p.
     *
     * @param p a <code>Permutation</code> value
     */
    public void assignDivide(Permutation p){
	assignDivide(p.getPermutation());
    }

    /**
     * Returns <code>j</code> such that {@link #getPermutation
     * <code>permutation[j]==i</code>}. Same as <code>{@link
     * #inverse}({@link #getPermutation permutation})[i]</code>. Beware that
     * permutations are not in the functional format but its inverse.
     *
     * @param i an <code>int</code> value
     * @return an <code>int</code> value
     */
    public int applyTo(int i) {
	return de.jtem.numericalMethods.algebra.group.Permutation.applyTo(s,i);
    }
    /**
     * Returns <code>permutation[i]</code>. Beware that
     * permutations are not in the functional format but its inverse.
     *
     * @param i an <code>int</code> value
     * @return an <code>int</code> value
     */
    public int applyToFun(int i) {
	return s[i];
    }

    /**
     * Returns this composed with p.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    public int[] times(int[] p){

	de.jtem.numericalMethods.algebra.group.Permutation.times(s,p,t);
	return (int[]) t.clone();
    }
    /**
     * Returns this composed with p.
     *
     * @param p A Permutation value
     * @return The composition of this with p.
     */
    public Permutation times(Permutation p){
	de.jtem.numericalMethods.algebra.group.Permutation.times
	    (s,p.getPermutation(),t);
	return new Permutation((int[]) t.clone());
    }

    /**
     * Assigns this to this.times(p).
     *
     * @param p an <code>int[]</code> value
     */
    public void assignTimes(int[] p){
	de.jtem.numericalMethods.algebra.group.Permutation.times
	    (s,p,t);
	exchangeST();
    }

    /**
     * Assigns this to this.times(p).
     *
     * @param p a <code>Permutation</code> value
     */
    public void assignTimes(Permutation p){
	de.jtem.numericalMethods.algebra.group.Permutation.times
	    (s,p.getPermutation(),t);
	exchangeST();
    }

    /**
     * Returns this composed with p as functions.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int[]</code> value
     * @param t an <code>int[]</code> of sufficient size.
     */
    public int[] timesFun(int[] p){

	 de.jtem.numericalMethods.algebra.group.Permutation.timesFun(s,p,t);
	 return (int[]) t.clone();
    }
    /**
     * Returns this composed with p as functions.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    public Permutation timesFun(Permutation p){
	de.jtem.numericalMethods.algebra.group.Permutation.timesFun
	    (s,p.getPermutation(),t);
	return new Permutation((int[]) t.clone());
    }

    /**
     * Assigns this to this.times(p) as functions.
     *
     * @param p an <code>int[]</code> value
     */
    public void assignTimesFun(int[] p){
	de.jtem.numericalMethods.algebra.group.Permutation.timesFun
	    (s,p,t);
	exchangeST();
    }

    /**
     * Assigns this to this.times(p).
     *
     * @param p a <code>Permutation</code> value
     */
    public void assignTimesFun(Permutation p){
	de.jtem.numericalMethods.algebra.group.Permutation.timesFun
	    (s,p.getPermutation(),t);
	exchangeST();
    }

   /**
     * Gives the inversion vector, that is, for each index, the number of
     * greater indices to its left in the permutation.
     *
     * @return an <code>int[]</code> value
     */
    public  int[] inversions(){
	 de.jtem.numericalMethods.algebra.group.Permutation.inversions(s,t);
	 return (int[]) t.clone();
    }
    /**
     * Assign this to the Permutation having inv as its inversions vector.
     *
     * @param inv an <code>int[]</code> value
     */
    public  void fromInversions(final int[] inv){
	de.jtem.numericalMethods.algebra.group.Permutation.fromInversions
	    (inv,s);
    }
    public  int[][][] youngTableaux(){
	return de.jtem.numericalMethods.algebra.group.Permutation.youngTableaux(s);
    }
     public void assignFromYoungTableaux(final int[][][] y){
	assign
	    (de.jtem.numericalMethods.algebra.group.Permutation.fromYoungTableaux
	     (y));
    }
    public String cyclesToString() {
	return de.jtem.numericalMethods.algebra.group.Permutation.cyclesToString(s);
    }

    /**
     * Describe <code>toString</code> method here.
     *
     * @return a <code>String</code> value
     */
    public String toString() {
	return de.jtem.numericalMethods.algebra.group.Permutation.toString(s);
    } 
}// Permutation
