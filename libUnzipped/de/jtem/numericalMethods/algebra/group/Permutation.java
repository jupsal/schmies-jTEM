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

package de.jtem.numericalMethods.algebra.group;


import java.io.StreamTokenizer;
import java.io.StringReader;


/**
 * Permutation.java 
 *
 * <p>This class implements some static algorithms that deal with the group
 * of permutations of finite sets. The permutation that sends
 * <code>{A,B,C,D}</code> to <code>{D,B,A,C}</code> is encoded as
 * <code>{3,1,0,2}</code> and its <b>inverse</b> <code>{2,1,3,0}</code> is
 * the function from index to index.</p> <p>Functions which are sensible to
 * that notation have the <code>Fun</code> form which adopt the functional
 * notation, for example {@link #times} and {@link #timesFun}.</p> <p>Methods
 * which require an intermediate array (whether <code>int[]</code> or
 * <code>boolean[]</code> of the same (or greater) length as the permutation)
 * are provided with a signature that allows an instance to be passed on so
 * that no intermediate instances are created. It is not thread safe and the
 * values in these arrays are destroyed! In particular don't use twice the
 * permutation array as the pair of parameters.</p>
 *
 * Created: Thu Nov 21 14:42:03 2002
 *
 * @author <a href="mailto:mercat@sfb288.math.tu-berlin.de">Christian Mercat</a>
 * @version
 */

public class Permutation {


    private static final long serialVersionUID = 1L;

    /**
     * Returns an identity permutation of the desired length.
     *
     * @param len an <code>int</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] identity(final int len){
	final int[] result = new int[len];
	identity(result);
	return result;
    }
    /**
     * Replaces the <code>int[]</code> parameter with the identity of the
     * same length.
     *
     * @param s an <code>int[]</code> value
     */
    static public void identity(final int[] s){
	for (int  i = s.length;  --i>=0 ; )  s[i]=i;
    }

    /**
     * Returns a random permutation on the letters
     * <code>0..len-1</code>.
     *
     * @param len an <code>int</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] random(final int len) {
	int [] s=new int[len];
	random(s);
	return s;    
    }
    /**
     * Returns a random permutation on the letters <code>0..len-1</code>. If
     * <code>derangement</code> is non zero, the result will be without fixed
     * point.
     *
     * @param len an <code>int</code> value
     * @param derangement an <code>int</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] random(final int len, final int derangement) {
	int [] s=new int[len];
	random(s,derangement);
	return s;    
    }
    /**
     * Replaces the parameter array with a random array  on the letters
     * <code>0..len-1</code> where <code>len=s.length</code>.
     *
     * @param s an <code>int[]</code> value
     */
    static public void random(final int[] s) {
	random(s,0);
    }
    /**
     * Replaces the parameter array with a random array on the letters
     * <code>0..len-1</code> where <code>len=s.length</code>. If the second
     * parameter is non zero, the result is a derangement, that is a
     * permutation with no fixed point.
     *
     * @param s an <code>int[]</code> value
     * @param derangement an <code>int</code> value
     */
    static public void random(final int[] s, final  int derangement) {
	
	final boolean d = (derangement!=0);
	identity(s);	
	
	final	java.util.Random rnd = new java.util.Random();
	for (int i = s.length; i > 1;)
	    {
		// if derangement, i can not be chosen and has to change.
		final int r = rnd.nextInt(d?--i:i--); 

		final int copy = s[i];
		s[i] = s[r];
		s[r] = copy;
	    }
	
    }
  
    /**
     * Returns whether the given <code>int[]</code> is a valid permutation of
     * the indices <code>0..length-1</code>.
     *
     * @param p an <code>int[]</code> value
     * @return a <code>boolean</code> value
     * @see #isPermutation(int[] p, boolean[] flags)
     */
    static public boolean isPermutation(final int [] p){
	final boolean[] flags = new boolean[p.length];
	return isPermutation(p,flags);
    }
    /**
 * Returns whether the given <code>int[]</code> is a valid permutation of the
 * indices <code>0..length-1</code>. Uses the <code>garbage int[]</code>
 * instance as a flag. Its values are altered and it has to be different
 * from <code>p</code>.
 *
 * @param p an <code>int[]</code> value
 * @param garbage an <code>int[]</code> of at least the same length.
 * @return a <code>boolean</code> value
 */
    static public boolean isPermutation(final int [] p, final int [] garbage){
	if (p.length>garbage.length) {
	    throw new IllegalArgumentException
		("mfc.de.jtem.numericalMethods.group.Permutation.isPermutation "
		 +"called with int[] of insufficient size.");
	}
	

	java.util.Arrays.fill(garbage,1);
	for (int i=p.length;--i>=0;) {
	    try {
		garbage[p[i]]-=1;
	    } catch (ArrayIndexOutOfBoundsException e) {
		return false;
	    }
	} 
	for (int i=p.length;--i>=0;) {
	    if (garbage[i]!=0) return false;
	} 
	return true;
    }
/**
 * Returns whether the given <code>int[]</code> is a valid permutation of
 * the indices <code>0..length-1</code>. Uses the <code>flags
 * boolean[]</code> instance. Its value is altered.
 *
 * @param p an <code>int[]</code> value
 * @param flags a <code>boolean[]</code> of at least the same length.
 * @return a <code>boolean</code> value
 */
static public boolean isPermutation(final int [] p, 
				    final boolean [] flags){
    if (p.length>flags.length) {
	throw new IllegalArgumentException
	    ("mfc.de.jtem.numericalMethods.group.Permutation.isPermutation "
	     +"called with int[] of insufficient size.");
    }
	

    java.util.Arrays.fill(flags,true);
    for (int i=p.length;--i>=0;) {
	try {
	    if (flags[p[i]]){
		flags[p[i]]=false;} else {
		    return false;
		}
	} catch (ArrayIndexOutOfBoundsException e) {
	    return false;
	}
	    
    } 
    for (int i=p.length;--i>=0;) {
	if (flags[i]) return false;
    } 
    return true;
}
    /**
     * Returns whether the given <code>int[]</code> is a valid derangement
     * (permutation without fixed point) of the indices
     * <code>0..length-1</code>.
     *
     * @param p an <code>int[]</code> value
     * @return a <code>boolean</code> value
     * @see #isDerangement(int[] p, boolean[] flags)
     */
    static public boolean isDerangement(final int [] p){
	final boolean[] flags = new boolean[p.length];
	return isDerangement(p,flags);
    }
    /**
 * Returns whether the given <code>int[]</code> is a valid derangement
 * (permutation without fixed point) of the indices
 * <code>0..length-1</code>. Uses the <code>garbage int[]</code> instance as
 * a flag. Its values are altered and it has to be different from
 * <code>p</code>.
 *
 * @param p an <code>int[]</code> value
 * @param garbage an <code>int[]</code> of at least the same length.
 * @return a <code>boolean</code> value
 */
    static public boolean isDerangement(final int [] p, final int [] garbage){
    for (int  i = p.length; --i>=0;) {
	if (p[i]==i) return false;
    }

    return isPermutation(p,garbage);    

    }
/**
 * Returns whether the given <code>int[]</code> is a valid derangement
 * (permutation without fixed point) of the indices
 * <code>0..length-1</code>. Uses the <code>flags boolean[]</code>
 * instance. Its value is altered.
 *
 * @param p an <code>int[]</code> value
 * @param flags a <code>boolean[]</code> of at least the same length.
 * @return a <code>boolean</code> value
 */
static public boolean isDerangement(final int [] p, 
				    final boolean [] flags){
    for (int  i = p.length; --i>=0;) {
	if (p[i]==i) return false;
    }

    return isPermutation(p,flags);    
}

/**
 * Returns the <code>cycles</code> of a given permutation, as
 * <code>int</code> arrays.
 *
 * @param s The <code>int[]</code> permutation.
    * @return The <code>int[][]</code> giving the array of cycles.
    */
    static public int[][] cycles(final int[] s) {

	final int n = s.length;
	// Whether it belongs already in a cycle.
	final boolean[] flag = new boolean[n];
	final int[][] c = new int[n][];
	final int[] grow = new int[n]; 

	int number = 0;

	for (int i=0; i<n; i++) {
	    if (!flag[i]){
		int len = 0;
		try {
		    for (int j=s[i]; !flag[j]; j=s[j]) {
			len++;
			grow[n-len]=j;
			flag[j]=true;
			} 
		    } catch (ArrayIndexOutOfBoundsException e) {
		    throw new IllegalArgumentException
		    ("Permutation.cycles called with a non valid "
		     +"permutation of its indices 0..n-1."+e.getMessage());
		    } 
		if (grow[n-len]!=i) {// Not a cycle!
		    throw new IllegalArgumentException
		    ("Permutation.cycles called with a non valid "
		     +"permutation of its indices 0..n-1.");
		    }
		c[number]=new int[len];
		System.arraycopy(grow,n-len,c[number],0,len);
		number++;
		}
	    }
	final int[][] result = new int[number][];
	System.arraycopy(c,0,result,0,number);

	return result;
	}
/**
 * Returns the <code>cycles</code> of a given permutation understood as a
 * function of the indices, as <code>int</code> arrays.
 *
 * @param s The <code>int[]</code> permutation.
    * @return The <code>int[][]</code> giving the array of cycles.
    */
    static public int[][] cyclesFun(final int[] s) {

	final int n = s.length;
	// Whether it belongs already in a cycle.
	final boolean[] flag = new boolean[n];
	final int[][] c = new int[n][];
	final int[] grow = new int[n]; 

	int number = 0;

	for (int i=0; i<n; i++) {
	    if (!flag[i]){
		int len = 0;
		try {
		    for (int j=i; !flag[j]; j=s[j]) {
			grow[len++]=j;
			flag[j]=true;
			} 
		    } catch (ArrayIndexOutOfBoundsException e) {
		    throw new IllegalArgumentException
		    ("Permutation.cycles called with a non valid "
		     +"permutation of its indices 0..n-1."+e.getMessage());
		    } 
		if (s[grow[len-1]]!=i) {// Not a cycle!
		    throw new IllegalArgumentException
		    ("Permutation.cycles called with a non valid "
		     +"permutation of its indices 0..n-1.");
		    }
		c[number]=new int[len];
		System.arraycopy(grow,0,c[number],0,len);
		number++;
		}
	    }
	final int[][] result = new int[number][];
	System.arraycopy(c,0,result,0,number);

	return result;
	}

    /**
     * Changes <code>result</code> to a permutation of the same length, built
     * from the given cycles. Cycles of length <code>1</code> don't have to
     * be included.
     *
     * @param c an <code>int[][]</code> value
     * @param result an <code>int[]</code> value
     * @exception ArrayIndexOutOfBoundsException if an error occurs
     */
    static public void fromCycles(final int[][] c, final int[] result) 
    throws ArrayIndexOutOfBoundsException{

	final int len = result.length;
	identity(result);

	for (int  i = 0;  i<c.length ; i++) {
	    final int last = c[i].length-1;
	    result[c[i][0]]=c[i][last];
	    for (int  k = last ; k>0;) {
		result[c[i][k]]=c[i][--k];
		}
	    } 
	}
    /**
     * Changes <code>result</code> to a permutation of the same length, built
     * from the given cycles, as a function of the indices. Cycles of length
     * <code>1</code> don't have to be included.
     *
     * @param c an <code>int[][]</code> value
     * @param result an <code>int[]</code> value
     * @exception ArrayIndexOutOfBoundsException if an error occurs
     */
    static public void fromCyclesFun(final int[][] c, final int[] result) 
    throws ArrayIndexOutOfBoundsException{

	final int len = result.length;
	identity(result);

	for (int  i = 0;  i<c.length ; i++) {
	    final int last = c[i].length-1;
	    result[c[i][last]]=c[i][0];
	    for (int  k = 0 ; k<last;) {
		result[c[i][k]]=c[i][++k];
		}
	    } 
	}
    /**
     * Returns a permutation built from its cycles. Fixed points have to be
     * included.
     *
     * @param c an <code>int[][]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] fromCycles(final int[][] c){
	int len=0;
	for (int  i = 0;  i<c.length ; i++) {
	    len+=c[i].length;
	    } 
	final int[] result = new int[len];

	try {
	    fromCycles(c,result);
	    } catch (ArrayIndexOutOfBoundsException e) {
	    throw new IllegalArgumentException
	    ("Permutation.fromCycles called with a non valid "
	     +"array of cycles.");
	    }
	
	return result; 
	}

     /**
     * Returns the inverse permutation. It is the
     * <code>int[] p</code> such that the index of <code>i</code> in
     * <code>s</code> is <code>p[i]</code>.
     *
     * @param s an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] inverse(int[] s){
	final int[] result = new int[s.length];
	inverse(s,result);
	return result;   
	} 
 
    /**
     * Changes <code>result</code> to the inverse permutation of
     * <code>s</code>. It is the <code>int[] result</code> such that the
     * index of <code>i</code> in <code>s</code> is
     * <code>result[i]</code>. No checks are performed to test whether the
     * parameter was indeed a permutation, but may throw
     * IllegalArgumentException's if it was not.
     *
     * @param s an <code>int[]</code> value
     * @param result an <code>int[]</code> of the same (or greater)
     * length. Has to be different from <code>s</code>.
     */
    static public void inverse(final int[] s, final int[] result) 
	throws IllegalArgumentException {

// 	java.util.Arrays.fill(result,-1);

	for (int  i = s.length ; --i>=0;) {
	    try {
		    result[s[i]]=i;
	    } catch (ArrayIndexOutOfBoundsException e) {
		throw new  IllegalArgumentException
		("Permutation.inverse called with a non valid "
		 +"permutation of its indices 0..n-1.");
	    }
		}

// 	for (int  i = s.length ; --i>=0;) {
// 	    if (result[i]==-1) {
// 		throw new IllegalArgumentException
// 		("Permutation.inverse called with a non valid "
// 		 +"permutation of its indices 0..n-1.");
// 		}
// 	    } 
	

	} 
 
    /**
     * Each permutation can be written as a product of transpositions. Each
     * cycle in the {@link #cycle} decomposition gives rise to an independant
     * product of transposition. The cycle {2,3,...,n,1} can be written as
     * the composition of n-1 transpositions: (1,n)(1,n-1)...(1,2). This
     * method returns the minimal total number of required transpositions.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int</code> value
     */
    static public int numTranspos(final int[] p){
	return numTranspos(p,new int[p.length]);	
	}
    /**
     * Same as {@link #numTranspos} but uses a preinstanciated
     * <code>int[]</code> of length at least the length of <code>p</code>. It
     * will be filled with <code>1</code>.
     *
     * @param p an <code>int[]</code> value
     * @param garbage an <code>int[]</code> such that <code>garbage.length
     * >= p.length</code>.
     * @return an <code>int</code> value
     */

    static public int numTranspos(final int[] p, final int[] garbage){
	java.util.Arrays.fill(garbage,0);

	int result = 0;
	for (int  i = p.length; --i>=0;){
	int j=i;
	if (garbage[i]==0) {
	do{
	    garbage[j]=1;
	    j=p[j];
	    result++;
	    } while (garbage[j]==0);
	result--;
	}
	}
	
	return result;	
	}
    /**
     * Describe <code>numTranspos</code> method here.
     *
     * @param p an <code>int[]</code> value
     * @param flags a <code>boolean[]</code> value
     * @return an <code>int</code> value
     */
    static public int numTranspos(final int[] p, final boolean[] flags){
	java.util.Arrays.fill(flags,true);

	int result = 0;
	for (int  i = p.length; --i>=0;){
	int j=i;
	if (flags[i]) {
	do{
	    flags[j]=false;
	    j=p[j];
	    result++;
	    } while (flags[j]);
	result--;
	}
	}
	
	return result;	
	}
   
    /**
     * Each permutation can be written as a product of transpositions of
     * <b>consecutive</b> elements. This method returns the minimal total
     * number of required transpositions, sum of the elements in the {@link
     * #inversions} vector.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int</code> value
     */
    static public int numInversions(final int[] p){
	int result=0;
	final int len = p.length;
	for (int  i = len;--i>=0;) {
	    for (int  j = 0;(p[j]!=i)&&(j<len);j++) {
		if (p[j]>i) result++;
		}
	    }
	
	return result;
	}
 
    /**
     * Gives the parity of the number of transpositions.
     *
     * @param p an <code>int[]</code> value
     * @return an <code>int</code> value
     * @see #numInversions
     */
    static public int parity(int[] p){ return numInversions(p)%2;}

    /**
     * Returns the order of the given permutation. It is the least number
     * <code>d</code> such that <code>p^d=Id</code>.
     *
     * @param p an <code>int[]</code> permutation.
     * @return Its order
     */
    static public int order(int[] p){
	return order(p,new boolean[p.length]);
	}

    /**
     * Same as {@link #order} but uses a preinstanciated
     * <code>int[]</code> of length at least the length of <code>p</code>. It
     * will be filled with <code>1</code>.
     *
     * @param p an <code>int[]</code> value
     * @param garbage an <code>int[]</code> such that <code>garbage.length
     * >= p.length</code>.
     * @return an <code>int</code> value
     */

    static public int order(final int[] p, final int[] garbage){
	java.util.Arrays.fill(garbage,0);

	int result = 1;
	for (int  i = p.length; --i>=0;){
	int j=i;

	if (garbage[i]==0) {
	int len = 0;
	do{
	    garbage[j]=1;
	    j=p[j];
	    len++;
	    } while (garbage[j]==0);
	result=lcm(result,len);
	}
	}
	
	return result;	
	}
    /**
     * Describe <code>order</code> method here.
     *
     * @param p an <code>int[]</code> value
     * @param flags a <code>boolean[]</code> value
     * @return an <code>int</code> value
     */
    static public int order(final int[] p, final boolean[] flags){
	java.util.Arrays.fill(flags,true);

	int result = 1;
	for (int  i = p.length; --i>=0;){
	int j=i;

	if (flags[i]) {
	int len = 0;
	do{
	    flags[j]=false;
	    j=p[j];
	    len++;
	    } while (flags[j]);
	result=lcm(result,len);
	}
	}
	
	return result;	
	}

    /**
     * Computes the least common multiple of two numbers. This should
     * definitely belong elsewhere!
     *
     * @param m An <code>int</code> value
     * @param n An <code>int</code> value
     * @return Their least common multiple.
     */
    static public int lcm(int m, int n) {
	return m*n/gcd(m,n);
	}
    /**
     * Computes the greatest common divisor of two numbers. This should
     * definitely belong elsewhere!
     *
     * @param m An <code>int</code> value
     * @param n An <code>int</code> value
     * @return The greatest integer that divides them both.
     */
    static public int gcd(int m, int n) {
	int t;

	if(m<0){m=-m;}
	if(n<0){n=-n;}
	while( m > 0 )
 { /* invariant: gcd(m,n)=g */
		if( n > m )
 {t = m; m = n; n = t; } /* swap */
		/* m >= n > 0 */
		m -= n;
		}
	return n;
	}
    /**
     * Computes the  least common multiple of several numbers. 
     *
     * @param nums an <code>int[]</code> value
     * @return an <code>int</code> value
     */
    static public int lcm(int[] nums) {
	int i = nums.length-1;
	int d = nums[i];
	for (;--i>=0;) {
	    d=lcm(d,nums[i]);
	    }
	return d;
	}
    /**
     * Computes the greatest common divisor of several numbers. 
     *
     * @param nums an <code>int[]</code> value
     * @return an <code>int</code> value
     */
    static public int gcd(int[] nums) {
	int i = nums.length-1;
	int d = nums[i];
	for (;--i>=0;) {
	    d=gcd(d,nums[i]);
	    }
	return d;
	}
    /**
     * <code>factorial(n)=n!=n*(n-1)*...*1</code>. It is the number of
     * permutations on <code>n</code> letters.
     *
     * @param n an <code>int</code> value
     * @return an <code>int</code> value
     */
    static public int factorial(final int n){
	int result = 1;
	for (int  i = n;  i>0;i--) {
	    result *=i;
	}
	return result;
    }	
    /**
     * <code>subFactorial(n+1)=!(n+1)=n*(!n+!(n-1))=[n!/E]</code>. It is the
     * number of deplacements on <code>n</code> letters.
     *
     * @param n an <code>int</code> value
     * @return an <code>int</code> value
     */
    static public int subFactorial(final int n){
// 	return (int) (.5+factorial(n)/java.lang.Math.E);
	int sNm1 = 0;
	int sN = (n>1)?1:0;
	for (int i=2; i<n;i++) {
	    final int sC = sN;
	    sN = i*(sN+sNm1);
	    sNm1 = sC;
	}
	return sN;
	
    }	
//     static public int subFactorial2(final int n){
// 	return (int) (.5+factorial(n)/java.lang.Math.E);
//     }	
    /**
     * Replaces a permutation with the next permutation in the cyclic
     * lexicographic order. If you want to iterate on all permutations of a
     * given size <code>n</code>, begin with a permutation and call
     * <code>{@link #factorial n!}</code> times this method.
     *
     * @param p an <code>int[]</code> value
     * @see #factorial
     */
    static	public	void next(int[] p){
	final int last = p.length-1;
	int current = p[last];
	for (int  i = last; --i>=0;) {
	    // Find the first descent starting from the end.
	    if (p[i]<current) {
		int pi = p[i];
		int j=p.length;
		// Find the smallest entry larger than p[i].
		while (p[--j]<pi) {}
		p[i]=p[j];
		p[j]=pi;
		// Flip the rest of the tail in an increasing order.
		final int crossing = last+i+1;
		for (j = (crossing>>1); j>i ; j--) {
		    final int reflect = crossing -j;
		    pi = p[j];
		    p[j]=p[reflect];
		    p[reflect]=pi;
		    }
		return;	
		}
	    current = p[i];
	    }
	// If the loop completes, then it was the last permutation. Return
	// the first one (identity). For a non cyclic order, we should
	// nullify p.
	identity(p);	
	return;
	}
    /**
     * Replaces a derangement (permutation with no fixed point) with the next
     * derangement in the cyclic lexicographic order. If you want to iterate
     * on all derangements of a given size <code>n</code>, begin with a
     * permutation and call <code>{@link #subFactorial !n}</code> times this
     * method.
     *
     * @param p an <code>int[]</code> derangement (not checked!).
     * @see #subFactorial
     * @see #firstDerangement
     */
    static	public	void nextDerangement(int[] p){
	final int last = p.length-1;
	int current = p[last];
	iLoop:
	for (int  i = last; --i>=0;) {
	    // Find the first descent starting from the end.
	    if (p[i]<current) {
		int pi = p[i];
		int j=p.length;
		// Find the smallest entry larger than p[i].
		while (p[--j]<pi) {}
		if (p[j]!=i) {
		    p[i]=p[j];
		    p[j]=pi;
		} else {
		    if (j>i+1) {
			p[j]=pi;
			p[i]=p[--j];
			p[j]=i;
		    } else {
			p[i]=p[j];
			p[j]=pi; 
			// It is no longer a derangement but in decreasing
			// order.
			current = p[i];
			continue iLoop;	
		    }
		    
		}
		
		// Flip the rest of the tail in an increasing order.
		final int crossing = last+i+1;
		for (j = (crossing>>1); j>i ; j--) {
		    final int reflect = crossing -j;
		    pi = p[j];
		    p[j]=p[reflect];
		    p[reflect]=pi;
		}
		// Enforce that it's a derangement.
		for (j = i; j<last ;j++) {
		    if (p[j]==j) {
			pi = p[j];
			p[j]=p[++j];
			p[j]=pi;
		    }
		}
		if (p[last]==last) {
		    p[last]=p[last-1];
		    p[last-1]=last;
		}

		return;	
	    }
	    current = p[i];
	}
	// If the loop completes, then it was the last permutation. Return
	// the first one. For a non cyclic order, we should
	// nullify p.
	firstDerangement(p);	
	return;
    }
    /**
     * Let the parameter be the <code>firstDerangement</code> (permutation
     * without fixed point) in the lexicographic order.
     *
     * @param p an <code>int[]</code> value
     */
    static public void firstDerangement(final int[] p){
	final int last = p.length-1;
	for (int  i = 0;i<last; i++){
	    p[i]=++i;
	    p[i]=i-1;
	}
	if (last%2==0) {
	    p[last]=last-2;
	    p[last-1]=last;
	    p[last-2]=last-1;
	}
    }
    /**
     * Flips the permutation left-right.
     *
     * @param s an <code>int[]</code> value
     */
    static public void flipLR(int[] s){
	final int last = s.length-1;
	for (int  i =(s.length>>1); i-->=0 ;) {
	    final int bucket = s[i];
	    final int flippedI = last-i;
	    s[i]=s[flippedI];
	    s[flippedI]=bucket;
	}
	
	} 
    /**
     * Flips the permutation up-down. <code>s[i]=s.length-1-s[i]</code>.
     *
     * @param s an <code>int[]</code> value
     */
    static public void flipUD(int[] s){
	final int last = s.length-1;
	for (int  i =s.length;  --i>=0 ;) {
	    s[i]=last-s[i];
	}
	
	} 
    /**
     * The previous one in the lexicographic order. Conjugates {@link #next}
     * by {@link #flipUD}.
     *
     * @param p an <code>int[]</code> value
     */
    static	public	void previous(int[] p){
	flipUD(p);
	next(p);
	flipUD(p);
    }
    /**
     * Returns <code>j</code> such that <code>s[j]==i</code>. Same
     * as <code>{@link #inverse}(s)[i]</code>. Beware that
     * permutations are not in the functional format but its inverse.
     *
     * @param s an <code>int[]</code> value
     * @param i an <code>int</code> value
     * @return an <code>int</code> value
     * @see #applyToFun
     */
    static public int applyTo(int[] s,int i) {
	for (int j=s.length;--j>=0;) {
	    if (s[j]==i) return j;
	    }
	throw new IllegalArgumentException
	("Permutation.applyTo called with a non valid "
	 +"index.");
	}
    /**
     * Returns <code>s[i]</code>. Beware that
     * permutations are not in the functional format but its inverse.
     *
     * @param s an <code>int[]</code> value
     * @param i an <code>int</code> value
     * @return an <code>int</code> value
     * @see #applyTo
     */
    static public int applyToFun(int[] s,int i){return s[i];} 

    static public int[] times(int[] p1, int[] p2){
	final  int[] result=new int[p1.length];
	times(p1,p2,result);
	return result;
    }
    /**
     * Composes two permutations. Gives back <code>p2[p1[i]]</code> and
     * <b>not</b> <code>p1[p2[i]]</code> since the permutation notation is
     * contra-functional. Consider {@link #timesFun} for the functional
     * notation composition.
     *
     * @param p1 an <code>int[]</code> permutation.
     * @param p2 an <code>int[]</code> permutation of the same length.
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
     * @see #timesInverse
     * @see #timesFun
     */
    static public void times(int[] p1, int[] p2, int[] result){
	if (p1.length!=p2.length) {
	    throw new IllegalArgumentException
	    ("Permutation.times called with permutations "
	     +"of different sizes.");
	    }
	
	for (int i = p1.length; --i>=0;) {
	    try {
	    result[i]=p2[p1[i]];
	    } catch (ArrayIndexOutOfBoundsException e) {
	    throw new IllegalArgumentException
	    ("Permutation.times called with invalid permutations.");
	    }
	    }
	return;
	}
     static public int[] timesInverse(int[] p1, int[] p2){
	final  int[] result=new int[p1.length];
	timesInverse(p1,p2,result);
	return result;
    }
   /**
     * Composes two permutations as functions and returns the permutation
     * (contra-functional) notation in <code>result</code>. It is the same as
     * <code>{@link #inverse}(result)</code> after <code>{@link
     * #times}({@link #inverse}(p1),{@link #inverse}(p2),result)</code>.
     *
     * @param p1 an <code>int[]</code> permutation.
     * @param p2 an <code>int[]</code> permutation of the same length.
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
     * @see #timesFun
     */
    static public void timesInverse(int[] p1, int[] p2, int[] result){
	if (p1.length!=p2.length) {
	    throw new IllegalArgumentException
	    ("Permutation.times called with permutations "
	     +"of different sizes.");
	    }
	
	for (int i = p1.length; --i>=0;) {
	    try {
	    result[p1[p2[i]]]=i;
	    } catch (ArrayIndexOutOfBoundsException e) {
	    throw new IllegalArgumentException
	    ("Permutation.timesInverse called with invalid permutations.");
	    }
	}
	return;
	}
    static public int[] timesFun(int[] p1, int[] p2){
	final  int[] result=new int[p1.length];
	timesFun(p1,p2,result);
	return result;
    }
    /**
     * Returns <code>p1[p2[i]]</code> in <code>result</code>.
     *
     * @param p1 an <code>int[]</code> permutation.
     * @param p2 an <code>int[]</code> permutation of the same length.
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
     */
    static public void timesFun(int[] p1, int[] p2, int[] result){
	if (p1.length!=p2.length) {
	    throw new IllegalArgumentException
	    ("Permutation.times called with permutations "
	     +"of different sizes.");
	    }
	
	for (int i = p1.length; --i>=0;) {
	    try {
	    result[i]=p1[p2[i]];
	    } catch (ArrayIndexOutOfBoundsException e) {
	    throw new IllegalArgumentException
	    ("Permutation.timesInverse called with invalid permutations.");
	    }
	}
	return;
	}

    static public int[] divide(int[] p1, int[] p2){
	final  int[] result=new int[p1.length];
	divide(p1,p2,result);
	return result;
    }
    /**
     * Puts <code>p1</code> composed with <code>{@link #inverse}(p2)</code>
     * into the third argument. Beware that the permutation notation is
     * contra-functional. Not well optimised.
     *
     * @param p1 an <code>int[]</code> permutation.
     * @param p2 an <code>int[]</code> permutation of the same length.
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
      * @see #divideFun
     */
    static public void divide(int[] p1, int[] p2, int[] result){
	if (p1.length!=p2.length) {
	    throw new IllegalArgumentException
	    ("Permutation.divide called with permutations "
	     +"of different sizes.");
	    }
	
	for (int i = p1.length; --i>=0;) {

		for (int  j = p1.length; --j>=0;){
		    if (p2[j]==p1[i]) {
			result[i]=j;
			break;
		    }
		    
		}
	    }
	return;
	}
   static public int[] divideFun(int[] p1, int[] p2){
	final  int[] result=new int[p1.length];
	divideFun(p1,p2,result);
	return result;
    }
    /**
     * Puts <code>p1[p2^-1[i]]</code> into the third argument.
     *
     * @param p1 an <code>int[]</code> permutation.
     * @param p2 an <code>int[]</code> permutation of the same length.
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
      * @see #divideFun
     */
    static public void divideFun(int[] p1, int[] p2, int[] result){
	if (p1.length!=p2.length) {
	    throw new IllegalArgumentException
	    ("Permutation.divideFun called with permutations "
	     +"of different sizes.");
	    }
	
	for (int i = p1.length; --i>=0;) {
	    try {
	    result[p2[i]]=p1[i];
	    } catch (ArrayIndexOutOfBoundsException e) {
	    throw new IllegalArgumentException
	    ("Permutation.divideFun called with invalid permutations.");
	    }
	    }
	return;
	}
    
    /**
     * Gives the inversion vector, that is, for each index, the number of
     * greater indices to its left in the permutation.
     *
     * @param s an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] inversions(final int[] s){
	final int len = s.length;

	final int[] result = new int[len];
	inversions(s,result);

	return result;
    }	
    /**
     * Gives the inversion vector, that is, for each index, the number of
     * greater indices to its left in the permutation.
     *
     * @param s an <code>int[]</code> value
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
     */
    static public void inversions(final int[] s,final int[] result){
	final int len = s.length;


	for (int  i = len; --i>=0;) {
	    for (int  j = 0;(s[j]!=i)&&(j<len);j++) {
		if (s[j]>i) result[i]++;
		}
	    } 
	return ;
    }	
    
    /**
     * Gives the inversion vector of a permutation in functional notation,
     * that is, for each index, the number of greater indices to its left in
     * the permutation.
     *
     * @param s an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] inversionsFun(final int[] s){
	final int len = s.length;

	final int[] result = new int[len];
	inversionsFun(s,result);

	return result;
    }
    /**
     * Gives the inversion vector of a permutation in functional notation,
     * that is, for each index, the number of greater indices to its left in
     * the permutation.
     *
     * @param s an <code>int[]</code> value
     * @param result an <code>int[]</code> of the same (or greater) length
     * that will hold the result.
     */
    static public void inversionsFun(final int[] s, final int[] result){
	final int len = s.length;


	for (int  i = len; --i>=0;) {
	    for (int  j = len;--j>i;) {
		if (s[j]<i) result[i]++;
		}
	    } 
	return;
    }	
	

    /**
     * Returns the permutation having inv as its inversions vector.
     *
     * @param inv an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] fromInversions(final int[] inv){
	final int[] s = new int[inv.length];
	fromInversions(inv,s);
	return s;
    }
   /**
     * Returns the permutation having inv as its inversions vector.
     *
     * @param inv an <code>int[]</code> value
     * @param s the <code>int[]</code> result, has to be of the same (or
     * greater) length.
     */
    static public void fromInversions(final int[] inv, final int[] s){
	final int len = inv.length;
	final java.util.ArrayList sL = new java.util.ArrayList(len);

	for (int  i = len; --i>=0;) {
	    sL.add(inv[i],new Integer(i));
	    } 

	for (int  i = len; --i>=0;) s[i]=((Integer) sL.get(i)).intValue();
	return;	
	}
    /**
     * Returns the permutation, in functional notation, having inv as its
     * inversions vector.
     *
     * @param inv an <code>int[]</code> value
     * @return an <code>int[]</code> value
     */
    static public int[] fromInversionsFun(final int[] inv){
	final int[] s = new int[inv.length];
	fromInversionsFun(inv,s);
	return s;
    }
    /**
     * Returns the permutation, in functional notation, having inv as its
     * inversions vector. It is the inverse permutation of {@link
     * fromInversions} with the same inversion vector.
     *
     * @param inv an <code>int[]</code> value
     * @param s the <code>int[]</code> result, has to be of the same (or
     * greater) length.
     */
    static public void fromInversionsFun(final int[] inv, final int[] s){
	final int len = inv.length;
	final java.util.ArrayList sL = new java.util.ArrayList(len);

	for (int  i = len; --i>=0;) {
	    sL.add(inv[i],new Integer(i));
	    } 

	for (int  i = len; --i>=0;) s[((Integer) sL.get(i)).intValue()]=i;
	return;	
	}
    /**
     * Gives the pair of <code>int[][]</code> Young tableaux associated with
     * the <code>int[]</code> permutation. Schensted's correspondance.
     *
     * @param s an <code>int[]</code> value
     * @return an <code>int[][][]</code> value
     */
    static public int[][][] youngTableaux(final int[] s){
	final int len = s.length;
	final java.util.ArrayList yP = new java.util.ArrayList(len);
	final java.util.ArrayList yQ = new java.util.ArrayList(len);


	for (int  i = len; --i>=0;) {
	    yP.add(0, new java.util.ArrayList(i));
	    yQ.add(0, new java.util.ArrayList(i));
	    } 
	// 	((java.util.ArrayList) yL.get(0)).add(new Integer(s[0]));


	for (int  i = 0; i<s.length ; i++) {
	    int si = s[i];
	    // 	    System.out.println("si="+si);
	    
	    jkloops:
	    for (int j = 0, limj = yP.size();  j<limj ; j++) {


		final java.util.ArrayList yPj = 
		(java.util.ArrayList) yP.get(j);
		final java.util.ArrayList yQj = 
		(java.util.ArrayList) yQ.get(j);
		// 		System.out.println("j="+j+"  yPj.size()="+yPj.size());
		
		boolean flag = false;

		for (int k = 0, limk = yPj.size(); k<limk ; k++) {
		    final int yPjk = ((Integer) yPj.get(k)).intValue();
		    if (yPjk>si) {
			flag=true;
			yPj.set(k,new Integer(si));
			si=yPjk;
			break;
			} 
		    
		    } 
		if (!flag) {
		    yPj.add(new Integer(si));
		    yQj.add(new Integer(i));
		    break jkloops;	
		    }
	
		} 
	    
	    }
	int dim=len;
	for (;--dim>=0;) {
	    if (((java.util.ArrayList) yP.get(dim)).size()!=0) break;
	    }
	
	final int[][][] y = new int[2][][];
	y[0] = new  int[++dim][];
	y[1] = new  int[dim][];
	
	for (int i=dim; --i>=0;) {
	final java.util.ArrayList yPI =((java.util.ArrayList) yP.get(i));
	final java.util.ArrayList yQI =((java.util.ArrayList) yQ.get(i));
	y[0][i] = new int[yPI.size()];
	y[1][i] = new int[yQI.size()];
	for (int  j = yPI.size(); --j>=0;) {
	    y[0][i][j]=	((Integer) yPI.get(j)).intValue();
	    y[1][i][j]=	((Integer) yQI.get(j)).intValue();
	    }
	}

	return y;	
	}
    /**
     * Gives the pair of <code>int[][]</code> Young tableaux associated with
     * the <code>int[]</code> permutation. Schensted's correspondance.
     *
     * @param inv an <code>int[]</code> value
     * @return an <code>int[][][]</code> value
     */
    static public int[] fromYoungTableaux(final int[][][] y){
	int len=0;
	final int dim = y[0].length;
	final java.util.ArrayList yP = new java.util.ArrayList(dim);
	
	for (int i=0; i<dim;i++) {
	    final int leni = y[0][i].length;
	    len+=leni;
	    final java.util.ArrayList yPI =new java.util.ArrayList(leni);

	    for (int  j = 0; j<leni;j++) {
	    yPI.add(new Integer(y[0][i][j]));
	    }
	    yP.add(yPI);


	    }

	final int[] result = new int[len];

	for (int  i = len; --i>=0;) {
	    
	    int j=0,k=0;

	    jkloops: // Find i in the second tableau.
	    for (int limj = yP.size();  j<limj ; j++) {
		k=0;
		for (int limk = y[1][j].length; k<limk ; k++) {
		    if (y[1][j][k]==i)  break jkloops;	
		    }
		}
	    // 	    System.out.println("i="+i+"  j="+j+"  k="+k);
	    
	    int value=0;
	    java.util.ArrayList yPj = (java.util.ArrayList) yP.get(j);

	    Integer valInt = (Integer) yPj.get(k);

	    value = valInt.intValue();

	    yPj.remove(k);

	    while(--j>=0){
	       yPj = (java.util.ArrayList) yP.get(j);
	       k=0;
	       final int limk=yPj.size();
	       if (limk>0) {
		   // Find the largest entry smaller than value.
		   for (; k < limk ; k++) {
		       if (((Integer)yPj.get(k)).intValue()>value) break;
		       } 
		   k--;
		   }
	   
	       final Integer valjk = (Integer)yPj.get(k);
	       yPj.set(k,valInt);
	       valInt=valjk;
	       value = valInt.intValue();
	       }
	    result[i]= value;
	    }
	

	return result;	
	}
    /**
     * Returns the <code>runs</code> of a given permutation, as
     * <code>int</code> arrays. It consists of increasing subsequences. I
     * don't know what it's good for.
     *
     * @param p an <code>int[]</code> value
     * @return The <code>int[][]</code> giving the array of runs.
     */
    static public int[][] runs(final int[] s) {
	final int n = s.length;
	final int[][] c = new int[n][];
	final int[] grow = new int[n]; 

	int number = 0;

	for (int i=0; i<n; i++) {
	    int len=0;
	    for (int j=i; j<n ; j++) {
		grow[len++]=s[j];
		}
	    c[number]=new int[len];
	    System.arraycopy(grow,0,c[number],0,len);
	    number++;
	    i+=len;
	    }
	
	final int[][] result = new int[number][];
	System.arraycopy(c,0,result,0,number);

	return result;
	}

    /**
     * Returns the {@link #cycles} view of the permutation as a {@link String}.
     *
     * @param s an <code>int[]</code> permutation
     * @return a <code>String</code> value
     */
    static public String cyclesToString(final int[] s) {
	return cyclesToString(cycles(s));
    }

    /**
     * Returns the {@link #cycles} view of the permutation as a {@link
     * String}. Doesn't check for validity, any <code>int[][]</code> will do.
     *
     * @param c an <code>int[][]</code>, supposedly the cycles of a
     * permutation.
     * @return a <code>String</code> value
     */
    static public String cyclesToString(final int[][] c) {

	StringBuffer sb=new StringBuffer().append('(');
	for (int j=0; j<c.length;j++) {
	    sb.append('(').append(' ');
	    for (int i=0; i<c[j].length; i++) {
		sb.append(c[j][i]+" ");
		}
	    sb.append(')');
	    } // end of for ()
	sb.append(')');

	return sb.toString();
	}

    /**
     * Hands back an <code>int[][]</code> coded by a string. () are
     * considered as separating the different arrays.  Creates a temporary
     * <code>int[]</code> value big enough to hold all the int's.  Describe
     * <code>stringToCycles</code> method here.
     *
     * @param s a <code>String</code> value
     * @return an <code>int[][]</code> value
     */
    static final public int[][] stringToCycles(final  String s){
	try{

	    StreamTokenizer st=new StreamTokenizer(new StringReader(s));

	    st.resetSyntax();      //sets all characters to be a seperate token

	    st.wordChars('0', '9');
	    st.wordChars('\u0000', '\u0020'); // Tab and space are letters.

	    st.whitespaceChars('(', ')'); // Parenthesis are separators.

	    int num = 0;
	    int maxLen = 0;
	    {// We actually do it twice to count
		StreamTokenizer stC=new StreamTokenizer(new StringReader(s));
	    stC.resetSyntax();
	    stC.wordChars('0', '9');
	    stC.wordChars('\u0000', '\u0020');
	    stC.whitespaceChars('(', ')');
	    while (stC.nextToken()!=stC.TT_EOF){
		    final int len = stC.sval.trim().length();
		if (len>maxLen) maxLen = len;   
		
		if(len>0)  num++; // There may be spaces between parenthesis.
	    }
	    }

	    {
		final	int[] basket = new int[maxLen];
		int[][] cycles = new int[num][];

		num=0;

		while (st.nextToken()!=st.TT_EOF) {
		    final String theS = st.sval.trim();
		    if(theS.length()>0) cycles[num++]=
					    stringToIntArray(theS, basket);

		}
		
		return cycles;
	    }

	}
	catch(java.io.IOException ex){ throw new Error(); }
    }

    /**
     * Hands back an integer array coded by a string. Space, tab, () ,;: {}
     * are considered as white spaces.  Creates a temporary
     * <code>int[]</code> value big enough to hold all the int's.

     *
     * @param s a <code>String</code> value
     * @return an <code>int[]</code> value
     * @see stringToIntArray(String,int[])
     */
    static public int[] stringToIntArray(final  String s) {
	return stringToIntArray(s, new int[s.length()]) ; // That's big enough.
    }
    /**
     * Hands back an integer array coded by a string. Space, tab, () ,;: {}
     * are considered as white spaces.
     *
     * @param s a <code>String</code> value
     * @param basket a temporary <code>int[]</code> value big enough to hold
     * all the int's.
     * @return an <code>int[]</code> value
     */
    static public int[] stringToIntArray(final  String s, final int[] basket) {
	try{

	    StreamTokenizer st=new StreamTokenizer(new StringReader(s));

	    st.resetSyntax();      //sets all characters to be a seperate token

	    st.whitespaceChars('\u0000', '\u002C'); // Unicode Tab, space, ',', ';', ':',  '(', ')'
	    st.whitespaceChars('\u003A', '\u003B'); // '{', '}' are white spaces.
	    st.whitespaceChars('{', '}'); 

	    st.parseNumbers(); //character combinations of 0 till 9
	    //are interpreted as one token, e.g.: 13420

	    int num = 0;
      
	    while (st.nextToken()!=st.TT_EOF) {

		basket[num++]=(int) st.nval;
	    }

	    final int[] theIntArray = new int[num];
	    System.arraycopy(basket,0,theIntArray,0,num);
	    return theIntArray;

	}
	catch(java.io.IOException ex){ throw new Error(); }
    }

    /**
     * Describe <code>toString</code> method here.
     *
     * @return a <code>String</code> value
     */
    static public String toString(final int[] s) {
	StringBuffer sb=new StringBuffer().append('(').append(' ');
	for (int i=0; i<s.length; i++) {
	    sb.append(s[i]+" ");
	    } // end of for ()

	sb.append(')');

	return sb.toString();
	} 
    }// Permutation
