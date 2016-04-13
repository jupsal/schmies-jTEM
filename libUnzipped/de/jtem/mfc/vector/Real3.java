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

package de.jtem.mfc.vector;

import java.io.Serializable;

import de.jtem.mfc.geometry.ComplexProjective1;

/**
 * <p> This class represents a real 3-vector <code>(x, y, z)</code>.</p>
 *
 * <p> All methods avoid the creation of temporarily used objects in their
 * internal computations unless stated in their descriptions. </p>
 *
 * @author Boris Springborn
 */

public class Real3 implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    /**
     * Threashold for zero tests.
     *
     * Default is 1.0e-14.
     */
    public static final double EPSILON = 1.0e-14;

    /**
     * The x-component.
     */
    public double x;

    /**
     * The y-component.
     */
    public double y;

    /**
     * The zcomponent.
     */
    public double z;

    /**
     * Create the zero vector.
     */
    public Real3() {
        assign(0.0, 0.0, 0.0);
    }

    /**
     * Create a vector with given entries.
     */
    public Real3(final double theX,
                 final double theY,
                 final double theZ) {
	assign(theX, theY, theZ);
    }

    /**
     * Create a copy of the given vector.
     */
    public Real3(final Real3 v) {
	assign(v);
    }

    /**
     * Get the value of x.
     * @return value of x
     */
    public double getX() {
        return x;
    }

    /**
     * Set the value of x.
     * @param v  value to assign to x
     */
    public void setX(final double  v) {
        this.x = v;
    }

    /**
     * Get the value of y.
     * @return value of y
     */
    public double getY() {
        return y;
    }

    /**
     * Set the value of y.
     * @param v  value to assign to y
     */
    public void setY(final double  v) {
        this.y = v;
    }

    /**
     * Get the value of z.
     * @return value of z
     */
    public double getZ() {
        return z;
    }

    /**
     * Set the value of z.
     * @param v  value to assign to z
     */
    public void setZ(final double  v) {
        this.z = v;
    }

    /**
     * Create a copy of <code>this</code> and return it.
     * @return a copy of <code>this</code>
     */
    public Real3 copy() {
	return new Real3(this);
    }

    /**
     * Assign <code>this</code> with zero.
     */
    public void assignZero() {
	assign(0.0, 0.0, 0.0);
    }

    /**
     * Assign <code>this</code> with the given entries.
     */
    public void assign(final double newX,
                       final double newY,
                       final double newZ) {
	this.x = newX;
        this.y = newY;
	this.z = newZ;
    }

    /**
     * Assign <code>this</code> this with the given vector.
     */
    public void assign(final Real3 v) {
	assign(v.x, v.y, v.z);
    }

    /**
     * Add a vector.
     *
     * @param v the vector to add
     */
    public void assignPlus(final Real3 v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    /**
     * Returns the sum of this vector with <code>v</code>.
     *
     * @param v the vector to add
     */
    public Real3 plus(Real3 v) {
        return new Real3( x+v.x, y+v.y, z+v.z );
    }

    /**
     * Returns a multiple of this vector.
     *
     * @param t the factor to multiply with
     */
    public Real3 times(double t) {
        return new Real3( t*x, t*y, t*z );
    }

    /**
     * Returns a multiple of this vector.
     *
     * @param t the factor to divide by
     */
    public Real3 divide(double t) {
        return new Real3( x/t, y/t, z/t );
    }
    
    /**
     * Subtract a vector.
     *
     * @param v the vector to subtract
     */
    public void assignMinus(final Real3 v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    /**
     * Subtracts vector <code>v</code> from <code>this</code>.
     *
     * @param Difference of <code>this</code> and <code>v</code>.
     */
    public Real3 minus(final Real3 v) {
    	return new Real3(x-v.x, y-v.y, z-v.z );
    }
    /**
     * Multiply by a scalar.
     *
     * @param a the scalar
     */
    public void assignTimes(final double a) {
        x *= a;
        y *= a;
        z *= a;
    }

    /**
     * Divide by a scalar.
     *
     * @param a the scalar
     */
    public void assignDivide(final double a) {
        x /= a;
        y /= a;
        z /= a;
    }

    /**
     * Assign <code>this</code> with a linear combination of two vectors.
     *
     * @param a a <code>double</code>, the factor of <code>v</code>
     * @param v a {@link Real3}
     * @param b a <code>double</code>, the factor of <code>w</code>
     * @param w a {@link Real3}
     */
    public void assignLinearCombination(final double a, final Real3 v,
                                        final double b, final Real3 w) {

        /* Use dummy variables, because v or w migh be this. */
        final double dummyX = a * v.x + b * w.x;
        final double dummyY = a * v.y + b * w.y;
        final double dummyZ = a * v.z + b * w.z;

        assign(dummyX, dummyY, dummyZ);
    }

    /**
     * Assign <code>this</code> with cross product of two vectors.
     *
     * @param v first factor in cross product.
     * @param w second factor in cross product.
     */
    public void assignCrossProduct(final Real3 v, final Real3 w) {

        /* Use dummy variables, because v or w migh be this. */
        final double dummyX = v.y * w.z - v.z * w.y;
        final double dummyY = v.z * w.x - v.x * w.z;
        final double dummyZ = v.x * w.y - v.y * w.x;

        assign(dummyX, dummyY, dummyZ);
    }


    /**
     * Computes cross product of <code>this</code> with presribed vector.
     *
     * @param v first factor in cross product.
     */
    public Real3 crossProduct(final Real3 v) {
    	Real3 result = new Real3();
    	result.assignCrossProduct( this, v );
    	return result;
    }
    
    /**
     * Return dot product of two vectors.
     *
     * @param v a {@link Real3}
     * @param w a {@link Real3}
     * @return <code>v.x * w.x + v.y * w.y + v.z * w.z</code>
     */
    final public static double dotProduct(final Real3 v, final Real3 w) {
        return v.x * w.x + v.y * w.y + v.z * w.z;
    }

    /**
     * Returns dot product of this with vector <code>v</code>.
     */
    final public double dot( final Real3 v ) {
    	return Real3.dotProduct( this, v );
    }
    
    /**
     * Return the square length of <code>this</code>.
     *
     * @return {@link #dotProduct(Real3, Real3) dotProduct(this, this)}
     * @see #norm
     * @see #normalize()
     */
    public double normSquared() {
        return dotProduct(this, this);
    }

    /**
     * Return the length of <code>this</code>.
     *
     * @return <code>&sqrt;(x * x + y * y + z * z)</code>
     * @see #normSquared()
     * @see #normalize()
     */
    public double norm() {
        return Math.sqrt(normSquared());
    }

    /**
     * Return the squared distance of <code>this</code> and <code>v</code>.
     */
    public double distSquared( Real3 v) {
        return (this.x-v.x)*(this.x-v.x) + 
        	   (this.y-v.y)*(this.y-v.y) + 
        	   (this.z-v.z)*(this.z-v.z) ;
    }

    /**
     *Return the distance of <code>this</code> and <code>v</code>.
     */
    public double dist( Real3 v) {
        return Math.sqrt(distSquared( v ));
    }

    /**
     * <p>Normalize <code>this</code>.</p>
     *
     * <p>Divide <code>this</code> by {@link #norm()}.</p>
     *
     * @see #normSquared()
     * @see #norm()
     * @see #normalize()
     */
    public void normalize() {
        assignDivide(norm());
    }

    /**
     * <p>Project a point on the unit sphere stereographically from the north
     * pole to complex projective space.</p>
     *
     * <p><code>(0, 0, 1)</code> is projected to the point with homogeneous
     * coordinates <code>(1, 0)</code>.<br> <code>(0, 0, -1)</code> is
     * projected to the point with homogenous coordinates <code>(0,
     * 1)</code>.</p>
     *
     * @param cp the {@link ComplexProjective1} used to hand over the
     * result. Must not be null.
     *
     * @throws IllegalArgumentException if {@link #normSquared()} differs from
     * <code>1.0</code> by more than {@link #EPSILON}
     */
    public void projectTo(final ComplexProjective1 cp)
        throws IllegalArgumentException {

        /* Check if the point lies on the sphere. */
        if (Math.abs(normSquared() - 1.0) > EPSILON) {
            throw new IllegalArgumentException("Norm must be 1.0. " +
                                               "Norm is " + norm() + ". ");
        }

        /* Check if the point lies above or below the xy-plane. */
        if (z > 0.0) {

            /* We might be close to the north pole. This would be bad,
             * numerically. So perform a 180 degree rotation around the
             * x-axis, project that point, and swap the entries of the
             * result.
             */

            /* Dummy variables */
            final double dummyRe;
            final double dummyIm;

            /* Rotate by 180 degree around the x-axis. */
            y = -y;
            z = -z;

            /* Project that point. Now, z < -EPSILON. */
            belowXYProjectTo(cp);

            /* Rotate back. */
            y = -y;
            z = -z;

            /* Swap the entries of the result.
             *
             * TODO: Maybe Complex2 should have a final method that swaps the
             * entries. This might speed up dereferencing the fields of cp,
             * because dereferencing fiels of other objects might be slower
             * than dereferencing fields of the same object. But actually, I
             * don't know if that is true or if it matters. */
            dummyRe = cp.aRe;
            dummyIm = cp.aIm;
            cp.aRe = cp.bRe;
            cp.aIm = cp.bIm;
            cp.bRe = dummyRe;
            cp.bIm = dummyIm;

        } else {

            /* Point lies below xy-plane so we use */
            belowXYProjectTo(cp);
        }
    }

    /* Performs the actual calculation. Is only used if point lies below the
     * xy-plane.
     */
    private final void belowXYProjectTo(ComplexProjective1 cp) {
//         final double oneMinusZ = 1.0 - z;
//         final double normalizingFactor = Math.sqrt(2.0 * oneMinusZ);
        final double bRe = 1.0 - z; // oneMinusZ / normalizingFactor;
        final double bIm = 0.0;
        final double aRe = x; // / normalizingFactor;
        final double aIm = y; // / normalizingFactor;
        cp.assign(aRe, aIm, bRe, bIm);
    }

    /**
     * Return a <code>String</code> representation of <code>this</code>.
     *
     * @return a <code>String</code> representation of <code>this</code>
     */
    public String toString() {
        return new String("(" + x + ", " + y + ", " + z + ")");
    }

}
