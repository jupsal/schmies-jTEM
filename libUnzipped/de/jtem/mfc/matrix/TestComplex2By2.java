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

package de.jtem.mfc.matrix;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.vector.Complex2;

/**
 * Test class for {@link Complex2By2}. If the test fails, there's definitely
 * something wrong. If it doesn't fail, there may be something wrong.
 *
 * @author Boris Springborn
 */
public class TestComplex2By2 {

    final static double EPS = 1.0e-14;
    
    static private double maxDiff(Complex2 v, Complex2 w) {
        final double d1 = Math.abs(v.aRe - w.aRe);
        final double d2 = Math.abs(v.aIm - w.aIm);
        final double d3 = Math.abs(v.bRe - w.bRe);
        final double d4 = Math.abs(v.bIm - w.bIm);
        final double dd1 = d1 > d2 ? d1 : d2;
        final double dd2 = d3 > d4 ? d3 : d4;
        return dd1 > dd2 ? dd1 : dd2;
    }

    static private double maxDiff(Complex z, Complex w) {
        final double d1 = Math.abs(z.re - w.re);
        final double d2 = Math.abs(z.im - w.im);
        return d1 > d2 ? d1 : d2;
    }

    public static void main(String[] args) {
        
        final Complex2By2 m1 = new Complex2By2();
        final Complex2    v1 = new Complex2();
        final Complex2    v2 = new Complex2();
        final Complex2    v3 = new Complex2();
        final Complex2    v4 = new Complex2();
        final Complex     z1 = new Complex();
        final Complex     z2 = new Complex();
        final Complex     z3 = new Complex();
        final Complex     z4 = new Complex();
        
        /*
         * Tests related to eigenvalues/eigenvectors.
         */

        /* Assign m1 by eigenvectors and eigenvalues. */
        z1.assign(2.5,  -3.2);
        z2.assign(-1.2, 10.0);
        v1.assign(1.0, 2.0, 3.0, 4.0);
        v2.assign(-2.0, -1.0, 0.0, 1.0);
        m1.assignByEigenvectors(z1, v1, z2, v2);

        /* Check if they're really the eigenvectors and eigenvalues. */
        v3.assignTimes(m1, v1);
        v4.assignTimes(z1, v1);
        if (maxDiff(v3, v4) > EPS) {
            throw new RuntimeException();
        }
        v3.assignTimes(m1, v2);
        v4.assignTimes(z2, v2);
        if (maxDiff(v3, v4) > EPS) {
            throw new RuntimeException();
        }

        /* Check if the eigenvalues are correctly calculated. */
        m1.getEigenValues(z3, z4);
        if (maxDiff(z1, z3) > EPS) {
            throw new RuntimeException();
        }
        if (maxDiff(z2, z4) > EPS) {
            throw new RuntimeException();
        }
    }
}
