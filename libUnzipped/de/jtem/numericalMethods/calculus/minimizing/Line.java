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

package de.jtem.numericalMethods.calculus.minimizing;

/**
 * This class represents a line in the n-dimensional space with first point and direction.
 * @author Markus Schmies, Vitali Lieder
 * @version 1.0
 */
public class Line implements java.io.Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * first point on line
     */
    final double [] point;
    /**
     * direction at the line
     */
    final double [] direction;
    /**
     * other point in the line. we need it for some package internal calculation so let it be.
     */
    final double [] otherPoint;
    /**
     * dimension of line.
     */
    final public int n;

    /**
     * Default Constructor with point and the direction in space.
     * @param point initial point of line
     * @param direction
     * @throw IllegalArgumentException if point and direction are not the same dimension.
     */
    public Line( double [] point, double [] direction ) {

        n = point.length;

        if( n != direction.length )
            throw new IllegalArgumentException
                (" dimension of direction and point do not coincide ");

        this.point     = point;
        this.direction = direction;

        otherPoint = new double[n];
    }

    /**
     * Calculates point+t*direction and write it in pointAtT.
     * @throw IllegalArgumentException if pointAtT are not the same dimension.
     * @param t scalar
     * @param pointAtT new point at the line.
     */
    final public void getPoint( final double t, final double [] pointAtT ) {
      if( n != pointAtT.length )
        throw new IllegalArgumentException
            (" dimension of pointAtT does not coincide ");
      for( int i = 0; i < n; i++)
        pointAtT[i] = point[i] + t * direction[i];
    }

}
