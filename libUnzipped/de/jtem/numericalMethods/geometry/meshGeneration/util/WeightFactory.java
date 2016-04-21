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

package de.jtem.numericalMethods.geometry.meshGeneration.util;


public class WeightFactory {

  WeightFactory() {
  }

  final double areaOfTriangle2D( final double [] domain, final int [] index, int face ) {
    final int p0 = 2*index[3*face];
    final int p1 = 2*index[3*face+1];
    final int p2 = 2*index[3*face+2];

    final double x0 = domain[p0  ];
    final double y0 = domain[p0+1];

    final double dx1 = domain[p1  ] - x0;
    final double dy1 = domain[p1+1] - y0;

    final double dx2 = domain[p2  ] - x0;
    final double dy2 = domain[p2+1] - y0;

    final double dot11 = dx1*dx1+dy1*dy1;
    final double dot12 = dx1*dx2+dy1*dy2;
    final double dot22 = dx2*dx2+dy2*dy2;

    return 0.5*Math.sqrt( dot11 * dot22 - dot12 * dot12 );
  }

  final double areaOfTriangle3D( final double [] point,  final int [] index, int face ) {
    final int p0 = 3*index[3*face];
    final int p1 = 3*index[3*face+1];
    final int p2 = 3*index[3*face+2];

    final double x0 = point[p0] ;
    final double y0 = point[p0+1];
    final double z0 = point[p0+2];

    final double dx1 = point[p1]   - x0;
    final double dy1 = point[p1+1] - y0;
    final double dz1 = point[p1+2] - z0;

    final double dx2 = point[p2]   - x0;
    final double dy2 = point[p2+1] - y0;
    final double dz2 = point[p2+2] - z0;

    final double dot11 = dx1*dx1+dy1*dy1+dz1*dz1;
    final double dot12 = dx1*dx2+dy1*dy2+dz1*dz2;
    final double dot22 = dx2*dx2+dy2*dy2+dz2*dz2;

    return 0.5*Math.sqrt( dot11 * dot22 - dot12 * dot12 );
  }

  /**
   * Creates the weights at the vertices as the ratio of the sum of the area
   * of the triangles in 3 space agacent to a vertex and the sum of those in the domain.
   * @param domain coords of domain (length: 2 * number of vertices)
   * @param point coords of immersion (length: 3 * number of vertices)
   * @param index indices of points which forms a triangle (length 3 * number
   * of triangles)
   * @return weight array with length equal to the number of vertices
   */
  public double [] createByAreaRatio( final double [] point, final double [] domain,
                                      final int [] index ) {
    double [] weight = new double[ point.length / 3 ];
    computeByAreaRatio( point, domain, index, weight, new double[ point.length / 3] );
    return weight;
  }

  /**
   * Computes the weights at the vertices as the ratio of the sum of the area
   * of the triangles in 3 space agacent to a vertex and the one in the domain.
   * @param domain coords of domain (length: 2 * number of vertices)
   * @param point coords of immersion (length: 3 * number of vertices)
   * @param index indices of points which forms a triangle (length 3 * number
   * of triangles)
   * @param weight result on output (length: number of vertices)
   * @param tmp for temp. used data (length: number of vertices)
   */
  public void computeByAreaRatio( final double [] point, final double [] domain,
                                  final int [] index,
                                  final double [] weight, final double [] tmp ) {

    final int nop = domain.length / 3;
    final int noe = index.length / 3;

    for( int i=0; i<weight.length; i++ )
      weight[i] = tmp[i] = 0;

    for( int i=0, I=0; i<noe; i++ ) {

      final double area2 = areaOfTriangle2D( domain, index, i );
      final double area3 = areaOfTriangle3D(  point, index, i );

      for( int j=0; j<3; j++ ) {
        weight[ index[I] ] += area3;
        tmp   [ index[I] ] += area2;
      }
    }

    for( int i=0; i<weight.length; i++ )
      weight[i] /= tmp[i];
  }
}
