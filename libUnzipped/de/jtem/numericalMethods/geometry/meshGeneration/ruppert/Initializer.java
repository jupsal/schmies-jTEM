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

package de.jtem.numericalMethods.geometry.meshGeneration.ruppert;

class Initializer implements java.io.Serializable 
{
  private static final long serialVersionUID = 1L;

  // rectangular domain [ with left lower corner x0,y0, right upper corner x1,y1 
  //   and xnum+1 points on the horizontal, ynum+1 points on the vertical edges ]
  public static final Ruppert init(double x0,double y0,double x1,double y1,int xnum,int ynum) {
    return init(x0,y0,x1,y1,xnum,ynum,new double[0],new int[0]);
  }

  // rectangular domain with circular holes
  // [ h={x_0,y_0,r_0 ... x_n,y_n,r_n} and ref={ref_0 ... ref_n} are interpreted as follows: "circle" i
  //   is centered at x_i,y_i, has radius r_i and is a regular polygon consisting of ref_i points ]
  public static final Ruppert init(double x0,double y0,double x1,double y1,int xnum,int ynum,
				   double[] h,int[] ref) {
    if (xnum<1 || ynum<1) 
      throw new IllegalArgumentException("arguments(s) corrupt");

    // test circles ?

    double[][] p=new double[1+h.length/3][];

    p[0]=new double[4*(xnum+ynum)];
    int c=0;
    double xplus=(x1-x0)/(double)xnum, yplus=(y1-y0)/(double)ynum;
    for (int i=0;i<xnum;i++) { p[0][c++]=x0+xplus*(double)i; p[0][c++]=y0; }
    for (int i=0;i<ynum;i++) { p[0][c++]=x1; p[0][c++]=y0+yplus*(double)i; }
    for (int i=0;i<xnum;i++) { p[0][c++]=x1-xplus*(double)i; p[0][c++]=y1; }
    for (int i=0;i<ynum;i++) { p[0][c++]=x0; p[0][c++]=y1-yplus*(double)i; }
  
    for (int i=0;i<h.length/3;i++) {
      double pplus=Math.PI*2/(double)(ref[i]);
      c=0;
      p[i+1]=new double[2*ref[i]];
      for (int j=0;j<ref[i];j++) {
	   p[i+1][c++]=h[3*i+2]*Math.cos(pplus*(double)j)+h[3*i];
	   p[i+1][c++]=h[3*i+2]*Math.sin(pplus*(double)j)+h[3*i+1];
      }
    }

    return new Ruppert(p);
  }
    
}






