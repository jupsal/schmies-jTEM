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

/** To triangulate a domain, we start with discrete line segments representing the boundaries.
  * The {@link Delaunay} algorithm yields a Delaunay triangulation for those boundary points, 
  * but we certainly want to add interior points. They should be arranged in a way that yields
  * a <q>good</q> triangulation, i.e. with prescribed maximal area and minimal angle for all triangles.
  * <p>
  * The algorithm introducing these points follows 
  * <a href="http://www-2.cs.cmu.edu/~quake/tripaper/triangle3.html">Ruppert's paper</a>.
  * <ul>
  * <li>Each <i>segment</i> is split in half if any point lies <i>inside</i> its <i>diametral</i> circle.
  * <li>In a second step, each <q>bad</q> triangle is split by inserting a
  * point at its <i>circumcenter</i> and maintaining the Delaunay property. 
  * <li>The circumcenter is <i>not</i> inserted if it produces bad triangles that include segments. 
  * In those cases, the according segments are split in half.
  * <li>New points are <i>not</i> inserted if they are too close to the given ones. 
  * This way we can maintain relatively coarse triangulations.
  * <li> 
  * Marked segments which are visible from the new point will be
  * included into the new triangulation. This way, segments cannot be replaced. 
  * </ul>
  * Triangles are <q>bad</q> if their area is larger or one of their angles is
  * smaller than prescribed.
  *  
  * Note: By default, the maximal number of triangles are set to 1000. 
  * Use {@link #setMaximalNumberOfTriangles(int i)} to alter this number.
  *
  * @see Delaunay
  */
public class Ruppert extends Delaunay
{
  private static final long serialVersionUID = 1L;

  double[] cosAngle = new double[3*maximalNumberOfFaces];
  double[] area = new double[maximalNumberOfFaces];

  double cosLeastAngle; // =0.90631;   cos30=0.866  cos25=0.90631  cos20=0.94
  double largestArea;  // 0.1;

  int maximalNumberOfTriangles = 100000;

  int numberOfBadTriangles=0;                    
  int maxNob=100000;
  int[] badTriangles=new int[maxNob];

  int numberOfRejectedBadTriangles;                     // number of bad triangles not refined in current iteration
  int[] rejectedBadTriangles=new int[maxNob];

  boolean refineFlag=false;

  double sqrLeastDistance=0;


  private boolean debug = false;


  /**
	 * A Delaunay triangulation for the given boundaries is created.
	 * <p>
	 * You may change the angle and/or area constraint and
	 * use the {@link refine} method to refine the triangulation according to these constraints.
	 * <br>
	 * The constraints can be altered using {@link setAngleConstraint(double a)} and/or 
	 * {@link setAngleConstraint(double a)}.
	 * <p>
	 * The array <code>p</code> is given as <code>x<sub>1,1</sub>,y<sub>1,1</sub>,
	 * x<sub>1,2</sub>,y<sub>1,2</sub>,...</code>
	 * <p>
	 * <b>Note:</b>
	 * <ul>
	 * <li> All boundaries are connected, e.g. <code>x<sub>1,n</sub>,y<sub>1,n</sub> 
	 * with <code>x<sub>1,1</sub>,y<sub>1,1</sub>
	 * <li> The first array <code>p[0]</code> must contain the exterior boundary
	 * </ul>
	 * <p>
	 * @param p the array given as <code>x<sub>1,1</sub>,y<sub>1,1</sub>,
	 * x<sub>1,2</sub>,y<sub>1,2</sub>,...</code> containing the boundary points
	 */    
  public Ruppert(final double[][] p)
  {
    super(p);

    computeAngleAndArea();
    refineFlag=true;

    if( debug ) System.out.println("checking segments");

    checkSegments();
    
    cosLeastAngle = cosLeastAngle();
    largestArea = largestArea() *1.01; 
  }


  /** You can prescribe the maximal number of triangles generated by the refinement process.
   */
  public void setMaximalNumberOfTriangles(int i) { maximalNumberOfTriangles = i; }

  /**
   * @return The prescribed maximal number of triangles generated by the refinement process.
   */
  public int getMaximalNumberOfTriangles() { return maximalNumberOfTriangles; }

  /**
   * You can prescribe a minimal angle for all triangles. 
   * It is, however, not reasonable to prescribe values much larger than 30 degrees.
   * <p>
   * <b>Note: </b>
   * This method does NOT trigger the refinement process. 
   * Call the {@link refine()} method to start this process.
   * @param a the prescribed angle constraint [unit: degree]
   */
  public void setAngleConstraint(final double a) { cosLeastAngle = Math.cos(a*Math.PI/180); }

  /**
   * @return The prescribed angle constraint [unit: degree] for the refinement process.
   */
  public double getAngleConstraint() { return Math.acos(cosLeastAngle)*180/Math.PI; }

  /**
   * You can prescribe a maximal area for all triangles. 
   * <p>
   * <b>Note: </b>
   * This method does NOT trigger the refinement process. 
   * Call the {@link refine()} method to start this process.
   * @param a the prescribed area constraint
   */
  public void setAreaConstraint(final double a) { largestArea = a; }

  /**
   * @return The prescribed area constraint for the refinement process.
   */
  public double getAreaConstraint() { return largestArea; }

  public void refine( int [] triangleList ) {
		  for( int i=0 ; i<triangleList.length; i++ )  {
			 refine(triangleList[i]);
		  }   
  }
	  
  
  /** Refines the current traingulation according to the given area/angle constraints.
   */
  public void refine()
  {
    sqrLeastDistance = 0.25 * Math.min( leastEdgeSqr(), largestArea );

    checkForBadTriangles();
    int oldNob=numberOfBadTriangles+1;

    if( debug ) System.out.println("start refining");

    while ( numberOfBadTriangles<oldNob )
    {
      oldNob=numberOfBadTriangles;
	
	  // We'll check the bad guys successively and do something about it IF WE CAN
	  for( int i=0 ; i<=numberOfBadTriangles; i++ )  {
	  	if( numberOfFaces > maximalNumberOfTriangles ) break;
		final int j=badTriangles[i];
		if (bad(j) ) refine(j);              // The triangle may have become a good triangle by now !
	  }   
 
      // So what to do with those triangles rejected on first glance ?
      // We' ll stack them on the bad pile again !
 	  numberOfBadTriangles=0;
      for (int k=0;k<numberOfRejectedBadTriangles;k++) {
      	 final int j=rejectedBadTriangles[k];
      	 if ( bad(j) )
      	 	badTriangles[numberOfBadTriangles++]=j; 
      }
    }

    if( numberOfFaces>=maximalNumberOfTriangles )
      if( debug ) System.out.println("reached maximal number of triangles");

    result();
  }

  /** Segments are edges that belong to the boundary.
   */
  public boolean[] getSegments()
  {
    boolean[] r=new boolean[3*numberOfFaces];
    System.arraycopy(segment,0,r,0,3*numberOfFaces);
    return r;
  }

  /**
   * Each triangle is referred to by an index as discribed in {@link Delaunay#getIndices}.
   * Accordingly, the index referencing the neighbors contains the indices of the neighboring triangles:
   * <code>... n<sub>3i</sub>, n<sub>3i+1</sub>, n<sub>3i+2</sub> ...</code> refer to the neighbors 
   * adjacent to edges 1 through 3 of triangle number <code>i</code>.
   * @return The index list of all faces' neighbors.
   */
  public int[] getNeighbors()
  {
    int[] r=new int[3*numberOfFaces];
    System.arraycopy(neighbor,0,r,0,3*numberOfFaces);
    return r;
  }

  // An edge is illegal if either of the adjacent triangles has the following property:
  // Inside its circumcircle lies an point of the triangulation
  boolean edgeIllegal(final int f,final int e)
  {
    if (refineFlag) {
      final int fp=3*f+e;
      if (segment[fp]) return false; // We have to keep the edges on the boundary !
      final int n=neighbor[fp];
      if (n<0) return false;             // There is no neighboring triangle !
      int np=2*face[getNeighborPtr(n,f)];
      return pointInCircle(point[np++],point[np],f); 
    }
    else return super.edgeIllegal(f,e);
  }

  void checkSegments() { for (int i=0;i<numberOfFaces;i++) checkSegments(i); }

  void checkSegments(final int f) { for (int i=0;i<3;i++) checkSegment(f,i); }

  // Is point x,y inside a circle with the edge e's diameter ?
  // (Which means the same as: Is it too close to the edge ?)
  boolean encroached(final int f,final int e,final double x,final double y)
  {
    final int fp=3*f;
    int p1=2*face[fp+(e+1)%3], p2=2*face[fp+(e+2)%3];
    final double x1=point[p1++],y1=point[p1],x2=point[p2++],y2=point[p2];
    circleX=(x1+x2)*0.5;
    circleY=(y1+y2)*0.5;
    final double xr=(x2-x1)*0.5, yr=(y2-y1)*0.5, xd=(x-circleX), yd=(y-circleY);
    return ( xd*xd+yd*yd + 1e-4  < xr*xr+yr*yr);
  }

  int encroached(final int f,final double x,final double y)
  {
    int fp=3*f;
    for (int i=0;i<3;i++,fp++)
      if (segment[fp]) { if (encroached(f,i,x,y)) return i; }
    return -1;
  }

  boolean encroached(final int f,final int e)
  {
    int fp=3*f+e;
    if (segment[fp]) {
      int p0=2*face[fp];
      return encroached(f,e,point[p0++],point[p0]);
    }
    else return false;
  }

  // is the point currently inserted too close to the boundary ?
  // Then we'll split the boundary segment ! 
  void checkSegment(final int f,final int e)
  {
	 
    if (encroached(f,e)) {
      int c=numberOfFaces;

      addPoint(circleX,circleY,f,e);

      checkSegment(f,e);
      checkSegment(c,0);
    }
  }

  // What does "legalize" mean in this context ?
  // Well: All new edges must be legal !
  void legalizeNewFaces()
  {
    if (refineFlag) {
      int f0=newFace[0], f1=newFace[1], f2=newFace[2], f3=newFace[3];

      super.legalizeNewFaces();

      checkSegments(f0);
      checkSegments(f1);
      if (f2>=0) checkSegments(f2);
      if (f3>=0) checkSegments(f3);
    }
    else super.legalizeNewFaces();
  }


  // Flipping the edge requires of course the adjustment of other information such as area, angle, etc.
  int flipEdge(final int f,final int e)
  {
    if (refineFlag) {

      final int fp=3*f, n=neighbor[fp+e], f2=fp+(e+2)%3;
      final int np=3*n, ne=getNeighbor(n,f), n2=np+(ne+2)%3;

      boolean sf=segment[f2];
      boolean sn=segment[n2];

      int r=super.flipEdge(f,e);

      segment[fp+e]=sn;
      segment[f2]=false;
      segment[np+ne]=sf;
      segment[n2]=false;

      computeAngleAndArea(f);
      computeAngleAndArea(n);

      return r;
    }
    else return super.flipEdge(f,e);
  }

  void splitTriangle(final int f)
  {
    if (refineFlag) {

      int fp=3*f,c=numberOfFaces,cf=3*numberOfFaces;
      boolean s0=segment[fp], s2=segment[fp+2];

      super.splitTriangle(f);

      segment[fp++]=false; segment[++fp]=false;
      segment[cf++]=false; segment[cf++]=false; segment[cf++]=s2;
      segment[cf++]=false; segment[cf++]=false; segment[cf]=s0;

      computeAngleAndArea(f);
      computeAngleAndArea(c++);
      computeAngleAndArea(c);
    }
    else super.splitTriangle(f);
  }

  void splitEdge(final int f,final int e)
  {
    if (refineFlag) {

      final int fp=3*f, f2=fp+(e+2)%3, n=neighbor[fp+e], np=3*n;
      int c=numberOfFaces, cf=3*numberOfFaces, ne=0, n1=0;
      final boolean sf=segment[fp+e], sf2=segment[f2];
      boolean sn1=false;

      if (n>=0) {
	    ne=getNeighbor(n,f);
	    n1=np+(ne+1)%3;
	    sn1=segment[n1];
      }

      super.splitEdge(f,e);

      segment[f2]=false;
      segment[cf++]=sf;
      segment[cf++]=false;
      segment[cf]=sf2;

      computeAngleAndArea(f);
      computeAngleAndArea(c);

      if (n>=0) {
	   segment[n1]=false;
	   segment[++cf]=false;
	   segment[++cf]=sf;
	   segment[++cf]=sn1;

	   computeAngleAndArea(n);
	   computeAngleAndArea(c+1);
      }
    }
    else super.splitEdge(f,e);
  }


	// This is the (literally) decisive method:
	// HOW to refine the triangle containing the new point ?!?
  	void refine(final int f)
  	{
    	computeCircumCircle(f);                            // The circumcenter (i.e. the new point) is stored in  circleX,circleY

    	final int i=findTriangle(circleX,circleY);     // The new point is located inside triangle  i

		if (i>=0)                                                      // The point lies INSIDE a triangle (.... and not on an edge)
		{
			final double x=circleX, y=circleY;
       		int j=encroached(i,x,y);
       		
       		if (j<0)                                                    // The new point is NOT too close to the boundary
       		{
	    		if (!pointIsToClose(x,y,i))                   // The new point is NOT too close to any othe point
	    		{
	       			if (pointOnEdge<0) addPoint(x,y,i);             // The generic case: We add the new point
	    			else checkNeighborForSplit(f,i,x,y);            // The new point is too close to an edge
				}                                                                      //      so we check the neighboring triangle as well
				else rejectBadTriangle(f);
      		}
      		else                                                       // The new point is too close to the boundary
      		{
				if (!pointIsToClose(circleX,circleY,i))             
	  				addPoint(circleX,circleY,i,j);                        // The generic case again: We add the new point
				else rejectBadTriangle(f);
      		}
    	}
    	else                                                           // The point lies on an existing edge
      	{
      		int e=largestEdge(f);
      		if (segment[3*f+e]) 
      		{
				computeMean(f,e);
				if (!pointIsToClose(circleX,circleY,f)) 
					addPoint(circleX,circleY,f,e);                      // The generic case again: We add the new point
				else rejectBadTriangle(f);
      		}
      		else rejectBadTriangle(f);
    	}
	}


    // The refine-method has found the new point is too close to an edge,
    // so we check the neighboring triangle as well
	void checkNeighborForSplit(final int f,final int i,final double x,final double y)
	{
		final int n=neighbor[3*i+pointOnEdge];
		if (n>=0)                                                            // The neighboring triangle must be checked 
		{                
			if (!pointIsTooClose(x,y,n,getNeighbor(n,i))) 
			{
				final int k=encroached(n,x,y);
				if (k<0) addPoint(x,y,i,pointOnEdge);              // The generic case: We add the new point
				else 
				{
					if (pointIsToClose(circleX,circleY,n))
						addPoint(circleX,circleY,n,k);           // We split the boundary segment in neighboring triangle

					else rejectBadTriangle(f);
				}
			}
			else rejectBadTriangle(f);                     
		}
		else addPoint(x,y,i,pointOnEdge);                   // There is no neighboring triangle
	}                                                                           // Hence, we can insert the point on the edge
	

  	// Is point x,y to close to an existing point ? 
  	boolean pointIsToClose(final double x,final double y,final int f)
  	{
    	for (int i=0;i<3;i++) 
    		if (pointIsTooClose(x,y,f,i)) 
    			return true;
    	return false;
  	}

  	// Is point x,y to close to an existing point ? 
  	boolean pointIsTooClose(final double x,final double y,final int f,final int e)
  	{
    	// weighted area constraint support
    	return distanceSqr( x, y, f, e ) < sqrLeastDistance;
  	}


  	void computeMean(final int f,final int e)         // ... of an edge
  	{
    	final int fp=3*f;
    	int p1=2*face[fp+(e+1)%3], p2=2*face[fp+(e+2)%3];
    	circleX=0.5*(point[p1++]+point[p2++]);
    	circleY=0.5*(point[p1]+point[p2]);
  	}

  	int largestEdge(final int f)
  	{
    	final int fp=3*f;
    	if (cosAngle[fp]<cosAngle[fp+1]) 
    	{
      		if (cosAngle[fp]<cosAngle[fp+2]) return 0;
      		else return 2;
    	} else {
      		if (cosAngle[fp+1]<cosAngle[fp+2]) return 1;
      		else return 2;
    	}
  	}


  	void computeAngleAndArea(final int f)
  	{
    	int fp=3*f, p0=2*face[fp++], p1=2*face[fp++], p2=2*face[fp];
    	final double x0=point[p0++], y0=point[p0], x1=point[p1++], y1=point[p1], x2=point[p2++], y2=point[p2];
    	final double e0x=x1-x0, e0y=y1-y0, e1x=x2-x1, e1y=y2-y1, e2x=x0-x2, e2y=y0-y2;
    	final double el0=Math.sqrt(e0x*e0x+e0y*e0y);
    	final double el1=Math.sqrt(e1x*e1x+e1y*e1y);
    	final double el2=Math.sqrt(e2x*e2x+e2y*e2y);
    	fp=3*f;
    	cosAngle[fp++]=-(e0x*e2x+e0y*e2y)/(el0*el2);
    	cosAngle[fp++]=-(e0x*e1x+e0y*e1y)/(el0*el1);
    	cosAngle[fp]=-(e1x*e2x+e1y*e2y)/(el1*el2);

   		// weighted area contraint support
    	area[f] = area(f);

    	checkForBadTriangles(f);
  	}

  	void computeAngleAndArea() { for (int i=0;i<numberOfFaces;i++) computeAngleAndArea(i); }

	// A triangle is bad if it it has a small angle or a large area
  	public boolean bad(final int f) {
  		return bad( f, cosLeastAngle, largestArea );
  	}
  	
  	public boolean bad(final int f, final double cosLeastAngle, final double largestArea )
  	{
    	int fp=3*f;
    	if (cosAngle[fp++]>cosLeastAngle) return true;
    	if (cosAngle[fp++]>cosLeastAngle) return true;
    	if (cosAngle[fp]>cosLeastAngle) return true;
    	if (area[f]>largestArea) return true;
    	return false;
  	}

  	void checkForBadTriangles() { numberOfBadTriangles=0; for (int i=0;i<numberOfFaces;i++) checkForBadTriangles(i); }

  	void checkForBadTriangles(final int f) { if (bad(f)) addToBadTriangles(f); }

  	void addToBadTriangles(final int f)
  	{
    	if (numberOfBadTriangles==maxNob) 
    	{ 
    		badTriangles=doubleSize(badTriangles); 
    		rejectedBadTriangles=doubleSize(rejectedBadTriangles); 
    		maxNob*=2;
    	}
    	badTriangles[numberOfBadTriangles++]=f;
  	}


  	void rejectBadTriangle(final int f)
  	{
    	if (numberOfRejectedBadTriangles==maxNob) 
    	{ 
    		badTriangles=doubleSize(badTriangles); 
    		rejectedBadTriangles=doubleSize(rejectedBadTriangles); 
    		maxNob*=2; 
    	}
    	rejectedBadTriangles[numberOfRejectedBadTriangles++]=f;
  	}


  	void checkFaceArray()
  	{
    	super.checkFaceArray();

    	if (area==null) area=new double[maximalNumberOfFaces];
    	if (cosAngle==null) cosAngle=new double[3*maximalNumberOfFaces];

    	if (maximalNumberOfFaces > area.length ) 
    	{
      		area=doubleSize(area);
      		cosAngle=doubleSize(cosAngle);
    	}
  	}

  	double leastEdgeSqr()
  	{
    	double r=4*xyBound*xyBound;
    	for (int i=0,c=0;i<numberOfFaces;i++,c+=3)
      		for (int j=0;j<3;j++)
				if (neighbor[c+j]<i) 
				{
 					// weighted area constraint support
	  				int p1 = face[c+(j+1)%3];
         			 int p2 = face[c+(j+2)%3];
          			double d = distanceSqr( p1, p2 );
          			if (r > d) { r = d; }
				}
		return r;
  	}

  	double cosLeastAngle()
  	{
    	double r=0;
    	final int imax=3*numberOfFaces;
    	for (int i=0;i<imax;i++) if (cosAngle[i]>r) r=cosAngle[i];
    	return r;
  	}

  	double largestArea() 
  	{
   		double r=0;
    	for (int i=0;i<numberOfFaces;i++) if (area[i]>r) r=area[i];
    	return r;
  	}

  	void result()
  	{
    	if( debug )
    	System.out.println(	"least angle = "+(180*Math.acos(cosLeastAngle())/Math.PI)+
										"  largest area = "+(largestArea())+
										"  number of faces = "+numberOfFaces);
  	}

  	public boolean isDebug() {
    	return debug;
  	}
  	
	public void setDebug(boolean debug) {
    	this.debug = debug;
  	}


  	// -----------------------------------------------------------
  	// Weight constraint support

  	/** Get the array of weight constraints indexed to the points.
      	  Returns null if no weight constraints have been set with
      	  setWeights().
       */
  	public double [] getWeight() {
   		if (this.weight == null) { return null; }
    	double [] w = new double [numberOfPoints];
    	System.arraycopy( weight, 0, w, 0, numberOfPoints );
    	return w;
  	}
  	
  	/** Set an array of weight constraints indexed to the points.
    	  The length of w[] must be equal to getNumPoints().
   		  w=null disables weight constraining (default).
	   */
	public void setWeight( double [] w ) {
    	if (w == null) { weight = null; return; }
    	if (w.length != numberOfPoints) 
      		throw new IllegalArgumentException("number of weights != number of points" );
    	weight = new double [maximalNumberOfPoints];
    	System.arraycopy( w, 0, weight, 0, numberOfPoints );
  	}



  	/** Overload of Delaunay.addPoint() for weighted area support.
     	  Interpolates the weight at the new point from the weights
          at the vertices of the encompassing triangle.
 	   */
  	void addPoint(final double x,final double y,final int f) {
    	if (weight != null) addWeight(x, y, f, numberOfPoints); 
    	super.addPoint(x, y, f);
  	}
  	
  	/** Overload of Delaunay.addPoint() for weighted area support.
     	  Interpolates the weight at the new point from the weights
     	  at the vertices of the encompassing triangle.
	   */
  	void addPoint(final double x,final double y,final int f,final int e) {
    	if (weight != null) addWeight(x, y, f, numberOfPoints); 
    	super.addPoint(x, y, f, e);
  	}
  	
  	/** Overloads of Delaunay.checkPointArray() for weighted area support.
     	  Reallocates the weight[] array.
     	  Uses knowledge of Delaunay.checkPointArray().
  	   */
	void checkPointArray() {
    	int oldMaxNop = maximalNumberOfPoints;
    	super.checkPointArray();
    	if (weight != null) 
      		if ( oldMaxNop != maximalNumberOfPoints ) 
      			weight = doubleSize(weight);
	}
	
  	/** Computes the weight at new point (x,y) with index n before it
     	  is added to the triangulation via Delaunay.addPoint().
     	  Interpolates the weight at the new point from the weights
     	  at the vertices of the encompassing triangle.
       */
	void addWeight( double x, double y, int f, int n ) {
    	if (weight != null) 
      		weight[n] = interpolatedWeight(x, y, f);
	}
 
 	/** Returns the weight at (x,y) interpolated from the weights
    	 at the vertices of the triangle indexed by f.
    	 The interpolation is computed by baricentric coordinates.
       */
  	double interpolatedWeight( double x, double y, int f ) {
    	int p0 = face[3*f];
    	int p1 = face[3*f+1];
    	int p2 = face[3*f+2];

    	double x0 = point[2*p0];
    	double y0 = point[2*p0+1];
    	double x1 = point[2*p1];
    	double y1 = point[2*p1+1];
    	double x2 = point[2*p2];
    	double y2 = point[2*p2+1];
    	double d = det(x0,y0,x1,y1,x2,y2);

    	return ( 	weight[p0] * det(x,y,x1,y1,x2,y2) +
						weight[p1] * det(x0,y0,x,y,x2,y2) +
						weight[p2] * det(x0,y0,x1,y1,x,y) ) / d;
  	}
  	
	/** Returns the square of the weighted distance beweeen the points indexed by p0, p1.
  	   */
	double distanceSqr( int p0, int p1 ) {
    	double x0 = point[2*p0];
    	double y0 = point[2*p0+1];
    	double x1 = point[2*p1];
    	double y1 = point[2*p1+1];
    	double d = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
    	if (weight != null) 
      		d *= 0.5 * ( weight[p0] + weight[p1] );
		return d;
	}
	
  	/** Returns the square of the weighted distance beweeen the point (x,y)
     	 and point e in face f.
     	If weight==null, returns the Euclidean square of the distance.
  	  */
	double distanceSqr( double x, double y, int f, int e ) {
    	int p = face[3*f+e];
    	double x0 = point[2*p];
    	double y0 = point[2*p+1];
    	double d = (x0 - x)*(x0 - x) + (y0 - y)*(y0 - y);
    	if (weight != null)
			d *= 0.5 * ( weight[p] + interpolatedWeight( x, y, f ) );
        return d;
	}
	
  	/** Returns twice the weighted area of the triangle indexed by f.
     	  If weight==null, returns twice the Euclidean area.
   	   */
	protected double area( int f ) {
    	int p0 = face[3*f];
    	int p1 = face[3*f+1];
    	int p2 = face[3*f+2];
    	double x0 = point[2*p0];
    	double y0 = point[2*p0+1];
   		double x1 = point[2*p1];
    	double y1 = point[2*p1+1];
    	double x2 = point[2*p2];
    	double y2 = point[2*p2+1];
    	double a = det( x0, y0, x1, y1, x2, y2 ) / 2.0;
     	if (weight != null) 
    		a *= ( 	weight[p0]*weight[p0] +
             			weight[p1]*weight[p1] +
            			weight[p2]*weight[p2] ) / 3.0;
		return a;
  	}
  	
	/** Returns twice the area of the triangle with vertices (x0,y0), (x1,y1), (x2,y2).
		  Equal to Det[{{x0,y0,1},{x1,y1,1},{x2,y2,1}}].
  		*/
	static double det(
    	double x0, double y0,
    	double x1, double y1,
    	double x2, double y2 ) {
    	return (y1 - y0)*(x0 - x2) - (x1 - x0)*(y0 - y2);
  	}

  	double [] weight = null;

}
