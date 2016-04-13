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

package de.jtem.mfc.set;

import java.beans.PropertyEditorSupport;

/**
 * Text based editor for Interval. 
 */
public final class IntervalEditor extends PropertyEditorSupport {
    Interval value;
    
    public IntervalEditor() {
	value = new Interval(1);
    }

    public String getAsText() {
	return value.toString();
    }

    public synchronized void setAsText( String s ) {	
	setValue( new Interval( s ) );
    }

    public Object getValue() {
	return new Interval( value );
    }

    public void setValue( Object o ) {
	Interval newVal = (Interval)o;

	if(!value.equals( newVal )) {
	    value.assign( newVal );
	    firePropertyChange();
	}
    }
    public String getJavaInitializationString()
    {
        double [] min=value.min;
        double [] max=value.max;
        if(min.length==1)
            return "new de.jtem.mfc.set.Interval("+min[0]+", "+max[0]+")";
        StringBuffer sb=new StringBuffer();
        sb.append("new de.jtem.mfc.set.Interval(new double[] {");
        for(int i=0, n=min.length; i<n; i++) sb.append(min[i]).append(", ");
        sb.append("}, new double[] {");
        for(int i=0, n=max.length; i<n; i++) sb.append(max[i]).append(", ");
        sb.append("})");
        return sb.toString();
    }
}
