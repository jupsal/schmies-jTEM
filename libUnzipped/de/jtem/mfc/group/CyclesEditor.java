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



import java.beans.PropertyEditorSupport;


/**
 * CyclesEditor.java 
 *
 * <p>A textual editor for
 * {@link de.jtem.mfc.group.Permutation} as a cycle.</p> 
 * 
 * Created: Thu Nov 21 14:42:03 2002
 *
 * @author <a href="mailto:mercat@sfb288.math.tu-berlin.de">Christian Mercat</a> * @version 1.0
 */
public final class CyclesEditor extends PropertyEditorSupport
{

    /**
     * The cycles value held by this editor. 
     *
     */
    private int[][] value;
//     private CyclesEditorPanel editor;


    public CyclesEditor() {
	value=new int[0][];
    }
    public CyclesEditor(int[][] v) {
	setValue(v);
    }
    public CyclesEditor(int[] p) {
	setValue(de.jtem.numericalMethods.algebra.group.Permutation.cycles(p));
    }
    public CyclesEditor(Permutation v) {
	setValue(v.getCycles());
    }
    /**
     * Returns a String that is encoding the value held by this editor. It is
     * very similar to toString().
     *
     * @return a <code>String</code> value
     * @see Permutation#toString
     */
    public String getAsText(){
	return de.jtem.numericalMethods.algebra.group.Permutation.cyclesToString
	    (value);
    }
    /**
     * The core of the editor, takes a String representing cycles, cooks up a
     * Permutation object associated with it and transfers it to the value
     * held by the editor.
     *
     * @param s a <code>String</code> value
     */
    public void setAsText(final  String s) throws IllegalArgumentException{
		setValue
		    (de.jtem.numericalMethods.algebra.group.Permutation.stringToCycles
		     (s));
	 
    }    

    /**
     * Hands out a copy of the value held by this editor.
     *
     * @return an <code>Object</code> value
     */
    public Object getValue()
    {
	final int[][] copy= new int[value.length][];
	for (int  i = value.length;  --i>=0 ; ) {
	    final int len = value[i].length;
	    copy[i]= new int[len];
	    System.arraycopy(value[i],0,copy[i],0,len);
	}
	   
	return copy;
    }

    /**
     * Sets the value associated with this editor. Fires a PropertyChange if
     * it differs from the presently held value.
     *
     * @param o an <code>Object</code> value, cast to an int[][]. Will
     * throw an exception if the runtime type is incompatible.
     */
    public void setValue(final Object o)
    {
	final int[][] newVal= (int[][])o; 
	// Casts the object to an int[][].
//     if(editor!=null)
//       editor.setValue(newVal);
{

	boolean equal = true;
	if(value.length!=newVal.length) {
	    equal = false;
	} else {
	    iLoop:
	    for (int  i = value.length;  --i>=0 ; ) {
		if (value[i].length!=newVal[i].length) {
		    equal = false;
		    break;
		} else {
		    for (int  j = value[i].length;  --j>=0 ; ) {
			if (value[i][j]!=newVal[i][j]) {
			    equal = false;
			    break iLoop;
			}
		    }
		}
	    }
	}

	if(!equal)
	    {
// 		super.setValue(o);
		
		value=newVal; // Sets the value.
		firePropertyChange();   // Lets the listeners know.
	    }
    	
    }

    }
    public void setValue(final int[] newVal){
	setValue(de.jtem.numericalMethods.algebra.group.Permutation.cycles(newVal));
    }
    

//   public Component getCustomEditor()
//   {
//     if(editor!=null)
//       return editor;
//     editor=new CyclesEditorPanel(this);
//     editor.setValue(value);
//     return editor;
//   }

  public boolean supportsCustomEditor()
  {
      return false;
//     return true;
  }

    //     public  void main(String[] args){ // To use the java debugger.
    // 	setAsText("");
    //     }

}// CyclesEditor

