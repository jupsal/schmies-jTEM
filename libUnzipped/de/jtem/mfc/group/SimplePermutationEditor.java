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

import de.jtem.numericalMethods.algebra.group.Permutation;
// import java.awt.Component;
/**
 * SimplePermutationEditor.java 
 *
 * <p>A textual editor for the <code>int[]</code> permutation property of a
 * {@link de.jtem.mfc.group.Permutation}.</p> 
 * 
 * Created: Thu Nov 21 14:42:03 2002
 *
 * @author <a href="mailto:mercat@sfb288.math.tu-berlin.de">Christian Mercat</a>
 * @version 1.0
 */
public final class SimplePermutationEditor extends PropertyEditorSupport
{

    /**
     * The Permutation value held by this editor. We use only the int[].
     *
     */
    private int[] value;
//     private PermutationEditorPanel editor;

    /**
     * To hold temporary things. Its size may increase if needed.
     *
     * @see #enlarge
     */
    private int[] basket = new int[20];

    public SimplePermutationEditor() {
	value=new int[0];
    }
    public SimplePermutationEditor(int[] v) {
	setValue(v);
    }
    public SimplePermutationEditor(de.jtem.mfc.group.Permutation v) {
	setValue(v.getPermutation());
    }
    /**
     * Returns a String that is encoding the value held by this editor. It is
     * very similar to toString().
     *
     * @return a <code>String</code> value
     * @see de.jtem.mfc.group.Permutation#toString
     */
    public String getAsText(){

	return Permutation.toString(value);
    }
    /**
     * The core of the editor, takes a String, cooks up a Permutation object
     * associated with it and transfers it to the value held by the editor.
     *
     * @param s a <code>String</code> value
     */
    public void setAsText(final  String s) throws IllegalArgumentException{
	final int[] theIntArray = stringTo(s);
	
	    if (Permutation.isPermutation
		(theIntArray)) {
		setValue(theIntArray);
	    } 	else
	  throw new IllegalArgumentException(
                      "Syntax: This is not a permutation of 0..n.");
	 
    }    
    public int[] stringTo(final  String s)
    {
	try{
	    return Permutation.stringToIntArray
		(s,basket);
	    }
	catch(ArrayIndexOutOfBoundsException ex){
	    enlarge();
	    return stringTo(s);
	}
    }

    /**
     * Doubles our basket's size.
     *
     */
    private void enlarge(){
	final int[] larger = new int[basket.length<<1];
	System.arraycopy(basket,0,larger,0,basket.length);
	basket = larger;
    }

    /**
     * Hands out a copy of the value held by this editor.
     *
     * @return an <code>Object</code> value
     */
    public Object getValue()
    {
	
	final int len = value.length;
	final int[] copy= new int[len];
	System.arraycopy(value,0,copy,0,len);

	return copy;
    }

    /**
     * Sets the value associated with this editor. Fires a PropertyChange if
     * it differs from the presently held value.
     *
     * @param o an <code>Object</code> value, cast to a Permutation. Will
     * throw an exception if the runtime type is incompatible.
     */
    public void setValue(final Object o)
    {

	final int[] newVal=  ((int[])o); 
	// Casts the object to an int[].
//     if(editor!=null)
//       editor.setValue(newVal);
	boolean equal = true;
	if(value.length!=newVal.length) {
	    equal = false;
	} else {
	    for (int  i = value.length;  --i>=0 ; ) {
		if (value[i]!=newVal[i]) {
	    equal = false;
	    break;
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

//   public Component getCustomEditor()
//   {
//     if(editor!=null)
//       return editor;
//     editor=new PermutationEditorPanel(this);
//     editor.setValue(value);
//     return editor;
//   }

//   public boolean supportsCustomEditor()
//   {
//     return true;
//   }

    //     public  void main(String[] args){ // To use the java debugger.
    // 	setAsText("");
    //     }

}// PermutationEditor

