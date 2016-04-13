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

package de.jtem.mfc.field;

import java.awt.Component;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyEditorSupport;

/**
 * Graphic editor for Quaternion.
 */
public class QuaternionEditor extends PropertyEditorSupport
    implements PropertyChangeListener
{
  QuaternionEditorPanel editor;
  final Quaternion value;
  private Quaternion editorPanelValue=new Quaternion();
  public QuaternionEditor()
  {
    value=new Quaternion();
  }
  public void setValue(Object obj)
  {
    Quaternion newValue=(Quaternion)obj;
    if(!newValue.equals(value))
    {
      value.assign(newValue);
      if(editor!=null)
        editor.setQuaternion(value);// implies sending an event
      else
        firePropertyChange();
    }
  }
  public Object getValue()
  {
    return new Quaternion(value);
  }
  public Component getCustomEditor()
  {
    if(editor!=null)
      return editor;
    editor=new QuaternionEditorPanel();
    editor.setQuaternion(value);
    editor.addPropertyChangeListener(
        QuaternionEditorPanel.QUATERNION_PROPERTY, this);
    return editor;
  }
  public boolean supportsCustomEditor()
  {
    return true;
  }
  public String getJavaInitializationString()
  {
    return "new de.jtem.mfc.field.Quaternion("+value.re+", "+value.x+", "
      +value.y+", "+value.z+")";
  }

  /**
   * This QuaternionEditor listens to used QuaternionEditorPanel and
   * this method is called, if a the "quaternion" property of
   * QuaternionEditorPanel is changed.
   */
  public void propertyChange(PropertyChangeEvent ev)
  {
    editor.getQuaternion(editorPanelValue);
    value.assign(editorPanelValue);
    firePropertyChange();
  }
}
