/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
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

package de.jtem.blas;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexEditor;

class ComplexVectorTableModel extends AbstractVectorTableModel
{
  ComplexEditor editor;

  public ComplexVectorTableModel()
  {
    vector=new ComplexVector(1);
    editor=new ComplexEditor();
  }

  public synchronized Object getValueAt(int row, int column)
  {
    Complex value = ((ComplexVector)vector).get(row);
    editor.setValue( value );
    //return editor.getAsText();
    return value;
  }
  public Class getColumnClass(int columnIndex)
  {
    return Complex.class;
  }
  public synchronized
         void setValueAt(Object aValue, int rowIndex, int columnIndex)
  {
    //editor.setAsText((String)aValue);
    Complex value = (Complex)aValue;
    editor.setValue( value );
    //((ComplexVector)vector).set(rowIndex,(Field.Complex)editor.getValue());
    ((ComplexVector)vector).set(rowIndex, value );
    fireTableCellUpdated(rowIndex, 0);
  }
}



