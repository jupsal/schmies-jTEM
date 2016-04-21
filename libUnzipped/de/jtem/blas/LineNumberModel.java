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

/**
 * An instance of this class is a TableModel object used for one columned
 * tables which cell content is restricted to integer values.
 */
class LineNumberModel extends javax.swing.table.AbstractTableModel
{
  int numRows;
  /**
   * Returns the type of class which instances can be content of a cell of the
   * table, in this case only integer values.
   * 
   * @return <code>Integer.class</code>
   */
  public Class getColumnClass(int columnIndex)
  {
    return Integer.class;
  }
  /**
   * Returns the number of colums the table consists of, in this case only one
   * column.
   * 
   * @return <code>1</code>
   */
  public int getColumnCount()
  {
    return 1;
  }
  /**
   * @return <code>Line</code>
   */
  public String getColumnName(int columnIndex)
  {
    return "Line";
  }
  /**
   * Returns the number of rows the table consists of.
   */
  public int getRowCount()
  {
    return numRows;
  }
  /**
   * Sets the number of rows the table consists of to the specified int value.
   */
  public void setRowCount(int rows)
  {
    if(rows==numRows) return;
    int old=numRows;
    numRows=rows;
    if(rows>old)
      fireTableRowsInserted(old, rows);
    else
      fireTableRowsDeleted(rows, old);
  }
  /**
   * Returns the value for the cell at columnIndex and rowIndex.
   */
  public Object getValueAt(int rowIndex, int columnIndex)
  {
    return new Integer(rowIndex);
  }
}
