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
import java.beans.BeanDescriptor;
import java.beans.IndexedPropertyDescriptor;
import java.beans.IntrospectionException;
import java.beans.PropertyDescriptor;
import java.beans.SimpleBeanInfo;
/**
 * BeanInfo for ComplexMatrix.
 */
public class ComplexMatrixBeanInfo extends SimpleBeanInfo
{
  private static final Class beanClass=ComplexMatrix.class;
  private final PropertyDescriptor[] pd
    =new PropertyDescriptor[7];
  private final BeanDescriptor bd=new BeanDescriptor(beanClass,ComplexMatrixCustomizer.class);
  public ComplexMatrixBeanInfo() throws IntrospectionException
  {
    try
    {
      System.out.println("ComplexMatrixBeanInfo.<init>");
      PropertyDescriptor p;
      java.lang.reflect.Method get=
        beanClass.getMethod("getNumCols", new Class[0]);
      p=new PropertyDescriptor("numCols",get,null);
      p.setDisplayName("number of columns");
      pd[0]=p;
      get=
        beanClass.getMethod("getNumRows", new Class[0]);
      p=new PropertyDescriptor("numRows",get,null);
      p.setDisplayName("number of rows");
      pd[1]=p;
      get=
        beanClass.getMethod("getNumEntries", new Class[0]);
      p=new PropertyDescriptor("numEntries",get,null);
      p.setDisplayName("number of entries");
      pd[2]=p;
      p=new IndexedPropertyDescriptor("row",null,null,
        beanClass.getMethod("getRow", new Class[] { int.class }),
        beanClass.getMethod("setRow", new Class[] { int.class, ComplexVector.class }));
      p.setDisplayName("row as vec");
      pd[3]=p;
      p=new IndexedPropertyDescriptor("col",null,null,
        beanClass.getMethod("getCol", new Class[] { int.class }),
        beanClass.getMethod("setCol", new Class[] { int.class, ComplexVector.class }));
      p.setDisplayName("column as vec");
      pd[4]=p;
      p=new PropertyDescriptor("square",beanClass,"isSquared",null);
      p.setDisplayName("square matrix ?");
      pd[5]=p;
      p=new PropertyDescriptor("symmetric",beanClass,"isSymmetric",null);
      p.setDisplayName("symmetric matrix ?");
      pd[6]=p;
      bd.setDisplayName("Field.Complex AbstractMatrix");
      bd.setShortDescription("A matrix of complex numbers.");
    }
    catch(NoSuchMethodException ex)
    {
      System.out.println(ex);
      throw new IntrospectionException("method missing");
    }
    catch(IntrospectionException ex)
    {
      System.out.println(ex);
      throw ex;
    }

  }

  public PropertyDescriptor[] getPropertyDescriptors()
  {
    return pd;
  }

  public int getDefaultPropertyIndex()
  {
    return -1;
  }

  public BeanDescriptor getBeanDescriptor()
  {
    return bd;
  }
}

