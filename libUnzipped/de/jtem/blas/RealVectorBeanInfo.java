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
 * BeanInfo for RealVector.
 */
public class RealVectorBeanInfo extends SimpleBeanInfo
{
  private static final Class beanClass=RealVector.class;
  private final PropertyDescriptor[] pd
    =new PropertyDescriptor[2];
  private final BeanDescriptor bd=new BeanDescriptor(beanClass,RealVectorCustomizer.class);
  public RealVectorBeanInfo() throws IntrospectionException
  {
    System.out.println("RealVectorBeanInfo.<init>");
    PropertyDescriptor p=new PropertyDescriptor("size",beanClass,"size",null);
    pd[0]=p;
    try
    {
      java.lang.reflect.Method set=
        beanClass.getMethod("set",new Class[] { int.class, double.class });
      java.lang.reflect.Method get=
        beanClass.getMethod("get",new Class[] { int.class });
      java.lang.reflect.Method getA=
        beanClass.getMethod("toArray",new Class[0] );
      p=new IndexedPropertyDescriptor("element",getA,null,get,set);
    System.out.println("RealVectorBeanInfo.<init>: all methods found");
    }
    catch(NoSuchMethodException ex)
    {
    System.out.println("RealVectorBeanInfo.<init>: not all methods found");
      throw 
        new IntrospectionException("missing method defined in this beaninfo");
    }
    p.setDisplayName("vector elements");
    pd[1]=p;
    bd.setDisplayName("Double AbstractVector");
    bd.setShortDescription("A AbstractVector of Double numbers.");
    System.out.println("RealVectorBeanInfo.<init>: ok");
  }

  public PropertyDescriptor[] getPropertyDescriptors()
  {
    return pd;
  }

  public int getDefaultPropertyIndex()
  {
    return 1;
  }

  public BeanDescriptor getBeanDescriptor()
  {
    return bd;
  }
}

