// JTEM - Java Tools for Experimental Mathematics
// Copyright (C) 2001 JEM-Group
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

package de.jtem.riemann.schottky;

import java.beans.BeanDescriptor;
import java.beans.IntrospectionException;
import java.beans.SimpleBeanInfo;

public class SchottkyBeanInfo extends SimpleBeanInfo {

    private static final Class beanClass = Schottky.class;
    
    private final BeanDescriptor bd = new BeanDescriptor
	( beanClass,  SchottkyCustomizer.class );


    public SchottkyBeanInfo() throws IntrospectionException {
    }


    public int getDefaultPropertyIndex () {
	return -1;
    } 

    public BeanDescriptor getBeanDescriptor () {
	return bd;
    } 
}



