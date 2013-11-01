/*
 *    Copyright (C) 2013 Geobeyond Srl
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.geotools.filter.function;

import org.geotools.filter.function.FluxoFilterFunction;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.geotools.feature.NameImpl;
import org.geotools.filter.FunctionFactory;
import org.opengis.feature.type.Name;
import org.opengis.filter.capability.FunctionName;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Function;
import org.opengis.filter.expression.Literal;

public class FluxoFunctionFactory implements FunctionFactory {    
    
    public List<FunctionName> getFunctionNames() {
        List<FunctionName> functionList = new ArrayList<FunctionName>();
        functionList.add(FluxoFilterFunction.NAME);        
        return Collections.unmodifiableList( functionList );
    }    
    public Function function(String name, List<Expression> args, Literal fallback) {
        return function(new NameImpl(name), args, fallback);
    }
    public Function function(Name name, List<Expression> args, Literal fallback) {
        if( FluxoFilterFunction.NAME.getFunctionName().equals(name)){
            return new FluxoFilterFunction( args, fallback );
        }
        return null; // we do not implement that function
    }
}
