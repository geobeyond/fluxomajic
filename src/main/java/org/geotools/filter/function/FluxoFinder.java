/*  
 *  
 *  Copyright (C) 2013  Geobeyond
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * 
 */

package org.geotools.filter.function;

import java.util.Iterator;
import java.util.Set;

import org.geotools.factory.FactoryCreator;
import org.geotools.factory.FactoryFinder;
import org.geotools.factory.FactoryRegistry;
import org.geotools.factory.Hints;
import org.geotools.filter.FunctionFactory;
import org.geotools.resources.LazySet;
import org.opengis.filter.FilterFactory;

public class FluxoFinder extends FactoryFinder {

    private static FactoryCreator registry;
    
    private static FactoryRegistry getServiceRegistry() {
        assert Thread.holdsLock(FluxoFinder.class);
        if (registry == null) {
            Class<?> categories[] = new Class<?>[] { FunctionFactory.class };
            registry = new FactoryCreator( categories);
        }
        return registry;
    }
    
    /**
     * Returns a set of all available implementations for the {@link FilterFactory} interface.
     *
     * @param  hints An optional map of hints, or {@code null} if none.
     * @return Set of available filter factory implementations.
     */
    public static synchronized Set<FunctionFactory> getFilterFactories(Hints hints) {
        hints = mergeSystemHints(hints);
        Iterator<FunctionFactory> serviceProviders = getServiceRegistry().getServiceProviders(
                FunctionFactory.class, null, hints);
        return new LazySet<FunctionFactory>(serviceProviders);
    }
        
    /** Allow the classpath to be rescanned */
    public static synchronized void scanForPlugins() {
        if (registry != null) {
            registry.scanForPlugins();
        }
    }
}
