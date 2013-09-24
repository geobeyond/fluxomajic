/*
 *    GeoTools - The Open Source Java GIS Toolkit
 *    http://geotools.org
 *
 *    (C) 2011, Open Source Geospatial Foundation (OSGeo)
 *    
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation;
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 */
package org.geotools.filter.function;

import com.vividsolutions.jts.algorithm.distance.DistanceToPoint;
import com.vividsolutions.jts.algorithm.distance.PointPairDistance;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import static org.geotools.filter.capability.FunctionNameImpl.*;
import org.geotools.filter.FunctionExpressionImpl;
import org.geotools.filter.capability.FunctionNameImpl;
import org.opengis.filter.capability.FunctionName;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geomgraph.Position;
import com.vividsolutions.jts.operation.buffer.BufferOp;
import com.vividsolutions.jts.operation.buffer.BufferParameters;
import com.vividsolutions.jts.operation.buffer.OffsetCurveBuilder;
import com.vividsolutions.jts.operation.linemerge.LineMerger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.geotools.geometry.DirectPosition2D;
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.referencing.CRS;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Literal;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

/**
 *
 * @author Francesco Bartoli (Geobeyond)
 */
public class FluxoFilterFunction extends FunctionExpressionImpl implements
        GeometryTransformation {

    private static int quadrantSegments = 16;
    private static double mitreLimit = 10.0;
    private static int OGC_DPI = 90;
    private static int RIGHT = 0;
    private static int LEFT = 1;
    public static FunctionName NAME = new FunctionNameImpl("fluxo", Geometry.class,
            parameter("geometry", Geometry.class),
            parameter("offset", Double.class),
            parameter("width", Double.class),
            parameter("driveMode", Integer.class),
            parameter("quadseg", Integer.class),
            parameter("endcap", Integer.class),
            parameter("join", Integer.class),
            parameter("outputCRS", CoordinateReferenceSystem.class),
            parameter("outputWidth", Integer.class),
            parameter("outputHeight", Integer.class),
            parameter("outputBBOX", ReferencedEnvelope.class));

    public FluxoFilterFunction() {
        super(NAME);
    }

    public FluxoFilterFunction(List<Expression> parameters, Literal fallback) {
        super(NAME);
        setParameters(parameters);
        setFallbackValue(fallback);
    }

    @Override
    public int getArgCount() {
        return 11;
    }

    @Override
    public Object evaluate(Object feature) {

        Geometry geom = getExpression(0).evaluate(feature, Geometry.class);

        Double offsetPx = getExpression(1).evaluate(feature, Double.class);
        if (offsetPx == null) {
            offsetPx = 0d;
        } else {
            offsetPx = Math.abs(offsetPx);
        }
        Double widthPx = getExpression(2).evaluate(feature, Double.class);
        if (widthPx == null) {
            widthPx = 0d;
        } else {
            widthPx = Math.abs(widthPx);
        }
        
        Integer dMode = getExpression(3).evaluate(feature, Integer.class);
        if (dMode == null) {
            dMode = RIGHT;
        } else if (dMode == 0) {
            dMode = RIGHT;
        } else if (dMode == 1) {
            dMode = LEFT;
        } else {
            dMode = RIGHT;
        }
        
        Integer quadseg = getExpression(4).evaluate(feature, Integer.class);
        if (quadseg == null) {
            quadseg = quadrantSegments;
        }

        Integer endcap = getExpression(5).evaluate(feature, Integer.class);
        if (endcap == null) {
            endcap = BufferParameters.CAP_ROUND;
        } else if (endcap == 2) {
            endcap = BufferParameters.CAP_FLAT;
        } else if (endcap == 3) {
            endcap = BufferParameters.CAP_SQUARE;
        } else {
            endcap = BufferParameters.CAP_ROUND;
        }
        
        Integer join = getExpression(6).evaluate(feature, Integer.class);
        if (join == null) {
            join = BufferParameters.JOIN_ROUND;
        } else if (join == 2) {
            join = BufferParameters.JOIN_MITRE;
        } else if (join == 3) {
            join = BufferParameters.JOIN_BEVEL;
        } else {
            join = BufferParameters.JOIN_ROUND;
        }
        
        CoordinateReferenceSystem outCRS = getExpression(7).evaluate(feature, CoordinateReferenceSystem.class);
        if (outCRS == null) {
            outCRS = DefaultGeographicCRS.WGS84;
        }
        
        //Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "outCRS:{0}", outCRS.toString());
        //Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "identificativo del sistema di coordinate di output SRS:{0}", outCRS.getCoordinateSystem().getIdentifiers().toString());
        
        Integer wmsWidth = getExpression(8).evaluate(feature, Integer.class);
        if (wmsWidth == null) {
            wmsWidth = 0;
        }
        
        Integer wmsHeight = getExpression(9).evaluate(feature, Integer.class);
        if (wmsHeight == null) {
            wmsHeight = 0;
        }
        
        Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore letto dalla geometria del db CRS:{0}", ((CoordinateReferenceSystem) geom.getUserData()).toString());
        
        Geometry finalGeom = null;
        try {
            finalGeom = transfGeom(geom, outCRS);
        } catch (NoSuchAuthorityCodeException ex) {
            Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
        } catch (FactoryException ex) {
            Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
        } catch (MismatchedDimensionException ex) {
            Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
        } catch (TransformException ex) {
            Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        //Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore SRID letto dalla finalGeom SRS:{0}", ((CoordinateReferenceSystem) finalGeom.getUserData()).toString());
        
        ReferencedEnvelope outBBox = getExpression(10).evaluate(feature, ReferencedEnvelope.class);
        if (outBBox == null) {
            outBBox = new ReferencedEnvelope(finalGeom.getEnvelopeInternal().getMinX(),finalGeom.getEnvelopeInternal().getMaxX(),finalGeom.getEnvelopeInternal().getMinY(),finalGeom.getEnvelopeInternal().getMaxY(),outCRS);
        }
        Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, outBBox.toString());
        
        double offsetMt;
        offsetMt = pixelToMeter(outBBox, wmsWidth, wmsHeight, offsetPx, outCRS);
                
        double widthMt;
        widthMt = pixelToMeter(outBBox, wmsWidth, wmsHeight, widthPx, outCRS);
        //Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, Double.toString(widthMt));
 
        
        double offsetCrs = distanceInCrs(offsetMt, finalGeom, outCRS);
        double widthCrs = distanceInCrs(widthMt, finalGeom, outCRS);

        BufferParameters bufferparams = new BufferParameters();
        bufferparams.setSingleSided(true);//metti false qui e inserisci la dichiarazione direttamente nella funzione privata
        bufferparams.setEndCapStyle(endcap);
        bufferparams.setJoinStyle(join);
        bufferparams.setQuadrantSegments(quadseg);
        bufferparams.setMitreLimit(mitreLimit);

        Geometry ret = null;
        if (doTravelLeft(dMode)) {
            ret = bufferWithParams(offsetCurve(finalGeom, -offsetCrs, bufferparams, false, bufferparams.getQuadrantSegments()), widthCrs, false, bufferparams.getQuadrantSegments(), endcap, join, mitreLimit);
        } else {
            ret = bufferWithParams(offsetCurve(finalGeom, offsetCrs, bufferparams, false, bufferparams.getQuadrantSegments()), widthCrs, false, bufferparams.getQuadrantSegments(), endcap, join, mitreLimit);
        }
        return ret;
    }

    private Geometry offsetCurve(Geometry geometry, double d, BufferParameters parameters, Boolean roughOffsetCurve, Integer qS) {
        GeometryFactory gf = geometry.getFactory();
        // If "geometry" is a surface, process its boundary
        if (geometry.getDimension() == 2) {
            geometry = geometry.getBoundary();
        }

        Collection offsetCurves = new ArrayList();
        if (roughOffsetCurve) {
            addRoughOffsetCurves(offsetCurves, geometry, parameters, d);
        } else {
            addCleanOffsetCurves(offsetCurves, geometry, parameters, d, qS);
        }
        return gf.buildGeometry(offsetCurves);
    }

    private void addCleanOffsetCurves(Collection offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance, Integer qS) {
        parameters.setSingleSided(true);
        parameters.setQuadrantSegments(qS);
        Geometry sidedBuffer = new BufferOp(sourceCurve, parameters)
                .getResultGeometry(offsetDistance)
                .getBoundary();
        Collection offsetSegments = new ArrayList();
        // Segments located entirely under this distance are excluded
        double lowerBound = Math.abs(offsetDistance) * Math.sin(Math.PI / (4 * qS));
        // Segments located entirely over this distance are included
        // note that the theoretical approximation made with quadrantSegments
        // is offset*cos(PI/(4*quadrantSegments) but offset*cos(PI/(2*quadrantSegments)
        // is used to make sure to include segments located on the boundary
        double upperBound = Math.abs(offsetDistance) * Math.cos(Math.PI / (2 * qS));
        for (int i = 0; i < sidedBuffer.getNumGeometries(); i++) {
            Coordinate[] cc = sidedBuffer.getGeometryN(i).getCoordinates();
            PointPairDistance ppd = new PointPairDistance();
            DistanceToPoint.computeDistance(sourceCurve, cc[0], ppd);
            double dj = ppd.getDistance();
            for (int j = 1; j < cc.length; j++) {
                double di = dj;
                ppd = new PointPairDistance();
                DistanceToPoint.computeDistance(sourceCurve, cc[j], ppd);
                dj = ppd.getDistance();
                // segment along or touching the source geometry : eclude it
                if (Math.max(di, dj) < lowerBound || di == 0 || dj == 0) {
                    continue;
                } // segment along the buffer boundary : include it
                else if (Math.min(di, dj) > upperBound) {
                    LineString segment = sourceCurve.getFactory().createLineString(
                            new Coordinate[]{cc[j - 1], cc[j]});
                    offsetSegments.add(segment);
                } // segment entirely located inside the buffer : exclude it
                else if (Math.min(di, dj) > lowerBound && Math.max(di, dj) < upperBound) {
                    continue;
                } // segment with a end at the offset distance and the other
                // located within the buffer : divide it
                else {
                    // One of the coordinates is closed to but not on the source
                    // curve and the other is more or less closed to offset distance
                    divide(offsetSegments, sourceCurve, cc[j - 1], cc[j], di, dj, lowerBound, upperBound);
                }
            }
        }
        offsetCurves.addAll(merge(offsetSegments));
    }

    // Recursive function to split segments located on the single-side buffer
    // boundary, but having a part of them inside the full buffer.
    private void divide(Collection offsetSegments, Geometry sourceCurve,
            Coordinate c1, Coordinate c2, double d1, double d2, double lb, double ub) {
        // I stop recursion for segment < 2*lb to exclude small segments
        // perpendicular but very close to the boundary
        if (c1.distance(c2) < 2 * lb) {
            return;
        }

        Coordinate c = new Coordinate((c1.x + c2.x) / 2.0, (c1.y + c2.y) / 2.0);
        PointPairDistance ppd = new PointPairDistance();
        DistanceToPoint.computeDistance(sourceCurve, c, ppd);
        double d = ppd.getDistance();
        if (Math.max(d1, d) < lb) {
        } else if (Math.min(d1, d) > lb && Math.max(d1, d) < ub) {
        } else if (Math.min(d1, d) > ub) {
            LineString segment = sourceCurve.getFactory().createLineString(
                    new Coordinate[]{c1, c});
            offsetSegments.add(segment);
        } else {
            divide(offsetSegments, sourceCurve, c1, c, d1, d, lb, ub);
        }
        if (Math.max(d, d2) < lb) {
        } else if (Math.min(d, d2) > lb && Math.max(d, d2) < ub) {
        } else if (Math.min(d, d2) > ub) {
            LineString segment = sourceCurve.getFactory().createLineString(
                    new Coordinate[]{c, c2});
            offsetSegments.add(segment);
        } else {
            divide(offsetSegments, sourceCurve, c, c2, d, d2, lb, ub);
        }
    }

    private void addRoughOffsetCurves(Collection offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance) {

        OffsetCurveBuilder builder = new OffsetCurveBuilder(
                sourceCurve.getFactory().getPrecisionModel(), parameters);

        for (int i = 0; i < sourceCurve.getNumGeometries(); i++) {
            if (sourceCurve.getGeometryN(i) instanceof LineString) {
                LineString lineString = (LineString) sourceCurve.getGeometryN(i);
                Coordinate[] cc = lineString.getCoordinates();
                if (lineString.isClosed()) {
                    offsetCurves.add(lineString.getFactory().createLineString(
                            builder.getRingCurve(cc,
                            offsetDistance > 0 ? Position.LEFT : Position.RIGHT,
                            Math.abs(offsetDistance))));
                } else {
                    offsetCurves.add(lineString.getFactory().createLineString(
                            builder.getOffsetCurve(cc, offsetDistance)));
                }
            }
        }
    }

    private Collection merge(Collection linestrings) {
        LineMerger merger = new LineMerger();
        merger.add(linestrings);
        return merger.getMergedLineStrings();
    }

    /**
     * Returns a buffered geometry with old shapes in the center of new ones.
     * If the buffer is issued at single side then a negative offset renders the
     * shape on the left while a positive offset on the right
     */
    public static Geometry bufferWithParams(Geometry geometry, Double offset, Boolean singleSided, Integer quadrantSegments, Integer capStyle, Integer joinStyle, Double mitreLimit) {
        double d = 0.0D;
        if (offset != null) {
            d = offset.doubleValue();
        }
        Boolean ss = false;
        if (singleSided != null) {
            ss = singleSided;
        }
   
        BufferParameters bufferparameters = new BufferParameters();
        
        //Inserimento custom per disegnare solo sul lato dell'offset della curva
        bufferparameters.setSingleSided(ss);
        
        if (quadrantSegments != null) {
            bufferparameters.setQuadrantSegments(quadrantSegments.intValue());
        }
        if (capStyle != null) {
            bufferparameters.setEndCapStyle(capStyle.intValue());
        }
        if (joinStyle != null) {
            bufferparameters.setJoinStyle(joinStyle.intValue());
        }
        if (mitreLimit != null) {
            bufferparameters.setMitreLimit(mitreLimit.doubleValue());
        }

        return BufferOp.bufferOp(geometry, d, bufferparameters);
    }

    /**
     * Returns an translated rendering envelope if the offsets are not using
     * feature attributes. If the offsets are feature dependent the user will
     * have to expand the rendering area via the renderer buffer parameter
     */
    @Override
    public ReferencedEnvelope invert(ReferencedEnvelope renderingEnvelope) {
        Envelope bufferedEnvelope = JTS.toGeometry((Envelope) renderingEnvelope).getEnvelopeInternal();
        return new ReferencedEnvelope(bufferedEnvelope, renderingEnvelope.getCoordinateReferenceSystem());
    }

    private double distanceInCrs(double inMeters, Geometry geom, CoordinateReferenceSystem crs) {
        try {
            double[] refXY = {
                geom.getCentroid().getCoordinates()[0].x,
                geom.getCentroid().getCoordinates()[0].y
            };
            //CoordinateReferenceSystem crs = CRS.decode("EPSG:4326");
            return distanceInCrs(inMeters, refXY, crs);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return 0d;
    }

    private double distanceInCrs(double inMeters,
            double[] refXY,
            CoordinateReferenceSystem crs) throws TransformException {
        double dist = 0;

        // calculate the distance in meters of 0.01 * refY in the ref CRS
        double[] sp = {refXY[0], refXY[1]};
        double[] dp = {refXY[0], refXY[1] * 1.01};

        GeodeticCalculator gc = new GeodeticCalculator(crs);

        gc.setStartingPosition(new DirectPosition2D(crs, sp[0], sp[1]));
        gc.setDestinationPosition(new DirectPosition2D(crs, dp[0], dp[1]));

        double refY01InMeters = gc.getOrthodromicDistance();

        // now, calculate the CRS distance as a proportional of 0.01 * refY
        dist = inMeters * (refXY[1] * 0.01) / refY01InMeters;

        return dist;
    }
    private double pixelSize(ReferencedEnvelope outputEnv, int outputWidth, int outputHeight)
    {
        // error-proofing
        if (outputEnv.getWidth() <= 0) return 0;
        // assume view is isotropic
        return outputWidth / outputEnv.getWidth();
    }
    private double pixelToMeter(ReferencedEnvelope outputEnv, int outputWidth, int outputHeight, double pixel_distance, CoordinateReferenceSystem crs)
    {
        double pixel_distance_m;
        double pixel_diag_distance;
        double pixel_diag_distance_m; 
        
        pixel_diag_distance = Math.sqrt((outputWidth * outputWidth)
                                      + (outputHeight * outputHeight));
        pixel_diag_distance_m = getGeodeticSegmentLength(outputEnv.getMinX(), outputEnv.getMinY(), outputEnv.getMaxX(), outputEnv.getMaxY(), crs);
        pixel_distance_m = pixel_diag_distance_m * pixel_distance / pixel_diag_distance;
        return pixel_distance_m;
        
    }
    private static double getGeodeticSegmentLength(double minx, double miny, double maxx, double maxy, CoordinateReferenceSystem crs) {
        final GeodeticCalculator calculator = new GeodeticCalculator(DefaultGeographicCRS.WGS84);
        //final GeodeticCalculator calculator = new GeodeticCalculator(crs);
        double rminx = rollLongitude(minx);
        double rminy = rollLatitude(miny);
        double rmaxx = rollLongitude(maxx);
        double rmaxy = rollLatitude(maxy);
        calculator.setStartingGeographicPoint(rminx, rminy);
        calculator.setDestinationGeographicPoint(rmaxx, rmaxy);
        return calculator.getOrthodromicDistance();
    }
    protected static double rollLongitude(final double x) {
        double rolled = x - (((int) (x + Math.signum(x) * 180)) / 360) * 360.0;
        return rolled;
    }
    
    protected static double rollLatitude(final double x) {
        double rolled = x - (((int) (x + Math.signum(x) * 90)) / 180) * 180.0;
        return rolled;
    }
    
    /**
     * Returns a geometry based on the transformation from a source geometry CRS to a defined target CRS.
     * 
     */
    private Geometry transfGeom(Geometry g, CoordinateReferenceSystem outputCRS) throws NoSuchAuthorityCodeException, FactoryException, MismatchedDimensionException, TransformException {
        CoordinateReferenceSystem srcCRS = null;
        
        if (srcCRS == null) {
            try {
                srcCRS = (CoordinateReferenceSystem) g.getUserData();
                Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore srcCRS letto dalla geometria sorgente CRS:{0}", srcCRS.toString());
            } catch (Exception e) {
                // may not have a CRS attached
            }
        }
//        if (srcCRS == null && g.getSRID() > 0) {
//            try {
//                srcCRS = CRS.decode("EPSG:" + g.getSRID());
//                Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore srcCRS letto dalla geometria sorgente CRS:{0}", srcCRS.toString());
//            } catch (Exception e) {
//                // may not have a CRS attached
//            }
//        }
        //try to force getting source/dest CRS
        srcCRS = CRS.decode("EPSG:4326");
        Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore forzato del CRS sorgente:{0}", srcCRS.toString());
        outputCRS = CRS.decode("EPSG:900913");
        Logger.getLogger(FluxoFilterFunction.class.getName()).log(Level.INFO, "valore forzato del CRS di output:{0}", outputCRS.toString());
        MathTransform transform;
        transform = CRS.findMathTransform(srcCRS, outputCRS, false);
        Geometry trgGeom = JTS.transform(g, transform);
        
        return trgGeom;
    }
    private Boolean doTravelLeft(Integer i) {
        if (i == 1) return true; 
        else return false;
    }
}
