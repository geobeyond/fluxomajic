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
import org.geotools.geometry.jts.JTS;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.opengis.filter.expression.Expression;
import org.opengis.filter.expression.Literal;

public class FluxoFilterFunction extends FunctionExpressionImpl implements
        GeometryTransformation {

    private static int quadrantSegments = 16;
    private static double mitreLimit = 10.0;
    public static FunctionName NAME = new FunctionNameImpl("fluxo5", Geometry.class,
            parameter("geometry", Geometry.class),
            parameter("offset", Double.class),
            parameter("width", Double.class));

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
        return 3;
    }

    @Override
    public Object evaluate(Object feature) {

        Geometry geom = getExpression(0).evaluate(feature, Geometry.class);

        Double offset = getExpression(1).evaluate(feature, Double.class);
        if (offset == null) {
            offset = 0d;
        }
        Double width = getExpression(2).evaluate(feature, Double.class);
        if (width == null) {
            width = 0d;
        }

        BufferParameters bufferparams = new BufferParameters();
        bufferparams.setSingleSided(false);
        bufferparams.setEndCapStyle(BufferParameters.CAP_ROUND);
        bufferparams.setJoinStyle(BufferParameters.JOIN_ROUND);
        bufferparams.setQuadrantSegments(quadrantSegments);
        bufferparams.setMitreLimit(mitreLimit);

        Geometry ret = bufferWithParams(offsetCurve(geom, offset, bufferparams, false), width, quadrantSegments, BufferParameters.CAP_ROUND, BufferParameters.JOIN_ROUND, mitreLimit);
        return ret;
    }

    private Geometry offsetCurve(Geometry geometry, double d, BufferParameters parameters, Boolean roughOffsetCurve) {
        GeometryFactory gf = geometry.getFactory();
        // If "geometry" is a surface, process its boundary
        if (geometry.getDimension() == 2) {
            geometry = geometry.getBoundary();
        }

        Collection offsetCurves = new ArrayList();
        if (roughOffsetCurve) {
            addRoughOffsetCurves(offsetCurves, geometry, parameters, d);
        } else {
            addCleanOffsetCurves(offsetCurves, geometry, parameters, d);
        }
        return gf.buildGeometry(offsetCurves);
    }

    private void addCleanOffsetCurves(Collection offsetCurves, Geometry sourceCurve, BufferParameters parameters, Double offsetDistance) {
        parameters.setSingleSided(true);
        parameters.setQuadrantSegments(quadrantSegments);
        Geometry sidedBuffer = new BufferOp(sourceCurve, parameters)
                .getResultGeometry(offsetDistance)
                .getBoundary();
        Collection offsetSegments = new ArrayList();
        // Segments located entirely under this distance are excluded
        double lowerBound = Math.abs(offsetDistance) * Math.sin(Math.PI / (4 * quadrantSegments));
        // Segments located entirely over this distance are included
        // note that the theoretical approximation made with quadrantSegments
        // is offset*cos(PI/(4*quadrantSegments) but offset*cos(PI/(2*quadrantSegments)
        // is used to make sure to include segments located on the boundary
        double upperBound = Math.abs(offsetDistance) * Math.cos(Math.PI / (2 * quadrantSegments));
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

    public static Geometry bufferWithParams(Geometry geometry, Double offset, Integer quadrantSegments, Integer capStyle, Integer joinStyle, Double mitreLimit) {
        double d = 0.0D;
        if (offset != null) {
            d = offset.doubleValue();
        }
        BufferParameters bufferparameters = new BufferParameters();

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
}
