import org.apache.commons.math4.legacy.fitting.WeightedObservedPoints;
import org.apache.commons.math4.legacy.fitting.PolynomialCurveFitter;

public class CalculateValues {
    // Solving constants
    private static double x; // meters (m)
    private static double y; // meters (m)
    private static double z; // meters (m)
    private static double theta; // radians (rad)

    private static double Vx; // meters per second (m/s)
    private static double Vy; // meters per second (m/s)
    private static double Vz; // meters per second (m/s)

    private static double Ax; // meters per second squared (m/s^2)
    private static double Ay; // meters per second squared (m/s^2)
    private static double Az; // meters per second squared (m/s^2)

    private static double time; // seconds (s)

    // Given constants
    private static double Vi; // meters per second (m/s)

    // Other constants
    private static final double deltaH = 1; // meters (m)
    private static final double lambda = 45 * (Math.PI) / 180; // radians (rad)

    private static final double m = 0.215; // kilograms (kg)
    private static final double g = -9.8; // meters per second squared (m/s^2)

    private static final double c = 0.39; // dimensionless (drag coefficient)
    private static final double rho = 1.225; // kilograms per cubic meter (kg/m^3)
    private static final double area = 0.0177;// square meters (m^2)

    private static final double drag = 0.5 * c * rho * area; // drag constant

    private static final double dt = 0.01; // seconds (s)

    public static double calculateX(double Vi) {
        boolean hasReachedPeak = false;

        x = 0;
        y = 0;
        z = 0;
        theta = 0;

        Vx = Vi * Math.cos(lambda);
        Vy = Vi * Math.sin(lambda);
        Vz = 0;

        Ax = -(drag * Vx * Vx) / m;
        Ay = g - (drag * Vy * Vy * Math.signum(Vy)) / m;
        Az = 0;

        time = 0;

        while (!hasReachedPeak) {
            time += dt;

            Vx += Ax * dt;
            Vy += Ay * dt;

            Ax = -(drag * Vx * Vx) / m;
            Ay = g - (drag * Vy * Vy * Math.signum(Vy)) / m;

            x += Vx * dt;
            y += Vy * dt;

            hasReachedPeak = Vy <= 0;
        }

        while (y > deltaH) {
            time += dt;

            Vx += Ax * dt;
            Vy += Ay * dt;

            Ax = -(drag * Vx * Vx) / m;
            Ay = g - (drag * Vy * Vy * Math.signum(Vy)) / m;

            x += Vx * dt;
            y += Vy * dt;
        }

        return x;
    }

    public static void main(String[] args) {
        double ViMin = 5.0; // minimum initial velocity (m/s)
        double ViMax = 5.5; // maximum initial velocity (m/s)
        double step = 1; // step size in m/s

        int numPoints = (int) ((ViMax - ViMin) / step) + 1;
        double[] velocities = new double[numPoints];
        double[] deltaX = new double[numPoints];

        for (int i = 0; i < numPoints; i++) {
            double Vi = ViMin + i * step;
            velocities[i] = Vi;
            deltaX[i] = calculateX(Vi);
        }

        // Map Δx -> Vi for polynomial fitting
        WeightedObservedPoints points = new WeightedObservedPoints();
        for (int i = 0; i < numPoints; i++) {
            points.add(deltaX[i], velocities[i]); // Δx as x, Vi as y
        }

        // Fit cubic polynomial
        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(4);
        double[] coeff = fitter.fit(points.toList());

        // System.out.println("Quartic fit (x -> Vi):");
        // System.out.printf(
        //         "Vi ≈ %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3 + %.6f*x^4%n",
        //         coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);

        // double xTarget = 8.0;

        // double ViEstimate = coeff[0]
        //         + coeff[1] * xTarget
        //         + coeff[2] * xTarget * xTarget
        //         + coeff[3] * xTarget * xTarget * xTarget
        //         + coeff[4] * xTarget * xTarget * xTarget * xTarget;

        // System.out.println("Estimated Vi for x = " + xTarget + " m: " + ViEstimate + " m/s");

        // System.out.println("x for estimated vi: " + calculateX(ViEstimate));

        System.out.println(calculateX(0.1));
    }
}
