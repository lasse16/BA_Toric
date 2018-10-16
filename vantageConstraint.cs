using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class vantageConstraint
{
    private Vector3 _AB;
    private Vector3 A;
    private Vector3 B;
    private Vector3 _targetPosition;

	public vantageConstraint(Vector3 A , Vector3 B)
	{
        _AB = B-A;
	}

    public Dictionary<float, Interval> getPositionFromVantageOneTarget(Vector3 targetPosition, Vector3 v, float deviationAngle, float samplingRate = 0.01f)
    {
        _targetPosition = targetPosition;

        if (v.Equals(Vector3.zero)) throw new Exception("No Vantage constraint set");
        if (deviationAngle > 90) throw new Exception("choose a deviattion angle smaller 90°");


        Vector3 upOnPlane, Vflat, PHIflat, vectorOnPlane;

        float phiVector = 0;



        Plane targetlevel = new Plane(Vector3.up, _targetPosition);

        vectorOnPlane = targetlevel.ClosestPointOnPlane(_targetPosition + v) - _targetPosition;

        upOnPlane = Vector3.ProjectOnPlane(Vector3.up, _AB);

        Vector3 phiZero = Vector3.Cross(upOnPlane, _AB);



        Vflat = Vector3.ProjectOnPlane(v.normalized, _AB);
        PHIflat = Vector3.ProjectOnPlane(phiZero.normalized, _AB);
        phiVector = Vector3.SignedAngle(Vflat, PHIflat, -_AB) * Mathf.Deg2Rad;


        v = vectorOnPlane;


        Vector3 prefferedVantageAngle = v.normalized;


        float signedLambda = Vector3.SignedAngle(_AB, prefferedVantageAngle, upOnPlane) * Mathf.Deg2Rad;

        float lambda = Mathf.Abs(signedLambda);
        deviationAngle *= Mathf.Deg2Rad;


        Vector2 x, y, y2;
        float a, b, c;
        Interval phiIntervall;
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();
        List<Vector2> possibleIntersections = new List<Vector2>();







        Boolean greaterThanHalfPi = lambda > Mathf.PI / 2;
        if (greaterThanHalfPi) lambda = Mathf.PI - lambda;



        Cone vantageCone = new Cone(prefferedVantageAngle, deviationAngle, _targetPosition);
        float r = Mathf.Tan(lambda);



        switch (checkConicSection(lambda, vantageCone))
        {

            case -3: //circle
                Debug.Log("Circle");
                phiIntervall = new Interval(-Mathf.PI, Mathf.PI, samplingRate);

                foreach (float phi in phiIntervall.getEveryValue())
                {
                    float x1 = Mathf.Tan(deviationAngle) * Mathf.Cos(phi);
                    float y1 = Mathf.Tan(deviationAngle) * Mathf.Sin(phi);
                    res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x1, y1), new Vector2(-x1, -y1), greaterThanHalfPi));
                }

                break;
            case -2: //ellipse
                Debug.Log("Ellipse");

                //Ellipses too small, are treated like circles 
                float lambdaTest = lambda;
                if (lambda > Mathf.PI / 2) lambdaTest = Mathf.PI - lambda;
                if (deviationAngle > lambdaTest) goto case -3;


                float sinDeviation = Mathf.Sin(deviationAngle);
                float coslambda = Mathf.Cos(lambda);
                float majorDistance = (sinDeviation * Mathf.Cos(deviationAngle)) / (Mathf.Abs(Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2)));
                float minorDistance = sinDeviation / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2))));
                float midPointDistance = (Mathf.Sin(lambda) * coslambda) / Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2);








                a = Mathf.Pow(minorDistance, 2) - Mathf.Pow(majorDistance, 2);
                b = -2 * Mathf.Pow(minorDistance, 2) * Mathf.Abs(midPointDistance);
                c = Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) + Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) - Mathf.Pow(minorDistance, 2));




                x = solveQuadraticEquation(a, b, c);


                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;
                y2.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.y, 2)));
                y2.y = -y2.x;
                y = sortVector2(y);
                y2 = sortVector2(y2);


                possibleIntersections = combineAllXandYValues(x, y, y2);



                phiIntervall = appropiatePhiIntervalBounds(new Vector2(Mathf.Abs(possibleIntersections.First().x), Mathf.Abs(possibleIntersections.First().y)), samplingRate);


                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * (minorDistance / majorDistance) * Mathf.Sqrt(Mathf.Pow(majorDistance, 2) - Mathf.Pow(midPointDistance, 2));
                        y.y = Mathf.Sin(phi) * (minorDistance / majorDistance) * -Mathf.Sqrt(Mathf.Pow(majorDistance, 2) - Mathf.Pow(midPointDistance, 2));
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }
                    else
                    {
                        a = Mathf.Pow(minorDistance, 2) + Mathf.Pow(majorDistance, 2) * Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;
                        c = Mathf.Pow(minorDistance, 2) * (Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2));

                        x = solveQuadraticEquation(a, b, c);






                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.x)) x.x *= -1;
                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.y)) x.y *= -1;
                        y.x = Mathf.Tan(phi) * x.x;
                        y.y = Mathf.Tan(phi) * x.y;


                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }              
                }

                break;

            case -1: //conic section = par_ABola
                Debug.Log("Par_ABola");

                float cotlambda = 1 / Mathf.Tan(lambda);
                float h = (Mathf.Tan(lambda) - cotlambda) / 2;
                float f = cotlambda / 2;
                h *= Mathf.Sign(signedLambda);

                Debug.Log("Height : " + h);

                x = solveQuadraticEquation(-1, -4 * f, 4 * f * Mathf.Abs(h) + Mathf.Pow(r, 2));
                x.x *= Mathf.Sign(signedLambda);
                x.y *= Mathf.Sign(signedLambda);

                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;
                y2.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.y, 2)));
                y2.y = -y2.x;
                y = sortVector2(y);
                y2 = sortVector2(y2);
                possibleIntersections = combineAllXandYValues(x, y, y2);





                float intersectionX = possibleIntersections.First().x;
                float intersectionY = possibleIntersections.First().y;


                phiIntervall = appropiatePhiIntervalBounds(new Vector2(intersectionX, intersectionY), samplingRate);

                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * 2 * Mathf.Sqrt(-f * h);
                        y.y = -y.x;
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }
                    else
                    {
                        a = Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -4 * f;
                        c = 4 * f * h;
                        x = solveQuadraticEquation(a, b, c);

                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.x)) x.x *= -1;
                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.y)) x.y *= -1;

                        y.x = Mathf.Tan(phi) * x.x;
                        y.y = Mathf.Tan(phi) * x.y;
                        //DEBUG
                        Debug.Log(y);

                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }

                }

                break;

            case 0: //conic section = hyperbola      
                Debug.Log("Hyperbola");
                sinDeviation = Mathf.Sin(deviationAngle);
                coslambda = Mathf.Cos(lambda);
                majorDistance = (sinDeviation * Mathf.Cos(deviationAngle)) / (Mathf.Abs(Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2)));
                minorDistance = sinDeviation / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2))));
                midPointDistance = (Mathf.Sin(lambda) * coslambda) / Mathf.Pow(coslambda, 2) - Mathf.Pow(sinDeviation, 2);




                a = Mathf.Pow(minorDistance, 2) + Mathf.Pow(majorDistance, 2);
                b = -2 * Mathf.Pow(minorDistance, 2) * Mathf.Abs(midPointDistance);
                c = Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) + Mathf.Pow(minorDistance, 2));

                x = solveQuadraticEquation(a, b, c);


                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));

                y.y = -y.x;
                y2.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.y, 2)));
                y2.y = -y2.x;
                y = sortVector2(y);
                y2 = sortVector2(y2);


                possibleIntersections = combineAllXandYValues(x, y, y2);



                if (possibleIntersections.Count == 0) throw new Exception("No solution possible");

                intersectionX = possibleIntersections.First().x;
                intersectionY = possibleIntersections.First().y;


                phiIntervall = appropiatePhiIntervalBounds(new Vector2(Mathf.Abs(intersectionX), Mathf.Abs(intersectionY)), samplingRate);


                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * (minorDistance / majorDistance) * Mathf.Sqrt(Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2)); ;
                        y.y = -y.x;
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }
                    else
                    {
                        a = Mathf.Pow(minorDistance, 2) - Mathf.Pow(majorDistance, 2) * Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;


                        x = solveQuadraticEquation(a, b, Mathf.Pow(minorDistance, 2) * (Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2)));



                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.x)) x.x *= -1;
                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.y)) x.y *= -1;

                        y.x = Mathf.Tan(phi) * x.x;
                        y.y = Mathf.Tan(phi) * x.y;





                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }
                }

                break;

            case 1: //conic section = line
                Debug.Log("Line");
                res.Add(0, new Interval(Mathf.PI / 2, Mathf.PI / 2));
                break;

            case 2: //conic section = two lines
                Debug.Log("Two lines");
                phiIntervall = new Interval(-deviationAngle, deviationAngle);
                foreach (float phi in phiIntervall.getEveryValue())
                {
                    float plane = Mathf.PI * 1.5f;
                    if (signedLambda < 0) plane = Mathf.PI / 2;
                    res.Add(phi, new Interval(plane - deviationAngle, plane + deviationAngle));
                }
                break;
            default:
                throw new Exception();

        }
        foreach (float key in res.Keys.ToArray())
        {
            Interval beta;
            res.TryGetValue(key, out beta);
            res.Remove(key);
            res.Add(key + phiVector, beta);
        }


        return res;
    }


    private List<Vector2> combineAllXandYValues(Vector2 x, Vector2 y, Vector2 y2)
    {

        List<Vector2> res = new List<Vector2>();

        Debug.Log(x);
        Debug.Log(y);
        Debug.Log(y2);


        if (!float.IsNaN(y.x))
        {
            Vector2 x1y11 = new Vector2(x.x, y.x);
            Vector2 x1y12 = new Vector2(x.x, y.y);
            res.Add(x1y11); res.Add(x1y12);
        }

        if (!float.IsNaN(y2.x))
        {
            Vector2 x2y21 = new Vector2(x.y, y2.x);
            Vector2 x2y22 = new Vector2(x.y, y2.y);
            res.Add(x2y21); res.Add(x2y22);
        }

        return res;

    }

    private void Visualize(Vector3 __targetPosition, float r, float lambda, Cone vantageCone, List<Vector2> possibleIntersections, float majorDistance, float minorDistance)
    {

        float coneHeight = Mathf.Sqrt(1 + r * r);


        vantageCone.Draw(Color.red, coneHeight);
        Vector3 intersectionVantagePlane = __targetPosition + vantageCone._midline * coneHeight;
        Debug.DrawLine(__targetPosition, intersectionVantagePlane, Color.black, Mathf.Infinity);
        Debug.DrawLine(__targetPosition, __targetPosition + -_AB, Color.black, Mathf.Infinity);

        Vector3 intersection_ABPlane = __targetPosition + _AB.normalized * checkBetaForPlane(lambda);
        DrawPlane(intersection_ABPlane, -_AB.normalized);

        foreach (Vector2 item in possibleIntersections)
        {

            Ellipse intersectionC = new Ellipse(majorDistance, minorDistance, intersectionVantagePlane, -_AB, 0);
            Ellipse circlePhi = new Ellipse(r, r, intersection_ABPlane, -_AB);
            circlePhi.Draw(Color.blue);
            intersectionC.Draw(Color.red);
        }
    }

    private Interval appropiatePhiIntervalBounds(Vector2 x, float samplingRate)
    {
        if (x.x == float.NaN || x.y == float.NaN)
        {
            throw new Exception("no intersection between the vantage");
        }


        //RES 1 Thesis
        float phiBoundsRes1 = Mathf.Atan(x.y / x.x);




        return new Interval(phiBoundsRes1, -phiBoundsRes1, samplingRate);


    }

    private Interval appropiateBetaIntervalBounds(Vector2 beta1, Vector2 beta2, Boolean greaterHalfPi)
    {

        //RES 2 Paper
        float sqrtBeta1 = Mathf.Sqrt(beta1.x * beta1.x + beta1.y * beta1.y);
        float sqrtBeta2 = Mathf.Sqrt(beta2.x * beta2.x + beta2.y * beta2.y);

        float lowerBeta = Mathf.Atan(sqrtBeta1);
        float upperBeta = Mathf.Atan(sqrtBeta2);

        if (greaterHalfPi)
        {
            lowerBeta = Mathf.PI - lowerBeta;
            upperBeta = Mathf.PI - upperBeta;
        }

        if (_targetPosition.Equals(B))
        {
            lowerBeta = Mathf.PI - lowerBeta;
            upperBeta = Mathf.PI - upperBeta;
        }

        return new Interval(lowerBeta, upperBeta);
    }




    //private methods


    private void DrawPlane(Vector3 position, Vector3 normal)
    {
        Vector3 v3;

        if (normal.normalized != Vector3.forward)
            v3 = Vector3.Cross(normal, Vector3.forward).normalized * normal.magnitude;
        else
            v3 = Vector3.Cross(normal, Vector3.up).normalized * normal.magnitude; ;

        var corner0 = position + v3;
        var corner2 = position - v3;
        var q = Quaternion.AngleAxis(90, normal);
        v3 = q * v3;
        var corner1 = position + v3;
        var corner3 = position - v3;

        Debug.DrawLine(corner0, corner2, Color.green, Mathf.Infinity);
        Debug.DrawLine(corner1, corner3, Color.green, Mathf.Infinity);
        Debug.DrawLine(corner0, corner1, Color.green, Mathf.Infinity);
        Debug.DrawLine(corner1, corner2, Color.green, Mathf.Infinity);
        Debug.DrawLine(corner2, corner3, Color.green, Mathf.Infinity);
        Debug.DrawLine(corner3, corner0, Color.green, Mathf.Infinity);
        Debug.DrawRay(position, normal, Color.red, Mathf.Infinity);
    }



    private int checkConicSection(float beta, Cone vantageCone)
    {
        switch (checkBetaForPlane(beta))
        {
            case -1:
                if (beta == 0) return -3; //circle
                if (beta + vantageCone._deviationAngle < Mathf.PI / 2) return -2; //ellipse
                if (beta + vantageCone._deviationAngle == Mathf.PI / 2) return -1; // parabola
                else return 0; //hyperbola
            case 0:
                if (beta + vantageCone._deviationAngle == Mathf.PI / 2) return 1;
                else return 2;
            case 1:
                float differentEquals = beta - vantageCone._deviationAngle - (Mathf.PI / 2);

                if (beta == Mathf.PI) return -3; //circle
                if (beta - vantageCone._deviationAngle > Mathf.PI / 2) return -2; // ellipse
                if ((differentEquals < 0.00001f) && (Mathf.Abs(differentEquals) <= 0.00001f)) return -1; //parabola
                else return 0; //hyperbola
            default:
                return 3;

        }

    }

    private int checkBetaForPlane(float beta)
    {
        if (beta < Mathf.PI / 2)
        {
            return -1;
        }
        else if (beta == Mathf.PI / 2)
        {
            return 0;
        }

        return 1;
    }

    private Vector2 solveQuadraticEquation(float A, float B, float C)
    {

        Vector2 res = new Vector2();
        float discriminant = Mathf.Sqrt((B * B) - 4 * A * C);

        if (discriminant > 0)
        {
            res[0] = (-B + discriminant) / (2 * A);
            res[1] = (-B - discriminant) / (2 * A);
        }
        else
        {
            res[0] = float.NaN;
            res[1] = float.NaN;
        }

        return sortVector2(res);

    }

    private Vector2 sortVector2(Vector2 toSort)
    {
        float temp = toSort[0];
        if (!(float.IsNaN(temp) || float.IsNaN(toSort[1])))
        {
            toSort[0] = Mathf.Min(temp, toSort[1]);
            toSort[1] = Mathf.Max(temp, toSort[1]);
        }
        else
        {
            toSort[0] = float.NaN;
            if (!float.IsNaN(temp)) toSort[1] = temp;
        }

        return toSort;
    }
}
