
using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class ToricComputing
{

    private GameObject _target1;
    private GameObject _target2;
    private Camera _main;

    //target positions
    private Vector3 A;
    private Vector3 B;
    private Vector3 AB;

    private float Sx;
    private float Sy;

    private float _radius;
    private float _frameRadius = 0.2f;
    private float[] possiblePhiIntersect;
    private bool vantageSet;

    public ToricComputing(GameObject target1, GameObject target2)
    {
        _target1 = target1;
        _target2 = target2;
        A = _target1.transform.position;
        B = _target2.transform.position;
        AB = A - B;

    }

    /**
  * Main method for calculating the accepted alpha values
  * 
  * @param distanceToA a vector2 in the form [min. Distance to A, max. DIstance to A]
  * @param distanceToB a vector2 in the form [min. Distance to B, max. DIstance to B]
  * 
  * @return a Dictionary with an interval of accepted alpha values for each theta key in RADIANS
  */

    public Dictionary<float, Interval> getIntervalOfAcceptedAlpha(Vector2 distanceToA, Vector2 distanceToB)
    {
        float minDistanceToA = distanceToA[0];
        float maxDistanceToA = distanceToA[1];
        float minDistanceToB = distanceToB[0];
        float maxDistanceToB = distanceToB[1];

        Dictionary<float, Interval> IaA = getIntervalFromA(minDistanceToA, maxDistanceToA);

        Dictionary<float, Interval> IaB = getIntervalFromB(minDistanceToA, maxDistanceToA);
        Dictionary<float, Interval> Ia = new Dictionary<float, Interval>();


        foreach (float k in IaA.Keys)
        {
            Interval alphasForThetaA;
            Interval alphasForThetaB;
            IaA.TryGetValue(k, out alphasForThetaA);
            IaB.TryGetValue(k, out alphasForThetaB);

            Interval alphaDouble = alphasForThetaA.Intersect(alphasForThetaB);

            Ia.Add(k, alphaDouble);
        }

        Ia = ClearDictionaryValues(Ia);

        return Ia;
    }

  

    /**
     * Calculating the alpha interval for each theta value based on the distance to target B
     * 
     * later used in getIntervalOfAcceptedAlpha()
     */
    public Dictionary<float, Interval> getIntervalFromB(float minDistanceToB, float maxDistanceToB)
    {

        Interval possibleThetaValues = new Interval(1, 359, 1);
        Dictionary<float, Interval> IaB = new Dictionary<float, Interval>();

        foreach (float t in possibleThetaValues.getEveryValue())
        {

            Interval alphaTestMin = GetAlphaFromDistanceToB(minDistanceToB, t);
            Interval alphaTestMax = GetAlphaFromDistanceToB(maxDistanceToB, t);


            Interval alphaInterval = alphaTestMin - alphaTestMax;

            IaB.Add(t, alphaInterval);
        }

        return IaB;
    }

    /**
    * Calculating the alpha interval for each theta value based on the distance to target A
    * 
    * later used in main getIntervalOfAcceptedAlpha()
    * 
    */
    private Dictionary<float, Interval> getIntervalFromA(float minDistanceToA, float maxDistanceToA)
    {
        Interval possibleThetaValues = new Interval(1, 359, 1);
        Dictionary<float, Interval> IaA = new Dictionary<float, Interval>();

        foreach (float t in possibleThetaValues.getEveryValue())
        {
            float AlphaMinA = GetAlphaFromDistanceToA(minDistanceToA, t);
            float AlphaMaxA = GetAlphaFromDistanceToA(maxDistanceToA, t);
            AlphaMinA = Mathf.Clamp(AlphaMinA, 1, 359);
            AlphaMaxA = Mathf.Clamp(AlphaMaxA, 1, 359);
            Interval alphaInterval = new Interval(AlphaMinA, AlphaMaxA);
            IaA.Add(t, alphaInterval);
        }
        return IaA;
    }

    /**
 *    Returns an pair of theta/alpha for a distance to  target B and a theta input
 *    @param distance the distance for which to calculate alpha for
 *    @param theta the angle for which to calculate alpha for
 *    @return alpha for the specified distance and theta value
 *    
 *     
 *    //TODO find error in formula     
 *    
 */
    public Interval GetAlphaFromDistanceToB(float distance, float theta)
    {
        float[] res = new float[2];


        if (distance <= AB.magnitude)
        {
            //TODO theta = 0 after clamp?
            theta = Mathf.Clamp(theta, 1, 2 * Mathf.Asin(distance / AB.magnitude) * Mathf.Rad2Deg);
            double acosTest = Mathf.Clamp(AB.magnitude / distance * Mathf.Sin(theta * Mathf.Deg2Rad / 2), -1, 1);
            double acos = 0;
            if (acosTest != 0.0)
            {
                acos = Math.Acos(acosTest);
            }

            double alphaMINUS = Mathf.PI / 2 - acos;
            res[0] = (float)(alphaMINUS * Mathf.Rad2Deg);
            double alphaPLUS = Mathf.PI / 2 + acos;
            res[1] = (float)(alphaPLUS * Mathf.Rad2Deg);
        }
        else
        {
            float acosTest = Mathf.Clamp(AB.magnitude / distance * Mathf.Sin(theta * Mathf.Deg2Rad / 2), -1, 1);
            float acos = 0;
            if (acosTest != 0) acos = Mathf.Acos(acosTest);
            float alpha = Mathf.PI / 2 - acos;
            res[0] = alpha * Mathf.Rad2Deg;
            res[1] = (Mathf.PI - (theta * Mathf.Deg2Rad / 2)) * Mathf.Rad2Deg;
        }

        return Interval.fromFloatArray(res);
    }

    /**
  *    returns a theta/alpha pair calculated from the distance to A and a theta input
  *    
  *    @param distance the distance for which to calculate alpha for
  *    @param theta the angle for which to calculate alpha for
  *    @return alpha for the specified distance and theta value DEGREES
  *     
  */
    public float GetAlphaFromDistanceToA(float distance, float theta)
    {
        FixAngle thetaFix = new FixAngle(theta);
        float top1 = AB.magnitude * Mathf.Cos(thetaFix.toRad() / 2);
        float top = distance - top1;
        float bottom1 = top1 * distance * 2;
        float bottom2 = Mathf.Pow(AB.magnitude, 2) + Mathf.Pow(distance, 2);
        float bottom = Mathf.Sqrt(bottom2 - bottom1);

        return (Mathf.Acos(top / bottom) * Mathf.Rad2Deg);
    }

    public static string FloatArrayToString(float[] array)
    {
        string res = "";
        foreach (float f in array)
        {
            res = String.Concat(res, ";" + f);
        }
        return res;
    }
    /**
     * Calculates a distance interval given a Vector of a projected size interval and the radius of the target, approximated by a sphere.
     * @param sizeConstraint Vector2 in the form [min Size, max Size]
     * @param boundingSphereRadius the radius of the enclosing sphere
     * @return an interval of distance to the target
     * 
     */
    public Interval DistanceFromProjectedSize(Vector2 sizeConstraint, float boundingSphereRadius, GameObject target)
    {
        Camera.main.transform.LookAt(target.transform);

        float minSize = sizeConstraint[0];
        float maxSize = sizeConstraint[1];
        float _radius = boundingSphereRadius;

        Vector2 SxAndSy = ComputeScale();
        Sx = SxAndSy[0];
        Sy = SxAndSy[1];

        float minDistance = ExactDistanceByProjectedSize(minSize, _radius);
        float maxDistance = ExactDistanceByProjectedSize(maxSize, _radius);

        return new Interval(minDistance, maxDistance);
    }

    /**
 * helps compute qTrans depending on the main cameras field of view
 * https://gist.github.com/coastwise/5951291 (to get the horizontal fov from the vertical)
 */
    public static Vector2 ComputeScale()
    {
        Camera _main = Camera.main;
        /**
         * my version
        
        float VerticalfovAngleRAD = _main.fieldOfView * Mathf.Deg2Rad;

        float HorizontalfovAngleRAD = Mathf.Atan(Mathf.Tan(VerticalfovAngleRAD / 2) * _main.aspect);
        float Sx = 1 / Mathf.Tan(HorizontalfovAngleRAD);
        float Sy = 1 / Mathf.Tan(VerticalfovAngleRAD / 2);
        */

        //lino version
        float VerticalfovAngleRAD = _main.fieldOfView * Mathf.Deg2Rad;
        float tanY = Mathf.Tan(VerticalfovAngleRAD / 2);
        float tanX = tanY * _main.aspect;
        float Sx = 1 / tanX;
        float Sy = 1 / tanY;

        return new Vector2(Sx, Sy);
    }

    private float ExactDistanceByProjectedSize(float sizeToReach, float radius)
    {
        float radiand = (Mathf.PI * Sx * Sy) / (4 * sizeToReach);
        return radius * Mathf.Sqrt(radiand);
    }


    /**
     * 
     * returns the interval of possible theta alpha intersection based on phi
     * i.e. a theta interval in which the intersection is
     * 
     */
    public Dictionary<float, Interval> getThetaIntervallFromVantageBothTargets(Vector3 vantageA, float deviationA, Vector3 vantageB, float deviationB)
    {
        
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();
        Dictionary<float, Interval> phiBetaA = getPositionFromVantageOneTarget(1, vantageA, deviationA);
        Dictionary<float, Interval> phiBetaB = getPositionFromVantageOneTarget(2, vantageB, deviationB);
        vantageSet = true;

        float[] phiAKeys = phiBetaA.Keys.ToArray();
        float[] phiBKeys = phiBetaB.Keys.ToArray();



        possiblePhiIntersect = phiAKeys.Intersect(phiBKeys).ToArray();
        foreach (float phi in possiblePhiIntersect)
        {
            Interval phiIntervallA, thetaIntervallA;
            phiBetaA.TryGetValue(phi, out phiIntervallA);
            thetaIntervallA = new Interval(2 * phiIntervallA.getLowerBound(), phiIntervallA.getUpperBound() * 2);


            res.Add(phi, thetaIntervallA);

        }
        return res;
    }

    public Interval GetVantageAlphaInterval(float theta, Interval beta)
    {
       float alphaMin = Mathf.PI - theta / 2 - beta.getLowerBound();
       float alphaMax = Mathf.PI - theta / 2 - beta.getUpperBound();

        return new Interval(alphaMin, alphaMin);
    }

    private float CastBetaPrimeToTheta(float beta, float alpha)
    {
        return 2 * (Mathf.PI - alpha - beta);
    }

    public Interval getPhiInterval()
    {
        if (vantageSet)
        {
            return new Interval(Mathf.Min(possiblePhiIntersect) * Mathf.Rad2Deg, Mathf.Max(possiblePhiIntersect) * Mathf.Rad2Deg);
        }
        return new Interval(-180, 180);
        
    }


    /**
     * Should give an interval of beta for phi as a dictionary for all possible phi values
     * 
     * TODO Debug
     * 
     */
    public Dictionary<float, Interval> getPositionFromVantageOneTarget(float whichOne, Vector3 v, float deviationAngle)
    {
        Vector3 prefferedVantageAngle = v;
        FixAngle beta = new FixAngle(Vector3.Angle(AB, prefferedVantageAngle), 180);
        deviationAngle *= Mathf.Deg2Rad;

        if (deviationAngle > Mathf.PI / 2)
        {
            deviationAngle = Mathf.PI - deviationAngle;
            beta = new FixAngle(Mathf.PI - beta.angle());
        }
        //TODO check influence of early beta change on later functions || use exclusion /limit deviation angle

        Vector3 targetPosition;
        if (whichOne == 1) targetPosition = A;
        else targetPosition = B;




        Cone vantageCone = new Cone(prefferedVantageAngle, deviationAngle, targetPosition);
        float r = Mathf.Tan(beta.angle() * Mathf.Deg2Rad);





        Vector2 x, y;
        float a, b, c;
        Interval phiIntervall;
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();

        switch (checkConicSection(beta.angle(), vantageCone))
        {
            //TODO circle case
            case -2: //ellipse



                float sinDeviation = Mathf.Sin(deviationAngle * Mathf.Deg2Rad);
                float cosBeta = Mathf.Cos(beta.angle() * Mathf.Deg2Rad);
                float majorDistance = (sinDeviation * Mathf.Cos(deviationAngle * Mathf.Deg2Rad)) / (Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2)));
                float minorDistance = sinDeviation / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2))));
                float midPointDistance = (Mathf.Sin(beta.angle() * Mathf.Deg2Rad) * cosBeta) / Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2);
                Ellipse intersectionC = new Ellipse(majorDistance, minorDistance, midPointDistance * Vector2.up);


                a = Mathf.Pow(minorDistance, 2) - Mathf.Pow(majorDistance, 2);
                b = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;
                c = Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) + Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) - Mathf.Pow(minorDistance, 2));


                x = solveQuadraticEquation(a, b, c);
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                //DEBUG
                Debug.Log("Intersection points: " + x + y);



                phiIntervall = appropiatePhiIntervalBounds(new Vector2(x.x, y.x));

                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * (minorDistance / majorDistance) * Mathf.Sqrt(Mathf.Pow(majorDistance, 2) - Mathf.Pow(midPointDistance, 2));
                        y.y = Mathf.Sin(phi) * (minorDistance / majorDistance) * -Mathf.Sqrt(Mathf.Pow(majorDistance, 2) - Mathf.Pow(midPointDistance, 2));
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                    else
                    {
                        a = Mathf.Pow(minorDistance, 2) + Mathf.Pow(majorDistance, 2) * Mathf.Pow(Mathf.Tan(phi), 2);
                        c = Mathf.Pow(minorDistance, 2) * (Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2));
                        x = solveQuadraticEquation(a, b, c);

                        y.x = Mathf.Tan(phi) * x.x;

                        //TODO check for 1 or 2 intersection
                        y.y = Mathf.Tan(phi) * x.y;
                        //DEBUG
                        Debug.Log(y);

                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                }

                break;

            case -1: //conic section = parabola


                float cotBeta = 1 / Mathf.Tan(beta.toRad());
                float h = (Mathf.Tan(beta.toRad()) - cotBeta) / 2;
                float f = cotBeta / 2;


                x = solveQuadraticEquation(-1, -4 * f, 4 * f * h + Mathf.Pow(r, 2));
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                phiIntervall = appropiatePhiIntervalBounds(new Vector2(x.x, y.x));

                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * 2 * Mathf.Sqrt(-f * h);
                        y.y = -y.x;
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                    else
                    {
                        a = Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -4 * f;
                        c = 4 * f * h;
                        x = solveQuadraticEquation(a, b, c);

                        y.x = Mathf.Tan(phi) * x.x;

                        //TODO check for 1 or 2 intersection
                        y.y = Mathf.Tan(phi) * x.y;
                        //DEBUG
                        Debug.Log(y);

                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                }

                break;

            case 0: //conic section = hyperbola

                sinDeviation = Mathf.Sin(deviationAngle * Mathf.Deg2Rad);
                cosBeta = Mathf.Cos(beta.angle() * Mathf.Deg2Rad);
                majorDistance = (sinDeviation * Mathf.Cos(deviationAngle * Mathf.Deg2Rad)) / (Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2)));
                minorDistance = sinDeviation / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2))));
                midPointDistance = (Mathf.Sin(beta.angle() * Mathf.Deg2Rad) * cosBeta) / Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2);





                x = solveQuadraticEquation(Mathf.Pow(minorDistance, 2) + Mathf.Pow(majorDistance, 2), -2 * Mathf.Pow(minorDistance, 2) * midPointDistance, Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) + Mathf.Pow(minorDistance, 2)));
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                phiIntervall = appropiatePhiIntervalBounds(new Vector2(x.x, y.x));

                foreach (float phi in phiIntervall.getEveryValue())
                {
                    if (Mathf.Abs(phi) == Mathf.PI / 2)
                    {
                        x = Vector2.zero;
                        y.x = Mathf.Sin(phi) * (minorDistance / majorDistance) * Mathf.Sqrt(Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2)); ;
                        y.y = -y.x;
                        y = sortVector2(y);
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                    else
                    {
                        a = Mathf.Pow(minorDistance, 2) - Mathf.Pow(majorDistance, 2) * Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;
                        c = Mathf.Pow(minorDistance, 2) * (Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2));
                        x = solveQuadraticEquation(a, b, c);

                        y.x = Mathf.Tan(phi) * x.x;

                        //TODO check for 1 or 2 intersection
                        y.y = Mathf.Tan(phi) * x.y;
                        //DEBUG
                        Debug.Log(y);

                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y)));
                    }
                }

                break;

            case 1: //conic section = line
                res.Add(0, new Interval(Mathf.PI / 2, Mathf.PI / 2));
                break;

            case 2: //conic section = two lines
                //TODO check if correct||max Theta value
                phiIntervall = new Interval(-deviationAngle, deviationAngle);
                foreach (float phi in phiIntervall.getEveryValue())
                {
                    res.Add(phi, new Interval(0, 360));
                }
                break;
            default:
                throw new Exception();

        }

        return res;
    }

    private Interval appropiatePhiIntervalBounds(Vector2 x)
    {
        //TODO chech for positive x and y values | radians or degrees

        //RES 1 Thesis
        float phiBoundsRes1 = Mathf.Atan(x.x / x.y);

        return new Interval(phiBoundsRes1, -phiBoundsRes1);


    }

    private Interval appropiateBetaIntervalBounds(Vector2 x)
    {

        //TODO proper B* cast
        //RES 2 Paper
        float beta = Mathf.Atan(Mathf.Sqrt(Mathf.Pow(x.x, 2) + Mathf.Pow(x.y, 2)));
        float betaMinus = -beta;
        Vector2 betaSort = sortVector2(new Vector2(beta, betaMinus));
        return new Interval(betaSort[0], betaSort[1]);
    }


    /**
    * returns ,-2-1,0,1,2 for ellipse,parabola,hyperbola,line, two lines
    * Lino 2013 p.158
    */
    private int checkConicSection(float beta, Cone vantageCone)
    {
        switch (checkBetaForPlane(beta))
        {
            case -1:
                if (beta + vantageCone.getDiviationAngle() < 180) return -2; //ellipse
                if (beta + vantageCone.getDiviationAngle() == 180) return -1; // parabola
                else return 0; //hyperbola
            case 0:
                if (beta + vantageCone.getDiviationAngle() == 180) return 1;
                else return 2;
            case 1:
                if (beta - vantageCone.getDiviationAngle() < 180) return -2; //hyperbola
                if (beta - vantageCone.getDiviationAngle() == 180) return -1; //parabola
                else return 0; // ellipse
            default:
                return -3;

        }

    }


    /**
     * returns -1,0,1 for beta <90,90,>90
     * 
     */
    private int checkBetaForPlane(float beta)
    {
        if (beta < 90)
        {
            return -1;
        }
        else if (beta == 90)
        {
            return 0;
        }

        return 1;
    }


    /**
     * Returns an interval of alpha values based on frames
     * inside which the targets can be projected.
     * 
     * @param desPosA the desired  onscreen position of target A around the frame is going to be constructed
     * @param desPosB the desired  onscreen position of target B around the frame is going to be constructed
     * @return a interval of accepted alpha values in degrees
     * 
     */

    public Interval getAlphaIntervalFromOnscreenPositions(Vector2 desPosA, Vector2 desPosB)
    {
        List<Vector2> verticesFrameTargetA = getFrameTarget(desPosA);
        List<Vector2> verticesFrameTargetB = getFrameTarget(desPosB);



        List<Vector3> VerticesInCamSpaceA = new List<Vector3>();
        List<Vector3> VerticesInCamSpaceB = new List<Vector3>();
        foreach (Vector2 vec in verticesFrameTargetA)
        {
            Vector3 inCamSpace = GetVectorInCameraSpace(vec);
            VerticesInCamSpaceA.Add(inCamSpace);

        }
        foreach (Vector2 vec in verticesFrameTargetB)
        {
            Vector3 inCamSpace = GetVectorInCameraSpace(vec);
            VerticesInCamSpaceB.Add(inCamSpace);
        }



        List<float> alphaValues = new List<float>();

        foreach (Vector3 pA in VerticesInCamSpaceA)
        {

            foreach (Vector3 pB in VerticesInCamSpaceB)
            {
                float alpha = Vector3.Angle(pA, pB);
                alphaValues.Add(alpha);
            }
        }



        return new Interval(Mathf.Min(alphaValues.ToArray()), Mathf.Max(alphaValues.ToArray()));
    }


    public static Vector3 GetVectorInCameraSpace(Vector2 vectorToCam)
    {
        Vector2 scaleFactors = ComputeScale();
        float Sx = scaleFactors[0];
        float Sy = scaleFactors[1];
        Vector3 vec;
        vec = new Vector3(vectorToCam.x / Sx, vectorToCam.y / Sy, 1);
        return vec;
    }

    //creates a list of vertices belonging to an octagon around a target 
    private List<Vector2> getFrameTarget(Vector2 desPosTarget)
    {
        List<Vector2> res = new List<Vector2>();
        for (int i = 0; i < 8; i++)
        {
            float angle = 2 * Mathf.PI / 8;
            Vector2 axis = rotateVector2(Vector2.up, angle * i);
            Vector2 onscreenVertex = desPosTarget + axis * _frameRadius;

            res.Add(new Vector2(Mathf.Clamp(onscreenVertex.x, -1, 1), Mathf.Clamp(onscreenVertex.y, -1, 1)));
        }

        return res;
    }

    /// <summary>
    /// returns null if no closed form solution is found
    /// </summary>
    /// <param name="samplingRateN"></param>
    /// <param name="vantageA"></param>
    /// <param name="deviationA"></param>
    /// <param name="vantageB"></param>
    /// <param name="deviationB"></param>
    /// <param name="distanceToA"></param>
    /// <param name="distanceToB"></param>
    /// <param name="desPosA"></param>
    /// <param name="desPosB"></param>
    /// <param name="visibilityValues"></param>
    /// <returns></returns>
    public Toricmanifold FinalConstraintCombination(float samplingRateN, Vector3 vantageA,float deviationA, Vector3 vantageB, float deviationB, Vector2 distanceToA, Vector2 distanceToB, Vector2 desPosA, Vector2 desPosB, Vector2 visibilityValues)
    {
        List<Toricmanifold> possiblePositions = new List<Toricmanifold>();
        Toricmanifold res;
        Dictionary<float,Interval> ThetaPhi = getThetaIntervallFromVantageBothTargets(vantageA,deviationA,vantageB,deviationB);
        Dictionary<float, Interval> AlphaTheta = getIntervalOfAcceptedAlpha(distanceToA,distanceToB);
        Dictionary<float, Interval> phiBetaB = getPositionFromVantageOneTarget(2, vantageB, deviationB);
        Interval phiInterval = getPhiInterval();

        float dPHI = 2 * Mathf.Pow(samplingRateN / 2,   1.0f / 3);
        float dTHETA = 4 * Mathf.Pow(Mathf.Pow(samplingRateN / 2,2), 1.0f / 3);
        float dALPHA = samplingRateN;

        phiInterval.setSamplingRate(dPHI);

        foreach (float phi in phiInterval.getEveryValue())
        {
            Interval thetaInvPhi;
            ThetaPhi.TryGetValue(phi, out thetaInvPhi);

            thetaInvPhi.setSamplingRate(dTHETA);

            foreach (float theta in thetaInvPhi.getEveryValue())
            {
                Interval alphaDISTTheta, alphaOSP,alphaVANTThetaPhi;
                AlphaTheta.TryGetValue(theta, out alphaDISTTheta);
                alphaOSP = getAlphaIntervalFromOnscreenPositions(desPosA,desPosB);
                Interval betaInvB;
                phiBetaB.TryGetValue(phi, out betaInvB);
                alphaVANTThetaPhi =  GetVantageAlphaInterval(theta, betaInvB);
                Interval alphaFINAL = alphaOSP.Intersect(alphaDISTTheta.Intersect(alphaVANTThetaPhi));

                alphaFINAL.setSamplingRate(samplingRateN);

                foreach (float alpha in alphaFINAL.getEveryValue())
                {
                    possiblePositions.Add(new Toricmanifold(alpha, theta, phi, _target1, _target2));
                }

            }
        }

        float[] visibilityScores;
        Dictionary<float, Toricmanifold> tmVis = new Dictionary<float, Toricmanifold>(); 
        foreach (Toricmanifold tm in possiblePositions)
        {
            float visibility = visibilityCheck(tm);
            tmVis.Add(visibility, tm);
        }
        visibilityScores = sliceArrayAtValue(tmVis.Keys.ToArray(), visibilityValues);
        
        float topVis = Mathf.Max(visibilityScores);

        tmVis.TryGetValue(topVis, out res);

        return res;
    }

    private float visibilityCheck(Toricmanifold tm)
    {
        Vector3 origin = tm.ToWorldPosition();
        float visibilityScore = 0;

        //TODO check whether collider mesh or renderer bounds
        Bounds b = _target1.GetComponent<Renderer>().bounds;

        float visScoreTargetA = RaycastIntersetionSingleTarget(b, origin);
        
        b = _target2.GetComponent<Renderer>().bounds;

        float visScoreTargetB = RaycastIntersetionSingleTarget(b, origin);

        visibilityScore = (visScoreTargetA + visScoreTargetB) / 18;

        return visibilityScore;
    }

    private float RaycastIntersetionSingleTarget(Bounds b, Vector3 origin)
    {
        float visibilityScore = 0;
        Vector3 vertice1 = (b.center + new Vector3(-b.size.x, -b.size.y, -b.size.z) * 0.5f);
        Vector3 vertice2 = (b.center + new Vector3(b.size.x, -b.size.y, -b.size.z) * 0.5f);
        Vector3 vertice3 = (b.center + new Vector3(b.size.x, -b.size.y, b.size.z) * 0.5f);
        Vector3 vertice4 = (b.center + new Vector3(-b.size.x, -b.size.y, b.size.z) * 0.5f);
        Vector3 vertice5 = (b.center + new Vector3(-b.size.x, b.size.y, -b.size.z) * 0.5f);
        Vector3 vertice6 = (b.center + new Vector3(b.size.x, b.size.y, -b.size.z) * 0.5f);
        Vector3 vertice7 = (b.center + new Vector3(b.size.x, b.size.y, b.size.z) * 0.5f);
        Vector3 vertice8 = (b.center + new Vector3(-b.size.x, b.size.y, b.size.z) * 0.5f);

        Vector3[] verticesTarget = { b.center, vertice1, vertice2, vertice3, vertice4, vertice5, vertice6, vertice7, vertice8 };

        foreach (Vector3 vertex in verticesTarget)
        {
            if (!Physics.Linecast(origin, vertex)) visibilityScore++;
        }

        return visibilityScore;
    }

    //Helper methods

    private float[] sliceArrayAtValue(float[] array , Vector2 cutValues){
        int priorLength = array.Length;   
        array = new float[priorLength + 2];
        array[priorLength] = cutValues.x;
        array[priorLength+1] = cutValues.y;
        Array.Sort(array);
        int lowerCut = Array.IndexOf(array,cutValues.x);
        int upperCut = Array.IndexOf(array, cutValues.y);
        float[] res = new float[upperCut - lowerCut];
        int y = 0;
        for (int i = lowerCut + y; lowerCut +y< upperCut; y++)
        {
            res[y] = array[lowerCut + y];
        }

        return res;
    }


    private Vector2 rotateVector2(Vector2 v, float angle)
    {
        float sin = Mathf.Sin(angle);
        float cos = Mathf.Cos(angle);

        float tx = v.x;
        float ty = v.y;
        v.x = (cos * tx) - (sin * ty);
        v.y = (sin * tx) + (cos * ty);
        return v;
    }

    private Vector2 solveQuadraticEquation(float A, float B, float C)
    {
        Vector2 res = new Vector2();
        res[0] = (-B + Mathf.Sqrt(Mathf.Pow(B, 2) - 4 * A * C)) / 2 * A;
        res[1] = (-B - Mathf.Sqrt(Mathf.Pow(B, 2) - 4 * A * C)) / 2 * A;
        return sortVector2(res);
    }

    private Vector2 sortVector2(Vector2 toSort)
    {
        float temp = toSort[0];
        toSort[0] = Mathf.Min(temp, toSort[1]);
        toSort[1] = Mathf.Max(temp, toSort[1]);
        return toSort;
    }

    private Dictionary<float, Interval> ClearDictionaryValues(Dictionary<float, Interval> ia)
    {
        List<float> keys = new List<float>(ia.Keys);
        foreach (float f in keys)
        {
            Interval output;
            if (ia.TryGetValue(f, out output) && output == null)
            {
                ia.Remove(f);
            }
        }
        return ia;
    }

    //Test methods
    public float testDistanceFromA(float distance, float theta)
    {
        return GetAlphaFromDistanceToA(distance, theta);

    }
}