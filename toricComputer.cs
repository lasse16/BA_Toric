
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
    private float _frameRadius = 0.3f;
    private float[] possiblePhiIntersect;
    private bool vantageSet;
    Vector3 targetPosition;

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

    public Dictionary<float, Interval> getIntervalOfAcceptedAlpha(Vector2 distanceToA, Vector2 distanceToB, float samplingRate = 0.05f)
    {
        float minDistanceToA = distanceToA[0];
        float maxDistanceToA = distanceToA[1];
        float minDistanceToB = distanceToB[0];
        float maxDistanceToB = distanceToB[1];

        Dictionary<float, Interval> IaA = getIntervalFromA(minDistanceToA, maxDistanceToA,samplingRate);
        Dictionary<float, Interval> IaB = getIntervalFromB(minDistanceToB, maxDistanceToB,samplingRate);
        Dictionary<float, Interval> Ia = new Dictionary<float, Interval>();

        Interval keys = Interval.fromFloatArray(IaA.Keys.ToArray()).Intersect(Interval.fromFloatArray(IaB.Keys.ToArray()));
        

        foreach (float k in keys.filterArray(IaA.Keys.ToArray()))
        {
            
            
            Interval alphasForThetaA, alphasForThetaB;
         
            IaA.TryGetValue(k, out alphasForThetaA);
            float kAdjust = mergeKey(k, IaB.Keys.ToArray(), samplingRate);
            IaB.TryGetValue(kAdjust, out alphasForThetaB);

            Interval alphaDouble = alphasForThetaA.Intersect(alphasForThetaB);

            if(alphaDouble !=null) Ia.Add(k, alphaDouble);
        }
        return Ia;
    }



    /**
     * Calculating the alpha interval for each theta value based on the distance to target B
     * 
     * later used in getIntervalOfAcceptedAlpha()
     */
    public Dictionary<float, Interval> getIntervalFromB(float minDistanceToB, float maxDistanceToB,float samplingRate = 0.05f)
    {

        Dictionary<float, Interval> invDistanceMin = GetAlphaFromDistanceToB(minDistanceToB,samplingRate);
        Dictionary<float, Interval> invDistanceMax = GetAlphaFromDistanceToB(maxDistanceToB,samplingRate);
        Dictionary<float, Interval> invDistanceAB = GetAlphaFromDistanceToB(AB.magnitude +0.01f, samplingRate);

        Dictionary<float, Interval> thetaIntersect = new Dictionary<float, Interval>();

        float[] keysMin, KeysMax;
        keysMin = invDistanceMin.Keys.ToArray();
        KeysMax = invDistanceMax.Keys.ToArray();

        if (minDistanceToB <= AB.magnitude && maxDistanceToB > AB.magnitude) keysMin = KeysMax;

        foreach (float theta in KeysMax.Intersect(keysMin))
        {
            Interval alphaFromMin, alphaFromMax;
            invDistanceMin.TryGetValue(theta, out alphaFromMin);
            invDistanceMax.TryGetValue(theta, out alphaFromMax);

            if(alphaFromMin ==null) invDistanceAB.TryGetValue(theta, out alphaFromMin);

            Interval alphaInterval = new Interval(alphaFromMin.getLowerBound(), alphaFromMax.getLowerBound());

            thetaIntersect.Add(theta, alphaInterval);
        }

        return thetaIntersect;
    }

    /**
    * Calculating the alpha interval for each theta value based on the distance to target A
    * 
    * later used in main getIntervalOfAcceptedAlpha()
    * 
    */
    public Dictionary<float, Interval> getIntervalFromA(float minDistanceToA, float maxDistanceToA, float samplingRate = 0.05f)
    {
        Interval possibleThetaValues = new Interval(0.01f, Mathf.PI * 2,samplingRate);
        Dictionary<float, Interval> IaA = new Dictionary<float, Interval>();

        foreach (float t in possibleThetaValues.getEveryValue())
        {
            float AlphaMinA = GetAlphaFromDistanceToA(minDistanceToA, t);
            float AlphaMaxA = GetAlphaFromDistanceToA(maxDistanceToA, t);
            
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
    public Dictionary<float,Interval> GetAlphaFromDistanceToB(float distance, float samplingRate)
    {
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();
          
        Interval thetaValues = new Interval(samplingRate, Mathf.PI * 2 - samplingRate, samplingRate);

        Boolean smallerThanAB = distance <= AB.magnitude;

        if (smallerThanAB)
        {
            thetaValues = new Interval(samplingRate, 2 * Mathf.Asin(distance / AB.magnitude), samplingRate);
        }

        foreach (float theta in thetaValues.getEveryValue())
        {
            float[] alphaRes = new float[2];

            float acos = Mathf.Acos(AB.magnitude / distance * Mathf.Sin(theta / 2));

            alphaRes[0] = (Mathf.PI / 2 - acos);

            alphaRes[1] = (Mathf.PI / 2 + acos);


        

        if (smallerThanAB) alphaRes[1] = Mathf.PI - theta/2;

            res.Add(theta, Interval.fromFloatArray(alphaRes)); 
        }
        return res;
    }

    /**
  *    returns a alpha calculated from the distance to A and a theta input
  *    
  *    @param distance the distance for which to calculate alpha for
  *    @param theta the angle for which to calculate alpha for
  *    @return alpha for the specified distance and theta value 
  *     
  */
    public float GetAlphaFromDistanceToA(float distance, float theta)
    {
        return Mathf.Acos((distance - AB.magnitude * Mathf.Cos(theta / 2)) /
            (Mathf.Sqrt((distance * distance) + (AB.magnitude * AB.magnitude) 
            - 2 * AB.magnitude * distance * Mathf.Cos(theta / 2))));
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
    public Dictionary<float, Interval> getThetaIntervallFromVantageBothTargets(Vector3 vantageA, float deviationA, Vector3 vantageB, float deviationB, float samplingRate)
    {

        Dictionary<float, Interval> res = new Dictionary<float, Interval>();
        Dictionary<float, Interval> phiBetaA = getPositionFromVantageOneTarget(1, vantageA, deviationA,samplingRate);
        Dictionary<float, Interval> phiBetaB = getPositionFromVantageOneTarget(2, vantageB, deviationB,samplingRate);
        vantageSet = true;

        float[] phiAKeys = phiBetaA.Keys.ToArray();
        float[] phiBKeys = phiBetaB.Keys.ToArray();


        

        Interval phiA = Interval.fromFloatArray(phiAKeys);
        Interval phiB =Interval.fromFloatArray(phiBKeys);


        Debug.Log("phi A: " + phiA);
        Debug.Log("phi B: " + phiB);

        Debug.Log("phi intersect" + phiA.Intersect(phiB));

        Interval phiInv = phiA.Intersect(phiB);
        

        possiblePhiIntersect = phiInv.filterArray(phiAKeys);


        foreach (float phi in possiblePhiIntersect)
        {
            Interval phiIntervallA, thetaIntervallA;
            phiBetaA.TryGetValue(phi, out phiIntervallA);
            thetaIntervallA = new Interval(2 * phiIntervallA.getLowerBound() , 2 * phiIntervallA.getUpperBound());

            res.Add(phi, thetaIntervallA);

        }
        return res;
    }

    public Interval GetVantageAlphaInterval(float theta, Interval beta)
    {
        //TODO check for correct beta
        float alphaMin = Mathf.PI - theta / 2 - (beta.getLowerBound());
        float alphaMax = Mathf.PI - theta / 2 - (beta.getUpperBound());

        return new Interval(alphaMin, alphaMax);
    }

    private float CastBetaPrimeToTheta(float beta, float alpha)
    {
        return 2 * (Mathf.PI - alpha - beta);
    }

    public Interval getPhiInterval()
    {
        if (vantageSet)
        {
            return Interval.fromFloatArray(possiblePhiIntersect);
        }
        return new Interval(-Mathf.PI, Mathf.PI);

    }


    public Dictionary<float, Interval> getPositionFromVantageOneTarget(float whichOne, Vector3 v, float deviationAngle, float samplingRate = 0.01f)
    {
        targetPosition = B;
        if (whichOne == 1) targetPosition = A;

        if (v.Equals(Vector3.zero)) throw new Exception("No Vantage constraint set");
        if (deviationAngle > 90 ) throw new Exception("choose a deviattion angle smaller 90°");


        Vector3 upOnPlane, Vflat, PHIflat , vectorOnPlane;

        float phiVector = 0;
        
            
            
            Plane targetlevel = new Plane(Vector3.up, targetPosition);
            
            vectorOnPlane = targetlevel.ClosestPointOnPlane(targetPosition + v) - targetPosition;
           
            upOnPlane = Vector3.ProjectOnPlane(Vector3.up, AB);

            Vector3 phiZero = Vector3.Cross(upOnPlane,-AB);
            


            Vflat = Vector3.ProjectOnPlane(v.normalized, AB);
            PHIflat = Vector3.ProjectOnPlane(phiZero.normalized, AB);
            phiVector = Vector3.SignedAngle(Vflat, PHIflat, AB) * Mathf.Deg2Rad;

            
            v = vectorOnPlane;
        

        Vector3 prefferedVantageAngle = v.normalized;
        
        
        float signedLambda = Vector3.SignedAngle(-AB, prefferedVantageAngle, upOnPlane) * Mathf.Deg2Rad;
        
        float lambda = Mathf.Abs(signedLambda);
        deviationAngle *= Mathf.Deg2Rad;


        Vector2 x, y, y2;
        float a, b, c;
        Interval phiIntervall;
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();
        List<Vector2> possibleIntersections = new List<Vector2>();

        



        

        Boolean greaterThanHalfPi = lambda > Mathf.PI / 2;
        if (greaterThanHalfPi) lambda = Mathf.PI - lambda;



        Cone vantageCone = new Cone(prefferedVantageAngle, deviationAngle, targetPosition);
        float r = Mathf.Tan(lambda);



        switch (checkConicSection(lambda, vantageCone))
        {
            
            case -3: //circle
                Debug.Log("Circle");
                phiIntervall = new Interval(-Mathf.PI, Mathf.PI,samplingRate);

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
                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x),greaterThanHalfPi));
                    }
                    else
                    {
                        a = Mathf.Pow(minorDistance, 2) + Mathf.Pow(majorDistance, 2) * Mathf.Pow(Mathf.Tan(phi), 2);
                        b = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;
                        c = Mathf.Pow(minorDistance, 2) * (Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2));

                        x = solveQuadraticEquation(a, b, c);
                        



                        //TODO check for 1 or 2 intersection

                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.x)) x.x *= -1;
                        if (Mathf.Sign(Mathf.Cos(phi)) != Mathf.Sign(x.y)) x.y *= -1;
                        y.x = Mathf.Tan(phi) * x.x;
                        y.y = Mathf.Tan(phi) * x.y;


                        res.Add(phi, appropiateBetaIntervalBounds(new Vector2(x.y, y.y), new Vector2(x.x, y.x), greaterThanHalfPi));
                    }
                    VisualizeVantageConstraint(targetPosition, r, lambda, vantageCone, possibleIntersections, majorDistance, minorDistance);
                }

                break;

            case -1: //conic section = parabola
                Debug.Log("Parabola");

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
                b = -2 * Mathf.Pow(minorDistance, 2) * Mathf.Abs( midPointDistance);
                c = Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) - Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) + Mathf.Pow(minorDistance, 2));

                x = solveQuadraticEquation(a,b,c);

               
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
                    VisualizeVantageConstraint(targetPosition, r, lambda, vantageCone, possibleIntersections, majorDistance, minorDistance);
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


    
    private Interval appropiatePhiIntervalBounds(Vector2 x, float samplingRate)
    {
        if (x.x == float.NaN || x.y == float.NaN)
        {
            throw new Exception("no intersection between the vantage");
        }
       

        //RES 1 Thesis
        float phiBoundsRes1 = Mathf.Atan(x.y / x.x);
        
        
        

        return new Interval(phiBoundsRes1,-phiBoundsRes1, samplingRate);


    }

    private Interval appropiateBetaIntervalBounds(Vector2 beta1, Vector2 beta2, Boolean greaterHalfPi)
    {

        //RES 2 Paper
        float sqrtBeta1 = Mathf.Sqrt(beta1.x * beta1.x + beta1.y * beta1.y);
        float sqrtBeta2 = Mathf.Sqrt(beta2.x * beta2.x + beta2.y * beta2.y);

        float lowerBeta = Mathf.Atan(sqrtBeta1);
        float upperBeta = Mathf.Atan(sqrtBeta2);

        if(greaterHalfPi)
        {
            lowerBeta = Mathf.PI - lowerBeta;
            upperBeta = Mathf.PI - upperBeta;
        }

        if (targetPosition.Equals(B))
        {
            lowerBeta = Mathf.PI - lowerBeta;
            upperBeta = Mathf.PI - upperBeta;
        }

        return new Interval(lowerBeta, upperBeta);
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
                if (beta == 0) return -3; //circle
                if (beta + vantageCone.getDiviationAngle() < Mathf.PI/2) return -2; //ellipse
                if (beta + vantageCone.getDiviationAngle() == Mathf.PI/2) return -1; // parabola
                else return 0; //hyperbola
            case 0:
                if (beta + vantageCone.getDiviationAngle() == Mathf.PI/2) return 1;
                else return 2;
            case 1:
                float differentEquals = beta - vantageCone.getDiviationAngle()-(Mathf.PI / 2);

                if (beta == Mathf.PI) return -3; //circle
                if (beta - vantageCone.getDiviationAngle() > Mathf.PI/2) return -2; // ellipse
                if ((differentEquals< 0.00001f) && (Mathf.Abs(differentEquals)<= 0.00001f)) return -1; //parabola
                else return 0; //hyperbola
            default:
                return 3;

        }

    }


    /**
     * returns -1,0,1 for beta <pi/2,PI/2,>PI/2
     * 
     */
    private int checkBetaForPlane(float beta)
    {
        if (beta < Mathf.PI/2)
        {
            return -1;
        }
        else if (beta == Mathf.PI/2)
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
                float alpha = Vector3.Angle(pA, pB) * Mathf.Deg2Rad;
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


    public Toricmanifold FinalConstraintCombination(float samplingRateN, Vector3 vantageA, float deviationA, Vector3 vantageB, float deviationB, Vector2 distanceToA, Vector2 distanceToB, Vector2 desPosA, Vector2 desPosB, Vector2 visibilityValues)
    {
        samplingRateN = Mathf.Clamp(samplingRateN, 1, samplingRateN);
        float dPHI = 1f / (2 * Mathf.Pow(samplingRateN / 2, 1.0f / 3));
        float dTHETA = 1f / (4 * Mathf.Pow(Mathf.Pow(samplingRateN / 2, 2), 1.0f / 3));
        float dALPHA = 1f / (samplingRateN);

        List<Toricmanifold> possiblePositions = new List<Toricmanifold>();

        Dictionary<float, Interval> ThetaPhi = getThetaIntervallFromVantageBothTargets(vantageA, deviationA, vantageB, deviationB, dPHI);
        Dictionary<float, Interval> AlphaTheta = getIntervalOfAcceptedAlpha(distanceToA, distanceToB,dTHETA);
        Dictionary<float, Interval> phiBetaB = getPositionFromVantageOneTarget(2, vantageB, deviationB, dPHI);

        foreach (float phi in possiblePhiIntersect)
        {
            Interval thetaInvPhi, betaInvB;
            ThetaPhi.TryGetValue(phi, out thetaInvPhi);
            
            
            float[] possibleThetas = thetaInvPhi.filterArray(AlphaTheta.Keys.ToArray());
           phiBetaB.TryGetValue(mergeKey(phi,phiBetaB.Keys.ToArray(),dPHI), out betaInvB);

            foreach (float theta in possibleThetas)
            {
                Interval alphaDISTTheta, alphaOSP, alphaVANTThetaPhi;
                AlphaTheta.TryGetValue(theta, out alphaDISTTheta);
                
                alphaOSP = getAlphaIntervalFromOnscreenPositions(desPosA, desPosB);

               
                alphaVANTThetaPhi = GetVantageAlphaInterval(theta, betaInvB);


                Debug.Log("AlphaOnScreenPosition: " + alphaOSP + "Alpha DistanceToTargets: " + alphaDISTTheta + "Alpha Vantage Constraint: " + alphaVANTThetaPhi);


                Interval alphaSecond = alphaOSP.Intersect(alphaDISTTheta);
                if (alphaSecond != null)
                {

                    Interval alphaFINAL = alphaSecond.Intersect(alphaVANTThetaPhi);


                    if (alphaFINAL != null)
                    {

                        alphaFINAL.setSamplingRate(dALPHA);
                        foreach (float alpha in alphaFINAL.getEveryValue())
                        {
                            possiblePositions.Add(new Toricmanifold(alpha * Mathf.Rad2Deg, theta * Mathf.Rad2Deg, phi * Mathf.Rad2Deg, _target1, _target2));
                        }

                        Debug.Log("alphaFinal: " + alphaFINAL);

                    }
                }
            }
        }

        if (possiblePositions.Count == 0) throw new Exception("No solution possible");
        float[] visibilityScores;
        
        List<KeyValuePair<float, Toricmanifold>> tmVis = new List<KeyValuePair<float, Toricmanifold>>();
        foreach (Toricmanifold tm in possiblePositions)
        {
            float visibility = visibilityCheck(tm);
            tmVis.Add(new KeyValuePair<float, Toricmanifold>(visibility,tm));

        }

        Lookup<float, Toricmanifold> TableVisTm = (Lookup<float, Toricmanifold>) tmVis.ToLookup((item) => item.Key, (item) => item.Value);

        visibilityScores = new Interval(visibilityValues.x, visibilityValues.y).filterArray(TableVisTm.Select(g => g.Key).ToArray());

        float topVis = Mathf.Max(visibilityScores);

        List<KeyValuePair<float, Toricmanifold>> bestVisTm = tmVis.Where(g => g.Key == topVis).ToList();



        return bestVisTm.First().Value;
    }

    static public float visibilityCheck(Toricmanifold tm)
    {
        Vector3 origin = tm.ToWorldPosition();
        float visibilityScore = 0;

        //TODO check whether collider mesh or renderer bounds
        Bounds b = tm.getTarget1().GetComponent<Renderer>().bounds;

        float visScoreTargetA = RaycastIntersetionSingleTarget(b, origin);

        b = tm.getTarget2().GetComponent<Renderer>().bounds;

        float visScoreTargetB = RaycastIntersetionSingleTarget(b, origin);

        visibilityScore = (visScoreTargetA + visScoreTargetB) / 18;

        return visibilityScore;
    }

    static private float RaycastIntersetionSingleTarget(Bounds b, Vector3 origin)
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

    private void VisualizeVantageConstraint(Vector3 targetPosition, float r, float lambda, Cone vantageCone, List<Vector2> possibleIntersections, float majorDistance,float minorDistance)
    {

        float coneHeight = Mathf.Sqrt(1 + r * r);

        
        vantageCone.draw(coneHeight);
        Vector3 intersectionVantagePlane = targetPosition + vantageCone.getMidLineOfCone()* coneHeight;
        Debug.DrawLine(targetPosition, intersectionVantagePlane, Color.black, Mathf.Infinity);
        Debug.DrawLine(targetPosition, targetPosition + -AB, Color.black, Mathf.Infinity);

        Vector3 intersectionABPlane = targetPosition + AB.normalized * checkBetaForPlane(lambda);
        DrawPlane(intersectionABPlane, -AB.normalized);

        foreach (Vector2 item in possibleIntersections)
        {
           
            Ellipse intersectionC = new Ellipse(majorDistance, minorDistance, intersectionVantagePlane, -AB, 0);
            Ellipse circlePhi = new Ellipse(r, r, intersectionABPlane, -AB);
            circlePhi.draw(Color.blue);
            intersectionC.draw(Color.red);
            }
        }
    


    //Helper methods


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


    private float mergeKey(float key, float[] keys2, float samplingRate)
    {
        Interval possibleMerge = new Interval(key - samplingRate / 2, key + samplingRate / 2);
        foreach (float key2 in keys2)
        {
            if (possibleMerge.IsInside(key2)) return key2;
             
        }
        return -1;
    }


    //Test methods
    public float testDistanceFromA(float distance, float theta)
    {
        return GetAlphaFromDistanceToA(distance, theta);

    }

    public float testVisibility(Toricmanifold tm)
    {
        return visibilityCheck(tm);
    }
}