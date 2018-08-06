
using System;
using System.Collections.Generic;
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

    public Dictionary<float, Intervall> getIntervalOfAcceptedAlpha(Vector2 distanceToA, Vector2 distanceToB)
    {
        float minDistanceToA = distanceToA[0];
        float maxDistanceToA = distanceToA[1];
        float minDistanceToB = distanceToB[0];
        float maxDistanceToB = distanceToB[1];

        Dictionary<float, Intervall> IaA = getIntervalFromA(minDistanceToA, maxDistanceToA);

        Dictionary<float, Intervall> IaB = getIntervalFromB(minDistanceToA, maxDistanceToA);
        Dictionary<float, Intervall> Ia = new Dictionary<float, Intervall>();


        foreach (float k in IaA.Keys)
        {
            Intervall alphasForThetaA;
            Intervall alphasForThetaB;
            IaA.TryGetValue(k, out alphasForThetaA);
            IaB.TryGetValue(k, out alphasForThetaB);

            Intervall alphaDouble = alphasForThetaA.Intersect(alphasForThetaB);

            Ia.Add(k, alphaDouble);
        }

        Ia = ClearDictionaryValues(Ia);

        return Ia;
    }

    private Dictionary<float, Intervall> ClearDictionaryValues(Dictionary<float, Intervall> ia)
    {
        List<float> keys = new List<float>(ia.Keys);
        foreach (float f in keys)
        {
            Intervall output;
            if (ia.TryGetValue(f, out output) && output == null)
            {
                ia.Remove(f);
            }
        }
        return ia;
    }

    /**
     * Calculating the alpha interval for each theta value based on the distance to target B
     * 
     * later used in getIntervalOfAcceptedAlpha()
     */
    public Dictionary<float, Intervall> getIntervalFromB(float minDistanceToB, float maxDistanceToB)
    {

        Intervall possibleThetaValues = new Intervall(1, 359, 1);
        Dictionary<float, Intervall> IaB = new Dictionary<float, Intervall>();

        foreach (float t in possibleThetaValues.getEveryValue())
        {

            Intervall alphaTestMin = GetAlphaFromDistanceToB(minDistanceToB, t);
            Intervall alphaTestMax = GetAlphaFromDistanceToB(maxDistanceToB, t);


            Intervall alphaInterval = alphaTestMin - alphaTestMax;

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
    private Dictionary<float, Intervall> getIntervalFromA(float minDistanceToA, float maxDistanceToA)
    {
        Intervall possibleThetaValues = new Intervall(1, 359, 1);
        Dictionary<float, Intervall> IaA = new Dictionary<float, Intervall>();

        foreach (float t in possibleThetaValues.getEveryValue())
        {
            float AlphaMinA = GetAlphaFromDistanceToA(minDistanceToA, t);
            float AlphaMaxA = GetAlphaFromDistanceToA(maxDistanceToA, t);
            AlphaMinA = Mathf.Clamp(AlphaMinA, 1, 359);
            AlphaMaxA = Mathf.Clamp(AlphaMaxA, 1, 359);
            Intervall alphaInterval = new Intervall(AlphaMinA, AlphaMaxA);
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
    public Intervall GetAlphaFromDistanceToB(float distance, float theta)
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

        return Intervall.fromFloatArray(res);
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
    public Intervall DistanceFromProjectedSize(Vector2 sizeConstraint, float boundingSphereRadius, GameObject target)
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

        return new Intervall(minDistance, maxDistance);
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
     * should later give possible values of theta and phi  
     * 
     * TODO Everything
     * 
     */
    public Intervall getPositionFromVantage(float whichOne, Vector3 v, float deviationAngle)
    {

        Vector3 targetPosition;
        if (whichOne == 1) targetPosition = A;
        else targetPosition = B;


        Vector3 normedAB = AB.normalized;

        Plane pBsmallerHalfPi = new Plane(AB, targetPosition - normedAB);
        Plane pBHalfPi = new Plane(AB, targetPosition);
        Plane pBGreaterHalfPi = new Plane(AB, targetPosition + normedAB);

        Vector3 prefferedVantageAngle = v;
        FixAngle beta = new FixAngle(Vector3.Angle(AB, prefferedVantageAngle), 180);
        Cone vantageCone = new Cone(prefferedVantageAngle, deviationAngle, targetPosition);
        float r = Mathf.Tan(beta.angle() * Mathf.Deg2Rad);

        Plane currentPlane = pBHalfPi;

        switch (checkBetaForPlane(beta.angle()))
        {
            case -1:
                currentPlane = pBsmallerHalfPi;
                break;
            case 0:
                currentPlane = pBHalfPi;
                break;
            case 1:
                currentPlane = pBGreaterHalfPi;
                break;
        }

        Vector2 x, y;
        switch (checkConicSection(beta.angle(), vantageCone))
        {
            case -2: //ellipse


                
                float sinDeviation = Mathf.Sin(deviationAngle * Mathf.Deg2Rad);
                float cosBeta = Mathf.Cos(beta.angle() * Mathf.Deg2Rad);
                float majorDistance = (sinDeviation * Mathf.Cos(deviationAngle * Mathf.Deg2Rad)) / (Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2)));
                float minorDistance = sinDeviation / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2))));
                float midPointDistance = (Mathf.Sin(beta.angle() * Mathf.Deg2Rad) * cosBeta) / Mathf.Pow(cosBeta, 2) - Mathf.Pow(sinDeviation, 2);
                Ellipse intersectionC = new Ellipse(majorDistance, minorDistance, midPointDistance * Vector2.up);

                
                float A = Mathf.Pow(minorDistance, 2) - Mathf.Pow(majorDistance, 2);
                float B = -2 * Mathf.Pow(minorDistance, 2) * midPointDistance;
                float C = Mathf.Pow(minorDistance, 2) * Mathf.Pow(midPointDistance, 2) + Mathf.Pow(majorDistance, 2) * (Mathf.Pow(r, 2) - Mathf.Pow(minorDistance, 2));
              
               
                x = solveQuadraticEquation(A, B, C);
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                //DEBUG
                Debug.Log("Intersection points: " + x + y);



                return appropiateIntervalBounds(new Vector2(x.x,y.x));
               
                

            case -1: //conic section = parabola

               
                float cotBeta = 1 / Mathf.Tan(beta.toRad());
                float h = (Mathf.Tan(beta.toRad()) - cotBeta) / 2;
                float f = cotBeta / 2;

                
                x = solveQuadraticEquation(-1, -4 * f, 4 * f * h + Mathf.Pow(r, 2));
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                return appropiateIntervalBounds(new Vector2(x.x, y.x));

            case 0: //conic section = hyperbola

                float sinDeviation2 = Mathf.Sin(deviationAngle * Mathf.Deg2Rad);
                float cosBeta2 = Mathf.Cos(beta.angle() * Mathf.Deg2Rad);
                float majorDistance2 = (sinDeviation2 * Mathf.Cos(deviationAngle * Mathf.Deg2Rad)) / (Mathf.Abs(Mathf.Pow(cosBeta2, 2) - Mathf.Pow(sinDeviation2, 2)));
                float minorDistance2 = sinDeviation2 / Mathf.Sqrt((Mathf.Abs(Mathf.Pow(cosBeta2, 2) - Mathf.Pow(sinDeviation2, 2))));
                float midPointDistance2 = (Mathf.Sin(beta.angle() * Mathf.Deg2Rad) * cosBeta2) / Mathf.Pow(cosBeta2, 2) - Mathf.Pow(sinDeviation2, 2);
                




                x = solveQuadraticEquation(Mathf.Pow(minorDistance2, 2) + Mathf.Pow(majorDistance2, 2), -2 * Mathf.Pow(minorDistance2, 2) * midPointDistance2, Mathf.Pow(minorDistance2, 2) * Mathf.Pow(midPointDistance2, 2) - Mathf.Pow(majorDistance2, 2) * (Mathf.Pow(r, 2) + Mathf.Pow(minorDistance2, 2)));
                y.x = Mathf.Sqrt(Mathf.Pow(r, 2) - (Mathf.Pow(x.x, 2)));
                y.y = -y.x;

                return appropiateIntervalBounds(new Vector2(x.x, y.x));
            case 1: //conic section = line
                //TODO
                return null;
            
            case 2: //conic section = two lines
                //TODO
                return null; 
            default:
                return null; 






        }
    }

    private Intervall appropiateIntervalBounds(Vector2 x)
    {
        //TODO chech for positive x and y values
        //RES 1 Thesis
        float phiBoundsRes1 = Mathf.Atan(x.x / x.y);
        //already theta?
        return new Intervall(phiBoundsRes1, -phiBoundsRes1);

        //RES 2 Paper
        //float beta = Mathf.Atan(Mathf.Sqrt(Mathf.Pow(x.x, 2) + Mathf.Pow(x.y, 2)));
        //float betaMinus = -beta;

        //return new Intervall(2 * Mathf.Min(beta, betaMinus),2 *  Mathf.Max(beta, betaMinus));
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
            case 3:
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

    public Intervall getAlphaIntervalFromOnscreenPositions(Vector2 desPosA, Vector2 desPosB)
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



        return new Intervall(Mathf.Min(alphaValues.ToArray()), Mathf.Max(alphaValues.ToArray()));
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

    private Vector2 solveQuadraticEquation(float A,float B, float C)
    {
        Vector2 res = new Vector2();
        res[0] = (-B + Mathf.Sqrt(Mathf.Pow(B, 2) - 4 * A * C)) / 2 * A;
        res[1] = (-B - Mathf.Sqrt(Mathf.Pow(B, 2) - 4 * A * C)) / 2 * A;

        return res;
    }



    //Test methods
    public float testDistanceFromA(float distance, float theta)
    {
        return GetAlphaFromDistanceToA(distance, theta);

    }
}