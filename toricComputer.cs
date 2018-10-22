
using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class ToricComputing
{

    private GameObject _target1;
    private GameObject _target2;


    //target positions
    private Vector3 A;
    private Vector3 B;
    private Vector3 AB;

    private float Sx;
    private float Sy;

    private float _radius;
    private float[] possiblePhiIntersect;
    Vector3 targetPosition;

    public readonly AlphaComputer _alphaComputer;
    public readonly vantageConstraint _vantageCons;


    public ToricComputing(GameObject target1, GameObject target2)
    {
        _target1 = target1;
        _target2 = target2;
        A = _target1.transform.position;
        B = _target2.transform.position;
        AB = B - A;
        _alphaComputer = new AlphaComputer(AB);
        _vantageCons = new vantageConstraint(A,B);
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

        float _radius = boundingSphereRadius;

        Vector2 SxAndSy = ComputeScale();
        Sx = SxAndSy[0];
        Sy = SxAndSy[1];

        float minDistance = ExactDistanceByProjectedSize(sizeConstraint[0], _radius);
        float maxDistance = ExactDistanceByProjectedSize(sizeConstraint[1], _radius);

        return new Interval(minDistance, maxDistance);
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
        Dictionary<float, Interval> phiBetaA = _vantageCons.getPositionFromVantageOneTarget(A, vantageA, deviationA, samplingRate);
        Dictionary<float, Interval> phiBetaB = _vantageCons.getPositionFromVantageOneTarget(B, vantageB, deviationB, samplingRate);


        float[] phiAKeys = phiBetaA.Keys.ToArray();
        float[] phiBKeys = phiBetaB.Keys.ToArray();

        Interval phiA = Interval.fromFloatArray(phiAKeys);
        Interval phiB = Interval.fromFloatArray(phiBKeys);

        Interval phiInv = phiA.Intersect(phiB);

        if (phiInv == null) throw new Exception("No intersection");
        possiblePhiIntersect = phiInv.filterArray(phiAKeys);


        foreach (float phi in possiblePhiIntersect)
        {
            Interval phiIntervallA, thetaIntervallA;
            phiBetaA.TryGetValue(phi, out phiIntervallA);
            thetaIntervallA = new Interval(2 * phiIntervallA.LOWERBOUND, 2 * phiIntervallA.UPPERBOUND);

            res.Add(phi, thetaIntervallA);

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
        Dictionary<float, Interval> AlphaTheta = _alphaComputer.getIntervalOfAcceptedAlpha(distanceToA, distanceToB, dTHETA);
        Dictionary<float, Interval> phiBetaB = _vantageCons.getPositionFromVantageOneTarget(B, vantageB, deviationB, dPHI);

        foreach (float phi in possiblePhiIntersect)
        {
            Interval thetaInvPhi, betaInvB;
            ThetaPhi.TryGetValue(phi, out thetaInvPhi);


            float[] possibleThetas = thetaInvPhi.filterArray(AlphaTheta.Keys.ToArray());
            phiBetaB.TryGetValue(mergeKey(phi, phiBetaB.Keys.ToArray(), dPHI), out betaInvB);

            foreach (float theta in possibleThetas)
            {
                Interval alphaDISTTheta, alphaOSP, alphaVANTThetaPhi;
                AlphaTheta.TryGetValue(theta, out alphaDISTTheta);

                alphaOSP = _alphaComputer.getAlphaIntervalFromOnscreenPositions(desPosA, desPosB);


                alphaVANTThetaPhi = _alphaComputer.GetVantageAlphaInterval(theta, betaInvB);


                Debug.Log("AlphaOnScreenPosition: " + alphaOSP * Mathf.Rad2Deg + "Alpha DistanceToTargets: " + alphaDISTTheta * Mathf.Rad2Deg + "Alpha Vantage Constraint: " + alphaVANTThetaPhi * Mathf.Rad2Deg);


                Interval alphaSecond = alphaOSP.Intersect(alphaDISTTheta);
                if (alphaSecond != null)
                {

                    Interval alphaFINAL = alphaSecond.Intersect(alphaVANTThetaPhi);


                    if (alphaFINAL != null)
                    {

                        alphaFINAL._samplingRate = dALPHA;
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
            tmVis.Add(new KeyValuePair<float, Toricmanifold>(visibility, tm));

        }

        Lookup<float, Toricmanifold> TableVisTm = (Lookup<float, Toricmanifold>)tmVis.ToLookup((item) => item.Key, (item) => item.Value);

        visibilityScores = new Interval(visibilityValues.x, visibilityValues.y).filterArray(TableVisTm.Select(g => g.Key).ToArray());

        float topVis = Mathf.Max(visibilityScores);

        List<KeyValuePair<float, Toricmanifold>> bestVisTm = tmVis.Where(g => g.Key == topVis).ToList();



        return bestVisTm.First().Value;
    }

    public static float visibilityCheck(Toricmanifold tm)
    {
        Vector3 origin = tm.ToWorldPosition();
        float visibilityScore = 0;

        //TODO check whether collider mesh or renderer bounds
        Bounds b = tm._target1.GetComponent<Renderer>().bounds;

        float visScoreTargetA = visibilityOneTarget(b, origin);

        b = tm._target2.GetComponent<Renderer>().bounds;

        float visScoreTargetB = visibilityOneTarget(b, origin);

        visibilityScore = (visScoreTargetA + visScoreTargetB) / 18;

        return visibilityScore;
    }

    private static float visibilityOneTarget(Bounds b, Vector3 origin)
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

    

    //private methods




    static private Vector2 ComputeScale()
    {
        Camera _main = Camera.main;

        float VerticalfovAngleRAD = _main.fieldOfView * Mathf.Deg2Rad;
        float tanY = Mathf.Tan(VerticalfovAngleRAD / 2);
        float tanX = tanY * _main.aspect;
        float Sx = 1 / tanX;
        float Sy = 1 / tanY;

        return new Vector2(Sx, Sy);
    }

    //Helper methods

    public static Vector3 FromWorldPosition(Vector3 CamPos, Vector3 target1, Vector3 target2)
    {
        Vector3 AB = target2-target1;
        Vector3 AP = CamPos-target1;

        Vector3 n = Vector3.Cross(AP,AB).normalized;
        Vector3 z;

        Vector2 n2 = new Vector2(AB.x, AB.z);

        float tmp = n2[0];
        n2[0] = n2[1];
        n2[1] = -tmp;

       z = new Vector3(n2.x,0 , n2.y).normalized;


        float beta = Vector3.Angle(AB, AP) * Mathf.Deg2Rad;
        float alpha = Vector3.Angle(AP,CamPos - target2) * Mathf.Deg2Rad;
        float theta = beta * 2 * Mathf.Deg2Rad;
        float phi = Vector3.SignedAngle(n, -AB, z) * Mathf.Deg2Rad - Mathf.PI / 2;

        return new Vector3(alpha, theta, phi);
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

    static public Vector3 GetVectorInCameraSpace(Vector2 vectorToCam)
    {
        Vector2 scaleFactors = ComputeScale();
        float Sx = scaleFactors[0];
        float Sy = scaleFactors[1];
        Vector3 vec;
        vec = new Vector3(vectorToCam.x / Sx, vectorToCam.y / Sy, 1);
        return vec;
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

   
}