using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;


/// <summary>
/// Responsible for everything regarding alpha inside the toric manifold constraints
/// </summary>
public class AlphaComputer
{
    private Vector3 _AB;
    private Vector3 _A;
    private Vector3 _B;
    private float _frameRadius = 0.5f;

    public AlphaComputer(Vector3 A, Vector3 B)
    {
        _A = A;
        _B = B;
        _AB = B-A;
    }


    /// <summary>
    /// Creates the interval of accepted alpha given by the distance constrints on both targets
    /// </summary>
    /// <param name="distanceToA"> Vector2 [minimal distance; max distance] to A</param>
    /// <param name="distanceToB">Vector2 [minimal distance; max distance] to B</param>
    /// <param name="samplingRate">sampling rate for the resulting dictionary keys</param>
    /// <returns>an interval of accepted alpha</returns>
    public Dictionary<float, Interval> getIntervalOfAcceptedAlpha(Vector2 distanceToA, Vector2 distanceToB, float samplingRate = 0.05f, Boolean visualize = false)
    {
        if (visualize)
        {
            _visualizeDistance(distanceToA, distanceToB);
        }



        Dictionary<float, Interval> IaA = getIntervalFromA(distanceToA[0], distanceToA[1], samplingRate);
        Dictionary<float, Interval> IaB = getIntervalFromB(distanceToB[0], distanceToB[1], samplingRate);
        Dictionary<float, Interval> Ia = new Dictionary<float, Interval>();

        Interval keys = Interval.fromFloatArray(IaA.Keys.ToArray()).Intersect(Interval.fromFloatArray(IaB.Keys.ToArray()));


        foreach (float k in keys.filterArray(IaA.Keys.ToArray()))
        {
            Interval alphasForThetaA, alphasForThetaB;

            IaA.TryGetValue(k, out alphasForThetaA);
            float kAdjust = mergeKey(k, IaB.Keys.ToArray(), samplingRate);
            IaB.TryGetValue(kAdjust, out alphasForThetaB);

            Interval alphaDouble = alphasForThetaA.Intersect(alphasForThetaB);

            if (alphaDouble != null) Ia.Add(k, alphaDouble);
        }
        return Ia;
    }

  

    /// <summary>
    /// The alpha interval from the vantage angle constraint for a given beta and theta
    /// </summary>
    /// <param name="theta">theta</param>
    /// <param name="beta">beta prime</param>
    /// <returns>an alpha interval</returns>
    public Interval GetVantageAlphaInterval(float theta, Interval beta)
    {
       
        float alphaMin = Mathf.PI - theta / 2 - (beta.LOWERBOUND);
        float alphaMax = Mathf.PI - theta / 2 - (beta.UPPERBOUND);

        return new Interval(alphaMin, alphaMax);
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

    public Interval getAlphaIntervalFromOnscreenPositions(Vector2 desPosA, Vector2 desPosB,Boolean visualize = false)
    {

        
        List<Vector2> verticesFrameTargetA = getFrameTarget(desPosA);
        List<Vector2> verticesFrameTargetB = getFrameTarget(desPosB);



        List<Vector3> VerticesInCamSpaceA = new List<Vector3>();
        List<Vector3> VerticesInCamSpaceB = new List<Vector3>();
        foreach (Vector2 vec in verticesFrameTargetA)
        {
            Vector3 inCamSpace = ToricComputing.GetVectorInCameraSpace(vec);
            VerticesInCamSpaceA.Add(inCamSpace);

        }
        foreach (Vector2 vec in verticesFrameTargetB)
        {
            Vector3 inCamSpace = ToricComputing.GetVectorInCameraSpace(vec);
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


        if (visualize)
        {
            _visualizeOnscreenPosition(verticesFrameTargetA, verticesFrameTargetB);
        }


        return new Interval(Mathf.Min(alphaValues.ToArray()), Mathf.Max(alphaValues.ToArray()));
    }

    






    //private methods

        //TODO adjust to vecAB instead of Vector3.up
    private void _visualizeDistance(Vector2 distanceToA, Vector2 distanceToB)
    {
        Ellipse minDistanceA = new Ellipse(distanceToA[0], distanceToA[0], _A, Vector3.up);
        minDistanceA.Draw(Color.red);

        Ellipse maxDistanceA = new Ellipse(distanceToA[1], distanceToA[1], _A, Vector3.up);
        maxDistanceA.Draw(Color.green);

        Ellipse minDistanceB = new Ellipse(distanceToB[0], distanceToB[0], _B, Vector3.up);
        minDistanceB.Draw(Color.red);

        Ellipse maxDistanceB = new Ellipse(distanceToB[1], distanceToB[1], _B, Vector3.up);
        maxDistanceB.Draw(Color.green);
    }


    //TODO
    private void _visualizeOnscreenPosition(List<Vector2> verticesFrameTargetA, List<Vector2> verticesFrameTargetB)
    {

        Vector3 previousPoint = _A;

        foreach (Vector2 item in verticesFrameTargetA)
        {
            Vector3 cornerpoint = new Vector3(_A.x + item.x, _A.y + item.y, _A.z);
            Debug.DrawLine(previousPoint, cornerpoint, Color.red, Mathf.Infinity);
            previousPoint = cornerpoint;
        }


         previousPoint = _B;

        foreach (Vector2 item in verticesFrameTargetB)
        {
            Vector3 cornerpoint = new Vector3(_B.x + item.x, _B.y + item.y, _B.z);
            Debug.DrawLine(previousPoint, cornerpoint, Color.blue, Mathf.Infinity);
            previousPoint = cornerpoint;
        }
    }


















    private Dictionary<float, Interval> getIntervalFromB(float minDistanceToB, float maxDistanceToB, float samplingRate = 0.05f)
    {

        Dictionary<float, Interval> invDistanceMin = GetAlphaFromDistanceToB(minDistanceToB, samplingRate);
        Dictionary<float, Interval> invDistanceMax = GetAlphaFromDistanceToB(maxDistanceToB, samplingRate);
        Dictionary<float, Interval> invDistanceAB = GetAlphaFromDistanceToB(_AB.magnitude + 0.01f, samplingRate);

        Dictionary<float, Interval> thetaIntersect = new Dictionary<float, Interval>();

        float[] keysMin, KeysMax;
        keysMin = invDistanceMin.Keys.ToArray();
        KeysMax = invDistanceMax.Keys.ToArray();

        if (minDistanceToB <= _AB.magnitude && maxDistanceToB > _AB.magnitude) keysMin = KeysMax;

        foreach (float theta in KeysMax.Intersect(keysMin))
        {
            Interval alphaFromMin, alphaFromMax;
            invDistanceMin.TryGetValue(theta, out alphaFromMin);
            invDistanceMax.TryGetValue(theta, out alphaFromMax);

            if (alphaFromMin == null) invDistanceAB.TryGetValue(theta, out alphaFromMin);

            Interval alphaInterval = new Interval(alphaFromMin.LOWERBOUND, alphaFromMax.LOWERBOUND);

            thetaIntersect.Add(theta, alphaInterval);
        }

        return thetaIntersect;
    }

    
    private Dictionary<float, Interval> getIntervalFromA(float minDistanceToA, float maxDistanceToA, float samplingRate = 0.05f)
    {
        Interval possibleThetaValues = new Interval(0.01f, Mathf.PI * 2, samplingRate);
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

    
    private Dictionary<float, Interval> GetAlphaFromDistanceToB(float distance, float samplingRate)
    {
        Dictionary<float, Interval> res = new Dictionary<float, Interval>();

        Interval thetaValues = new Interval(samplingRate, Mathf.PI * 2 - samplingRate, samplingRate);

        Boolean smallerThanAB = distance <= _AB.magnitude;

        if (smallerThanAB)
        {
            thetaValues = new Interval(samplingRate, 2 * Mathf.Asin(distance / _AB.magnitude), samplingRate);
        }

        foreach (float theta in thetaValues.getEveryValue())
        {
            float[] alphaRes = new float[2];

            float acos = Mathf.Acos(_AB.magnitude / distance * Mathf.Sin(theta / 2));

            alphaRes[0] = (Mathf.PI / 2 - acos);

            alphaRes[1] = (Mathf.PI / 2 + acos);

            if (smallerThanAB) alphaRes[1] = Mathf.PI - theta / 2;

            res.Add(theta, Interval.fromFloatArray(alphaRes));
        }
        return res;
    }

   
    private float GetAlphaFromDistanceToA(float distance, float theta)
    {
        return Mathf.Acos((distance - _AB.magnitude * Mathf.Cos(theta / 2)) /
            (Mathf.Sqrt((distance * distance) + (_AB.magnitude * _AB.magnitude)
            - 2 * _AB.magnitude * distance * Mathf.Cos(theta / 2))));
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

    //creates a list of vertices belonging to an octagon around a target 
    //TODO adjust _frameradius to screensize, instead of world space
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

}

